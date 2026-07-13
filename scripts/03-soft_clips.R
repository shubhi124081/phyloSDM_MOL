# Purpose: Generate soft clip rasters for amphibian species
# Output: Soft clip rasters saved to analysis/soft_clips/

rm(list = ls())

suppressPackageStartupMessages({
  library(terra)
})

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
  root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
  epath <- "/vast/palmer/pi/jetz/ss4224/env"
  message("Running on HPC")
  Sys.setenv(TMPDIR = "/vast/palmer/scratch/jetz/ss4224")
  terra::terraOptions(tempdir = "/vast/palmer/scratch/jetz/ss4224")
} else {
  root <- "~/phyloSDM_MOL"
  epath <- "~/env"
  message("Running locally")
}

# ---- Args ----
if (interactive()) {
  EXP_ROOT <- "v0"
  EXP_ID <- "test_DSQ"
  DATASET <- "amphibians"
  CLUSTER <- "Chal1"
  FSP <- "ALL"
  REPNO <- 1
  NREP <- 1
  MODEL_TYPE <- "TMB"
  MODEL_NAME <- "05_lgcp_correceted.cpp"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  EXP_ROOT <- args[1]
  EXP_ID <- args[2]
  DATASET <- args[3]
  CLUSTER <- args[4]
  FSP <- args[5]
  REPNO <- as.integer(args[6])
  NREP <- as.integer(args[7])
  MODEL_TYPE <- args[8]
  MODEL_NAME <- args[9]
}
message(sprintf("Processing cluster: %s", CLUSTER))

dpath <- file.path(root, "raw_data", DATASET)
expert_directory <- file.path(root, "expert_ranges")
analysis_directory <- file.path(root, "analysis")
output_dir <- file.path(analysis_directory, "soft_clips")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load cluster species list ----
load(file.path(dpath, "spList.Rdata")) # loads spList
if (!(CLUSTER %in% names(spList))) {
  stop(sprintf(
    "Cluster '%s' not found in spList. Available: %s",
    CLUSTER, paste(head(names(spList), 10), collapse = ", ")
  ))
}
sps <- spList[[CLUSTER]]
message(sprintf("Cluster has %d species", length(sps)))

# ---- Template raster / extent ----
template_file <- file.path(epath, "CHELSA_bio_1.tif")
if (!file.exists(template_file)) stop("Template raster not found: ", template_file)

extent_file <- file.path(dpath, CLUSTER, paste0(CLUSTER, "_extent.csv"))
if (!file.exists(extent_file)) stop("Extent file not found. Run 03-gen_data.R first for this cluster.")

ex_df <- read.csv(extent_file)
ex <- terra::ext(ex_df$xmin, ex_df$xmax, ex_df$ymin, ex_df$ymax)
message(sprintf(
  "Loaded cluster extent: [%.2f, %.2f] x [%.2f, %.2f]",
  ex_df$xmin, ex_df$xmax, ex_df$ymin, ex_df$ymax
))

template_rast <- terra::rast(template_file)
r <- terra::crop(template_rast, ex)
message(sprintf("Template raster: %d x %d cells", terra::nrow(r), terra::ncol(r)))

# ---- Compute low-res raster ONCE ----
fact <- 10
message(sprintf("Building low-res raster once (aggregate fact=%d)...", fact))
r_lowres <- terra::aggregate(r, fact = fact)
message(sprintf("Low-res raster: %d x %d cells", terra::nrow(r_lowres), terra::ncol(r_lowres)))

# ---- Logistic params ----
X0 <- 30000
L <- 1
Kp <- 0.0002

# ---- Process each species ----
for (i in seq_along(sps)) {
  sp <- sps[i]

  # Check if soft clip already exists
  out_file <- file.path(output_dir, paste0(sp, "_soft_clip.tif"))
  if (file.exists(out_file)) {
    message(sprintf("[%d/%d] Skipping %s (soft clip already exists)", i, length(sps), sp))
    next
  }

  message(sprintf("[%d/%d] Soft clipping: %s", i, length(sps), sp))

  tryCatch(
    {
      gpkg_file <- file.path(expert_directory, paste0(sp, ".gpkg"))
      if (!file.exists(gpkg_file)) {
        warning(sprintf("Expert range not found for %s, skipping", sp))
        next
      }

      range <- terra::vect(gpkg_file)
      range_ext <- terra::ext(range)
      # This is a little expensive
      range_buffer <- terra::buffer(range, width = 1000000) # 1000 km buffer
      # range_buffer_ext <- terra::ext(range_buffer)

      # THis is cheaper
      buf_km <- 1000
      buf_deg_lat <- buf_km / 111.0 # ~ degrees latitude
      cc <- terra::crds(terra::centroids(range), df = TRUE)[1, ]
      lat <- cc$y
      buf_deg_lon <- buf_km / (111.0 * cos(lat * pi / 180))
      range_ext <- terra::ext(range)

      range_buffer_ext <- terra::ext(
        range_ext$xmin - buf_deg_lon, range_ext$xmax + buf_deg_lon,
        range_ext$ymin - buf_deg_lat, range_ext$ymax + buf_deg_lat
      )

      # clamp to cluster extent
      range_buffer_ext <- terra::intersect(range_buffer_ext, ex)

      r_lowres_crop <- terra::crop(r_lowres, range_buffer_ext, snap = "out")



      # Ensure buffered extent doesn't exceed cluster extent
      range_buffer_ext <- terra::intersect(range_buffer_ext, ex)

      r_lowres_crop <- terra::crop(r_lowres, range_buffer_ext)

      # Distance on low-res raster (still the heavyweight, but now only once per species)
      distance_lowres <- terra::distance(r_lowres_crop, range)

      # ---- IMPORTANT CHANGE 2: do logistic as raster algebra (no values()/setValues()) ----
      # logistic = L / (1 + exp(-Kp * (d - X0)))
      # soft mask = 1 - logistic
      decay_lowres_raster <- 1 - (L / (1 + exp(-Kp * (distance_lowres - X0))))

      # Resample to high-res prediction raster
      decay_highres <- terra::resample(decay_lowres_raster, r, method = "bilinear")

      # Rescale so that max is exactly 1
      max_val <- terra::global(decay_highres, "max", na.rm = TRUE)[1, 1]
      if (!is.na(max_val) && max_val > 0) {
        decay_highres <- decay_highres / max_val
      }

      terra::writeRaster(decay_highres, out_file, overwrite = TRUE)
      message(sprintf("  Saved: %s", basename(out_file)))
    },
    error = function(e) {
      warning(sprintf("Error processing %s: %s", sp, e$message))
    },
    finally = {
      rm(list = intersect(
        c("range", "distance_lowres", "decay_lowres_raster", "decay_highres", "max_val"),
        ls()
      ))
      gc()
    }
  )
}

message(sprintf("Done! Soft clips saved to: %s", output_dir))
