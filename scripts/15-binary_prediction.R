# ====================== CLIM_RISK_PHYLOSDM | BINARY MAP GENERATION ======================
# Creates binary presence/absence maps from spatial predictions using computed thresholds
# ========================================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(terra)
})

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
    message("Running on HPC")
    Sys.setenv(TMPDIR = "/vast/palmer/scratch/jetz/ss4224")
    terra::terraOptions(
        tempdir = "/vast/palmer/scratch/jetz/ss4224",
        memfrac = 0.1,
        todisk = TRUE
    )
} else {
    root <- "~/phyloSDM_MOL"
    message("Running locally")
}

# ---- Command line arguments or manual set ----
if (interactive()) {
    EXP_ROOT <- "v0"
    EXP_ID <- "sub1000"
    DATASET <- "amphibians"
    CLUSTER <- "Kass1"
    FSP <- "ALL"
    REPNO <- 1
    NREP <- 1
    MODEL_TYPE <- "STAN"
    MODEL_NAME <- "LGCP_background"
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
message(sprintf("Binary map generation for cluster: %s, rep: %d", CLUSTER, REPNO))

# ---- Configuration ----
CROP_TO_EXTENT <- TRUE
EXTENT_BUFFER_DEG <- 4 # Buffer around expert range extent (degrees)

# ---- Paths ----
raw_dir <- file.path(root, "raw_data")
analysis_dir <- file.path(root, "analysis")
dpath <- file.path(raw_dir, DATASET)
data_dir <- file.path(dpath, CLUSTER)
expert_range_dir <- file.path(root, "expert_ranges")

# Input/output directories
pred_dir <- file.path(analysis_dir, EXP_ROOT, "spatial_pred")
eval_dir <- file.path(analysis_dir, EXP_ROOT, "eval")
threshold_dir <- file.path(eval_dir, "thresholds")
binary_dir <- file.path(analysis_dir, EXP_ROOT, "binary_maps")

dir.create(binary_dir, recursive = TRUE, showWarnings = FALSE)

# ---- File naming ----
base_prefix <- paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO)

# ---- Load threshold data ----
threshold_file <- file.path(threshold_dir, paste0(base_prefix, "_thresholds.csv"))
if (!file.exists(threshold_file)) {
    stop(sprintf("Threshold file not found: %s", threshold_file))
}
threshold_df <- read.csv(threshold_file)
message(sprintf("Loaded thresholds for %d species", nrow(threshold_df)))

# ---- Load model data to get species list ----
model_data_file <- file.path(data_dir, paste0(base_prefix, "_model_data.Rdata"))
if (!file.exists(model_data_file)) {
    stop(sprintf("Model data file not found: %s", model_data_file))
}
load(model_data_file) # loads 'model_data'
species_names <- model_data$species_names
J <- length(species_names)
message(sprintf("Processing %d species", J))

# ---- Initialize tracking vectors ----
success_species <- character()
missing_pred <- character()
missing_threshold <- character()
error_species <- character()

# ---- Main loop over species ----
message(sprintf("\n========== Generating binary maps for %d species ==========", J))

for (j in seq_len(J)) {
    sp_name <- species_names[j]
    message(sprintf("\n[%d/%d] Processing species: %s", j, J, sp_name))

    tryCatch(
        {
            # ---- Check for prediction raster ----
            pred_file <- file.path(pred_dir, paste0(base_prefix, "_", sp_name, "_relprob.tif"))
            if (!file.exists(pred_file)) {
                message(sprintf("  Prediction raster not found: %s", basename(pred_file)))
                missing_pred <- c(missing_pred, sp_name)
                next
            }

            # ---- Get threshold for this species ----
            sp_threshold_row <- which(threshold_df$species_name == sp_name)
            if (length(sp_threshold_row) == 0) {
                message(sprintf("  No threshold found for species: %s", sp_name))
                missing_threshold <- c(missing_threshold, sp_name)
                next
            }

            thr <- threshold_df$threshold[sp_threshold_row]
            auc <- threshold_df$AUC[sp_threshold_row]

            if (is.na(thr)) {
                message("  Threshold is NA, using 0")
                thr <- 0
            }

            message(sprintf("  Threshold: %.4f (AUC: %.3f)", thr, ifelse(is.na(auc), NA, auc)))

            # ---- Load prediction raster ----
            pred_rast <- terra::rast(pred_file)

            # ---- Create binary map ----
            binary_rast <- terra::ifel(pred_rast >= thr, 1, 0)

            # ---- Optionally crop to expert range extent ----
            if (CROP_TO_EXTENT) {
                # Try to find expert range shapefile
                range_file <- file.path(expert_range_dir, paste0(sp_name, ".gpkg"))
                if (!file.exists(range_file)) {
                    range_file <- file.path(expert_range_dir, sp_name, paste0(sp_name, ".shp"))
                }
                if (!file.exists(range_file)) {
                    range_file <- file.path(expert_range_dir, paste0(sp_name, ".shp"))
                }

                if (file.exists(range_file)) {
                    range_vect <- terra::vect(range_file)
                    range_ext <- terra::ext(range_vect)

                    # Add buffer to extent
                    buffered_ext <- terra::ext(
                        range_ext$xmin - EXTENT_BUFFER_DEG,
                        range_ext$xmax + EXTENT_BUFFER_DEG,
                        range_ext$ymin - EXTENT_BUFFER_DEG,
                        range_ext$ymax + EXTENT_BUFFER_DEG
                    )

                    binary_rast <- terra::crop(binary_rast, buffered_ext)
                    message("  Cropped to expert range extent (+buffer)")
                } else {
                    message("  No expert range found, keeping full extent")
                }
            }

            # ---- Write binary raster ----
            out_file <- file.path(binary_dir, paste0(base_prefix, "_", sp_name, "_binary.tif"))
            terra::writeRaster(binary_rast, filename = out_file, overwrite = TRUE)
            message(sprintf("  Wrote: %s", basename(out_file)))

            # ---- Print summary stats ----
            n_cells <- terra::global(binary_rast, "notNA")[[1]]
            n_presence <- terra::global(binary_rast, "sum", na.rm = TRUE)[[1]]
            pct_presence <- 100 * n_presence / n_cells
            message(sprintf(
                "  Cells: %d total, %d presence (%.1f%%)",
                n_cells, n_presence, pct_presence
            ))

            success_species <- c(success_species, sp_name)

            # Clean up
            rm(pred_rast, binary_rast)
            gc()
        },
        error = function(e) {
            message(sprintf("  ERROR: %s", e$message))
            error_species <- c(error_species, sp_name)
        }
    )
}

# ---- Summary report ----
message("\n========== Binary Map Generation Summary ==========")
message(sprintf("Successfully processed: %d species", length(success_species)))
message(sprintf("Missing prediction rasters: %d species", length(missing_pred)))
message(sprintf("Missing thresholds: %d species", length(missing_threshold)))
message(sprintf("Errors: %d species", length(error_species)))

if (length(missing_pred) > 0) {
    message("\nSpecies with missing predictions:")
    message(paste("  ", missing_pred, collapse = "\n"))
}

if (length(missing_threshold) > 0) {
    message("\nSpecies with missing thresholds:")
    message(paste("  ", missing_threshold, collapse = "\n"))
}

if (length(error_species) > 0) {
    message("\nSpecies with errors:")
    message(paste("  ", error_species, collapse = "\n"))
}

# ---- Save processing summary ----
summary_df <- data.frame(
    species_name = species_names,
    status = ifelse(species_names %in% success_species, "success",
        ifelse(species_names %in% missing_pred, "missing_prediction",
            ifelse(species_names %in% missing_threshold, "missing_threshold",
                ifelse(species_names %in% error_species, "error", "unknown")
            )
        )
    )
)

summary_file <- file.path(binary_dir, paste0(base_prefix, "_binary_summary.csv"))
write.csv(summary_df, summary_file, row.names = FALSE)
message(sprintf("\nSaved processing summary: %s", summary_file))

message(sprintf("\n========== Done! Binary maps saved to: %s ==========", binary_dir))
