# Purpose: Generate training data for a single amphibian cluster
# Input: cluster_name passed as command line argument or set manually
# Output: Processed data files saved to raw_data/amphibians/<cluster>/

# Clear workspace
rm(list = ls())

# Libraries
library(ggplot2)

# ---- Command line arguments or manual set ----
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

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    print("Running on HPC")
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
    epath <- "/vast/palmer/pi/jetz/ss4224/env"
    species_csv_path <- "/vast/palmer/pi/jetz/newt_scratch/standard_occurrences/v1.2/harmonized/amphibians"
    terra::terraOptions(tempdir = "/vast/palmer/scratch/jetz/ss4224")
} else {
    print("Running locally")
    root <- "~/phyloSDM_MOL"
    epath <- "~/env"
    species_csv_path <- file.path(root, "raw_data", DATASET, "species_csvs")
}

# Paths that depend on root/DATASET
dpath <- file.path(root, "raw_data", DATASET)
tree_path <- file.path(dpath, "02-phylogeny_pruned_to_working_species_range_x_tree.Rdata")
expert_range_path <- file.path(root, "expert_ranges")

# ---- Configuration ----
env_files <- c(
    "CHELSA_bio_1.tif", # mean annual temp
    "CHELSA_bio_4.tif", # temp seasonality
    "CHELSA_bio_13.tif", # precip of wettest month
    "CHELSA_bio_15.tif", # precip seasonailty
    "cloudCover.tif", # cloud cover
    "Annual_EVI.tif", # annual evi
    "TRI.tif", # topographic ruggedness index
    "elevation_1KMmean_SRTM.tif" # elevation
)
env_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# ---- User-specified configuration ----
# Extent widening parameter (degrees) - for cluster extent only
WIDENBY <- 5

# Buffer distance (km) around presence points to mask from BG sampling
BUFFER_KM <- 20

# Buffer to add to species-specific extent (degrees, ~300km)
EXTENT_BUFFER_DEG <- 5

# ---- Load cluster species list ----
load(file.path(dpath, "spList.Rdata")) # loads spList
if (!(CLUSTER %in% names(spList))) {
    stop(sprintf(
        "Cluster '%s' not found in spList. Available: %s",
        CLUSTER, paste(head(names(spList), 10), collapse = ", ")
    ))
}

# ---- Determine species to process ----
if (FSP == "ALL") {
    sps <- spList[[CLUSTER]]
} else {
    stop("This is a cluster-level script; FSP other than ALL not supported")
}

message(sprintf("Cluster has %d species", length(sps)))

# ---- Load tree ----
contents <- load(tree_path)
# This loads: tree_working
tree_full <- tree_working

# ---- Load species occurrence data ----
# Each species has its own CSV file
message("Loading species CSVs...")
raw_list <- list()
missing_species <- c()

for (sp in sps) {
    csv_file <- file.path(species_csv_path, paste0(sp, ".csv"))
    if (file.exists(csv_file)) {
        sp_data <- read.csv(csv_file)
        # Standardize column names
        colnames(sp_data) <- tolower(colnames(sp_data))
        # Ensure lat/lon columns exist
        if ("latitude" %in% colnames(sp_data)) {
            colnames(sp_data)[colnames(sp_data) == "latitude"] <- "lat"
        }
        if ("longitude" %in% colnames(sp_data)) {
            colnames(sp_data)[colnames(sp_data) == "longitude"] <- "lon"
        }
        if (!all(c("lat", "lon") %in% colnames(sp_data))) {
            warning(sprintf("Species %s CSV missing lat/lon columns, skipping", sp))
            missing_species <- c(missing_species, sp)
            next
        }
        sp_data <- sp_data[complete.cases(sp_data[, c("lat", "lon")]), ]
        if (nrow(sp_data) > 0) {
            sp_data$species <- sp
            raw_list[[sp]] <- sp_data[, c("lat", "lon", "species")]
        } else {
            missing_species <- c(missing_species, sp)
        }
    } else {
        message(sprintf("CSV not found for species: %s", sp))
        missing_species <- c(missing_species, sp)
    }
}

# Update species list removing missing
if (length(missing_species) > 0) {
    message(sprintf(
        "Removing %d species with missing/empty CSVs: %s",
        length(missing_species),
        paste(missing_species, collapse = ", ")
    ))
    sps <- setdiff(sps, missing_species)
}

if (length(sps) == 0) {
    stop("No species with valid occurrence data found!")
}
message(sprintf("Processing %d species with data", length(sps)))

# Combine into single data frame
raw <- do.call(rbind, raw_list)
raw <- as.data.frame(raw)

# ---- Helper functions ----
getExtentDf <- function(cood, NAMES = c("lat", "lon")) {
    lat <- cood[[NAMES[1]]]
    lon <- cood[[NAMES[2]]]
    c(
        xmin = min(lon, na.rm = TRUE),
        xmax = max(lon, na.rm = TRUE),
        ymin = min(lat, na.rm = TRUE),
        ymax = max(lat, na.rm = TRUE)
    )
}

cropDFCood <- function(df, ext) {
    df[df$lon >= ext[1] & df$lon <= ext[2] &
        df$lat >= ext[3] & df$lat <= ext[4], , drop = FALSE]
}

annotateCoods <- function(envpath, files, crs, COOD) {
    result <- matrix(NA_real_, nrow = nrow(COOD), ncol = length(files))
    colnames(result) <- gsub("\\.tif.*", "", files)
    pts <- terra::vect(COOD, geom = c("lon", "lat"), crs = crs)
    for (i in seq_along(files)) {
        rast_file <- file.path(envpath, files[i])
        if (file.exists(rast_file)) {
            r <- terra::rast(rast_file)
            result[, i] <- terra::extract(r, pts)[, 2]
        } else {
            warning(sprintf("Env file not found: %s", rast_file))
        }
    }
    as.data.frame(result)
}

process_x_matrix <- function(x) {
    x <- as.data.frame(x)
    x$Intercept <- 1
    x
}

# ---- Prune tree to available species ----
tree_tips <- tree_full$tip.label
sps_in_tree <- intersect(sps, tree_tips)

if (length(sps_in_tree) < length(sps)) {
    not_in_tree <- setdiff(sps, sps_in_tree)
    message(sprintf(
        "Removing %d species not in tree: %s",
        length(not_in_tree),
        paste(not_in_tree, collapse = ", ")
    ))
    sps <- sps_in_tree
    raw <- raw[raw$species %in% sps, ]
}

if (length(sps) < 2) {
    stop("Need at least 2 species in tree to continue")
}

tree <- ape::keep.tip(tree_full, tip = sps)
message(sprintf("Tree pruned to %d tips", length(tree$tip.label)))

# ---- Compute cluster extent ----
ex <- getExtentDf(raw, NAMES = c("lat", "lon"))
ex <- terra::ext(ex + (c(-1, 1, -1, 1) * WIDENBY)) # widen
message(sprintf(
    "Extent: xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f",
    ex[1], ex[2], ex[3], ex[4]
))

# ---- Create empty raster at extent ----
# Use first env file as template for resolution
template_file <- file.path(epath, env_files[1])
if (file.exists(template_file)) {
    template_rast <- terra::rast(template_file)
    template_res <- terra::res(template_rast)
} else {
    # Default to ~1km resolution
    template_res <- c(1 / 120, 1 / 120)
}

# ---- Rasterize occurrences per species ----
# Use species-specific extents to avoid memory issues with global-scale clusters
message("Rasterizing occurrences (species-specific extents)...")
data <- vector("list", length(sps))
NAMES <- c("lon", "lat")

for (i in seq_along(sps)) {
    sp <- sps[i]
    sp_raw <- raw[raw$species == sp, ]

    if (nrow(sp_raw) == 0) {
        data[[i]] <- data.frame()
        next
    }

    # Compute species-specific extent with buffer
    sp_ext <- c(
        xmin = min(sp_raw$lon, na.rm = TRUE) - WIDENBY,
        xmax = max(sp_raw$lon, na.rm = TRUE) + WIDENBY,
        ymin = min(sp_raw$lat, na.rm = TRUE) - WIDENBY,
        ymax = max(sp_raw$lat, na.rm = TRUE) + WIDENBY
    )
    sp_terra_ext <- terra::ext(sp_ext)

    # Create species-specific template raster
    er_sp <- terra::rast(sp_terra_ext, res = template_res, crs = env_crs)

    # Create spatial vector
    sp_vect <- terra::vect(sp_raw[, NAMES, drop = FALSE], geom = NAMES, crs = env_crs)

    # Rasterize, summing occurrences per cell
    er1 <- terra::rasterize(sp_vect, er_sp, fun = "sum")
    er1_df <- terra::as.data.frame(er1, xy = TRUE)

    if (nrow(er1_df) == 0) {
        data[[i]] <- data.frame()
        next
    }

    colnames(er1_df) <- c("x", "y", "sum")
    er1_df$sp <- sp

    # Subsample if species has more than 1000 sites
    MAX_SITES_PER_SPECIES <- 1000
    if (nrow(er1_df) > MAX_SITES_PER_SPECIES) {
        message(sprintf("  Species %s has %d sites, subsampling to %d", sp, nrow(er1_df), MAX_SITES_PER_SPECIES))
        set.seed(42 + i) # Reproducible subsampling per species
        subsample_idx <- sample(nrow(er1_df), MAX_SITES_PER_SPECIES, replace = FALSE)
        er1_df <- er1_df[subsample_idx, , drop = FALSE]
    }

    data[[i]] <- er1_df
}

names(data) <- sps

# Remove species with no rasterized data
empty_sps <- sps[vapply(data, nrow, 1L) == 0L]
if (length(empty_sps) > 0) {
    message(sprintf("Removing %d species with no rasterized data", length(empty_sps)))
    data <- data[!(names(data) %in% empty_sps)]
    sps <- setdiff(sps, empty_sps)
    tree <- ape::drop.tip(tree, empty_sps[empty_sps %in% tree$tip.label])
}

if (length(sps) == 0) {
    stop("No species remaining after rasterization!")
}

# ---- Wide site x species counts ----
message("Creating site x species matrix...")
data_df <- do.call(rbind, data)
data_df <- as.data.frame(data_df)

# Pivot to wide format using base R reshape
data_wide <- stats::reshape(
    data_df[, c("x", "y", "sp", "sum")],
    idvar = c("x", "y"),
    timevar = "sp",
    direction = "wide",
    v.names = "sum"
)
# Clean column names (remove "sum." prefix)
colnames(data_wide) <- gsub("^sum\\.", "", colnames(data_wide))
# Replace NAs with 0
data_wide[is.na(data_wide)] <- 0
data_sps <- data_wide

# ---- Extract y and cood ----
# Ensure sps only includes species present in data_sps
sps <- intersect(sps, colnames(data_sps))
y <- data_sps[, sps, drop = FALSE]
cood <- data.frame(lat = data_sps$y, lon = data_sps$x)
cood <- cood[, c("lat", "lon")]

message(sprintf("Sites: %d, Species: %d", nrow(y), length(sps)))

# After creating y but before BG generation
for (i in seq_along(sps)) {
    n_pres <- sum(y[, i] > 0)
    if (n_pres < 2) {
        message(sprintf("Species %s has only %d presence sites; these will be dropped", sps[i], n_pres))
    }
}

# ---- Annotate with environmental covariates ----
message("Annotating with environmental covariates...")
x <- annotateCoods(epath, env_files, env_crs, cood)

# ---- Drop incomplete rows ----
rmv <- which(!complete.cases(x))
if (length(rmv) > 0) {
    message(sprintf("Removing %d rows with missing env data", length(rmv)))
    x <- x[-rmv, , drop = FALSE]
    y <- y[-rmv, , drop = FALSE]
    cood <- cood[-rmv, , drop = FALSE]
}

# Handle non-finite values
x <- as.matrix(x)
nonfinite <- which(!is.finite(x))
if (length(nonfinite)) {
    message(sprintf("Replaced %d non-finite values with small eps", length(nonfinite)))
    x[nonfinite] <- 1e-6
}
x <- process_x_matrix(as.data.frame(x))

# ---- Synchronize rownames ----
rownames(x) <- rownames(cood) <- rownames(y) <- seq_len(nrow(y))

# ---- Background points per species ----
message("Generating background points...")
J <- length(sps)
y$bg <- 0L
y <- as.data.frame(y)

# Initialize as empty lists to collect BG data
x_bg_list <- list()
y_bg_list <- list()
cood_bg_list <- list()

# Load world land mask (sf object)
world_file <- file.path(root, "raw_data", "world.Rdata")
if (file.exists(world_file)) {
    load(world_file) # loads 'world' as sf object
    world_vect <- terra::vect(world)
    message("Loaded world land mask")
} else {
    warning("World land mask not found, ocean masking disabled")
    world_vect <- NULL
}

set.seed(42)

for (species_idx in seq_len(J)) {
    sp_name <- sps[species_idx]
    message(sprintf("BG for species %d/%d: %s", species_idx, J, sp_name))

    tryCatch(
        {
            pres_idx <- which(y[, species_idx] > 0)
            n_pres <- length(pres_idx)

            # Get presence coordinates for this species (may be empty or single point)
            pres_cood <- cood[pres_idx, , drop = FALSE]

            # Try to load expert range map
            gpkg_file <- file.path(expert_range_path, paste0(sp_name, ".gpkg"))
            has_expert_range <- file.exists(gpkg_file)

            # If no presence points and no expert range, skip
            if (n_pres == 0 && !has_expert_range) {
                warning("No presence points and no expert range, skipping BG")
                next
            }

            # Compute presence extent (only if we have presences)
            if (n_pres > 0) {
                pres_ext <- c(
                    xmin = min(pres_cood$lon, na.rm = TRUE),
                    xmax = max(pres_cood$lon, na.rm = TRUE),
                    ymin = min(pres_cood$lat, na.rm = TRUE),
                    ymax = max(pres_cood$lat, na.rm = TRUE)
                )
            } else {
                pres_ext <- NULL
            }

            # Load expert range if available
            if (has_expert_range) {
                expert_range <- terra::vect(gpkg_file)
                expert_ext <- terra::ext(expert_range)

                if (!is.null(pres_ext)) {
                    # Union of expert range and presence extents
                    combined_ext <- c(
                        xmin = min(pres_ext["xmin"], expert_ext[1]),
                        xmax = max(pres_ext["xmax"], expert_ext[2]),
                        ymin = min(pres_ext["ymin"], expert_ext[3]),
                        ymax = max(pres_ext["ymax"], expert_ext[4])
                    )
                    message(sprintf("  Loaded expert range, combined extent"))
                } else {
                    # Use expert range extent only
                    combined_ext <- c(
                        xmin = expert_ext[1],
                        xmax = expert_ext[2],
                        ymin = expert_ext[3],
                        ymax = expert_ext[4]
                    )
                    message(sprintf("  Using expert range extent only"))
                }
            } else {
                # No expert range - use presence extent only
                if (n_pres < 2) {
                    warning("Only 1 presence point and no expert range, skipping BG")
                    next
                }
                combined_ext <- pres_ext
                message(sprintf("  No expert range found, using presence extent only"))
            }

            # Add 300km buffer to extent (use unname to avoid named vector issues)
            sp_ext_wide <- c(
                xmin = unname(combined_ext["xmin"]) - EXTENT_BUFFER_DEG,
                xmax = unname(combined_ext["xmax"]) + EXTENT_BUFFER_DEG,
                ymin = unname(combined_ext["ymin"]) - EXTENT_BUFFER_DEG,
                ymax = unname(combined_ext["ymax"]) + EXTENT_BUFFER_DEG
            )

            message(sprintf(
                "  Sampling extent: lon [%.2f, %.2f], lat [%.2f, %.2f]",
                sp_ext_wide["xmin"], sp_ext_wide["xmax"],
                sp_ext_wide["ymin"], sp_ext_wide["ymax"]
            ))

            # Number of BG points = number of presence points (minimum 1)
            n_bg_per_species <- max(1, n_pres)

            # Create mask raster at ~5km resolution for efficiency
            mask_ext <- terra::ext(sp_ext_wide)
            mask_rast <- terra::rast(mask_ext, res = 0.05, crs = env_crs)
            terra::values(mask_rast) <- 1L # Initialize all cells as valid

            # Mask out ocean/sea (keep only land)
            if (!is.null(world_vect)) {
                world_crop <- terra::crop(world_vect, mask_ext)
                if (!is.null(world_crop) && length(world_crop) > 0) {
                    land_rast <- terra::rasterize(world_crop, mask_rast, field = 1, background = NA)
                    mask_rast <- terra::mask(mask_rast, land_rast)
                }
            }

            # Mask out buffer around presence points
            if (nrow(pres_cood) > 0) {
                pres_vect <- terra::vect(pres_cood, geom = c("lon", "lat"), crs = env_crs)
                # Buffer in meters, convert km to m
                pres_buffer <- terra::buffer(pres_vect, width = BUFFER_KM * 2000)
                buffer_rast <- terra::rasterize(pres_buffer, mask_rast, field = 1, background = NA)
                # Set buffered areas to NA
                mask_rast[!is.na(buffer_rast)] <- NA
            }

            # Get valid cell coordinates
            valid_cells <- terra::as.data.frame(mask_rast, xy = TRUE, na.rm = TRUE)
            message(sprintf("  Valid cells after masking: %d", nrow(valid_cells)))

            if (nrow(valid_cells) == 0) {
                warning("No valid cells after land/buffer masking")
                next
            }

            # Oversample only to compensate for missing env data (2x should suffice)
            OVERSAMPLE_FACTOR <- 2
            n_candidates <- min(nrow(valid_cells), n_bg_per_species * OVERSAMPLE_FACTOR)

            # Sample uniformly from valid cells
            if (nrow(valid_cells) <= n_candidates) {
                sampled_idx <- seq_len(nrow(valid_cells))
            } else {
                sampled_idx <- sample(nrow(valid_cells), n_candidates)
            }
            bg_cood_tmp <- data.frame(
                lat = valid_cells$y[sampled_idx],
                lon = valid_cells$x[sampled_idx]
            )
            message(sprintf("  Sampled %d candidate points", nrow(bg_cood_tmp)))

            # Annotate with env vars (single call - most expensive operation)
            x_bg1 <- annotateCoods(epath, env_files, env_crs, bg_cood_tmp)

            # Keep only complete cases
            complete_idx <- complete.cases(x_bg1)
            x_bg1 <- x_bg1[complete_idx, , drop = FALSE]
            bg_cood_tmp <- bg_cood_tmp[complete_idx, , drop = FALSE]

            if (nrow(x_bg1) == 0) {
                warning("No valid background points after env annotation")
                next
            }

            message(sprintf("  After env annotation: %d valid points", nrow(x_bg1)))

            # Subsample to target if we have more
            if (nrow(bg_cood_tmp) > n_bg_per_species) {
                idx <- sample(nrow(bg_cood_tmp), n_bg_per_species)
                bg_cood_tmp <- bg_cood_tmp[idx, , drop = FALSE]
                x_bg1 <- x_bg1[idx, , drop = FALSE]
            } else if (nrow(bg_cood_tmp) < n_bg_per_species) {
                warning(sprintf("Only got %d/%d BG points", nrow(bg_cood_tmp), n_bg_per_species))
            }

            x_bg1 <- process_x_matrix(x_bg1)

            # Create y_bg with zeros, marking which species this BG is for
            y_bg_tmp <- as.data.frame(matrix(0,
                nrow = nrow(x_bg1), ncol = ncol(y),
                dimnames = list(NULL, colnames(y))
            ))
            y_bg_tmp$bg <- species_idx

            # Store in lists
            x_bg_list[[species_idx]] <- x_bg1
            y_bg_list[[species_idx]] <- y_bg_tmp
            cood_bg_list[[species_idx]] <- bg_cood_tmp

            message(sprintf("  Generated %d BG points", nrow(x_bg1)))
        },
        error = function(e) {
            warning(sprintf("Error processing species %d: %s", species_idx, e$message))
        }
    )
}

# Combine all BG data
message(sprintf(
    "BG lists collected: x=%d, y=%d, cood=%d",
    length(x_bg_list), length(y_bg_list), length(cood_bg_list)
))
# Remove NULL entries from lists
x_bg_list <- x_bg_list[!vapply(x_bg_list, is.null, logical(1))]
y_bg_list <- y_bg_list[!vapply(y_bg_list, is.null, logical(1))]
cood_bg_list <- cood_bg_list[!vapply(cood_bg_list, is.null, logical(1))]
message(sprintf(
    "After removing NULLs: x=%d, y=%d, cood=%d",
    length(x_bg_list), length(y_bg_list), length(cood_bg_list)
))

if (length(x_bg_list) > 0) {
    x_bg <- do.call(rbind, x_bg_list)
    y_bg <- do.call(rbind, y_bg_list)
    cood_bg <- do.call(rbind, cood_bg_list)

    message(sprintf("BG data combined: %d rows", nrow(x_bg)))
    message(sprintf(
        "y_bg has bg column: %s, values: %s",
        "bg" %in% colnames(y_bg),
        paste(head(unique(y_bg$bg)), collapse = ", ")
    ))

    # Append BG to data
    y <- rbind(as.data.frame(y), as.data.frame(y_bg))
    x <- rbind(as.data.frame(x), as.data.frame(x_bg))
    cood <- rbind(as.data.frame(cood), as.data.frame(cood_bg))
    message(sprintf(
        "After appending BG: y has %d rows, bg column sum = %d",
        nrow(y), sum(y$bg > 0)
    ))
} else {
    warning("No background points generated for any species")
}

# ---- Add site key ----
x$site <- y$site <- cood$site <- seq_len(nrow(y))

# ---- Scale environmental covariates ----
message("Scaling environmental covariates...")
nenv_vars <- c("Intercept", "site")
env_vars <- setdiff(colnames(x), nenv_vars)
xenv <- x[, env_vars, drop = FALSE]

rex_means <- apply(xenv, 2, mean, na.rm = TRUE)
rex_sds <- apply(xenv, 2, sd, na.rm = TRUE)
scales_df <- data.frame(variable = colnames(xenv), mean = rex_means, sd = rex_sds)

for (i in seq_len(ncol(xenv))) {
    varname <- colnames(xenv)[i]
    mean_val <- rex_means[i]
    sd_val <- rex_sds[i]
    if (!is.finite(sd_val) || sd_val == 0) {
        sd_val <- 1
        message(sprintf("Variable %s has zero sd; setting to 1", varname))
    }
    xenv[, i] <- (xenv[, i] - mean_val) / sd_val
}
x[, env_vars] <- xenv

# ---- Convert to long format ----
message("Converting to long format...")
species_names <- sps
y_sponly <- as.matrix(y[, species_names, drop = FALSE])

nz <- which(y_sponly != 0, arr.ind = TRUE)
y_count <- data.frame(
    site    = nz[, 1],
    species = nz[, 2],
    count   = y_sponly[nz]
)
message(sprintf("Presence records in y_count: %d", nrow(y_count)))

# Add BG rows (count=0 for indicated species)
# bg column indicates which species index this BG point belongs to
bg_rows <- which(y$bg > 0)
message(sprintf("Background rows found: %d", length(bg_rows)))
if (length(bg_rows) > 0) {
    bg_df <- data.frame(
        site = bg_rows,
        species = y$bg[bg_rows],
        count = 0
    )
    y_count <- rbind(y_count, bg_df)
    message(sprintf("Added %d background records to y_count", nrow(bg_df)))
}
message(sprintf("Total records in y_count: %d", nrow(y_count)))

# ---- Species index map ----
species_index_map <- setNames(species_names, seq_along(species_names))

# ---- Save outputs ----
message("Saving outputs...")
out_dir <- file.path(dpath, CLUSTER)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Save extent
ex_df <- data.frame(
    xmin = ex[1], xmax = ex[2],
    ymin = ex[3], ymax = ex[4]
)
write.csv(ex_df, file.path(out_dir, paste0(CLUSTER, "_extent.csv")), row.names = FALSE)

# Save scales
write.csv(scales_df, file.path(out_dir, paste0(CLUSTER, "_env_scales.csv")), row.names = FALSE)

# Save missing species (if any)
if (length(missing_species) > 0) {
    missing_df <- data.frame(species = missing_species)
    write.csv(missing_df, file.path(out_dir, paste0(CLUSTER, "_missing_species.csv")), row.names = FALSE)
    message(sprintf("Wrote %d missing species to CSV", length(missing_species)))
}

# Save main data
store <- list(
    y = y_count,
    x = x,
    cood = cood,
    tree = tree,
    species_index_map = species_index_map,
    sps = sps
)
save(store, file = file.path(out_dir, paste0(CLUSTER, "_run_files.Rdata")))

message(sprintf("Done! Output saved to: %s", out_dir))
message(sprintf("Final: %d sites, %d species", nrow(x), length(sps)))

# Optional sanity checks (plots)
# Plot raw data
# raw_df <- raw
# raw_df$presence <- 1
# world <- rnaturalearth::ne_countries(returnclass = "sf")
# world_vect <- terra::vect(world)
# ggplot() +
#     geom_sf(data = world, fill = "lightgray", color = "white") +
#     geom_point(data = raw_df, aes(x = lon, y = lat, color = species), alpha = 0.5) +
#     theme_minimal() +
#     ggtitle(sprintf("Raw occurrences for cluster %s", CLUSTER))

# # Plot prepared presence/background data for one species at a time
# p <- list()
# for (i in seq_along(sps)) {
#     sp <- sps[i]
#     sp_idx <- which(species_names == sp)
#     sp_y <- y_count[y_count$species == sp_idx, ]
#     sp_y$presence <- ifelse(sp_y$count > 0, "Presence", "Background")
#     sp_cood <- cood[sp_y$site, , drop = FALSE]

#     plot_df <- data.frame(
#         lon = sp_cood$lon,
#         lat = sp_cood$lat,
#         presence = sp_y$presence
#     )

#     p[[i]] <- ggplot() +
#         geom_sf(data = world, fill = "lightgray", color = "white") +
#         geom_point(data = plot_df, aes(x = lon, y = lat, color = presence), alpha = 0.5) +
#         theme_minimal() +
#         ggtitle(sprintf("Prepared data for species %s in cluster %s", sp, CLUSTER))
# }

# ---- Clean up ----
rm_objs <- c("data", "data_df", "data_sps", "data_wide", "raw", "raw_list")
rm_objs <- rm_objs[rm_objs %in% ls()]
if (length(rm_objs) > 0) rm(list = rm_objs)
gc()

# End of script
