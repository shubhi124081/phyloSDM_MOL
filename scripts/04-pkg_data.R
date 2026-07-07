# Purpose: Package raw data for model fitting (creates stan_data objects)
# This script processes a single cluster for a specified replicate
# Input: cluster_name and rep number passed as command line argument or set manually
# Output: stan_data .Rdata files saved to raw_data/<DATASET>/<CLUSTER>/

# Clear workspace
rm(list = ls())

# Libraries
library(terra)

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
    epath <- "/vast/palmer/pi/jetz/ss4224/env"
    message("Running on HPC")
    Sys.setenv(TMPDIR = "/vast/palmer/scratch/jetz/ss4224")
    terra::terraOptions(tempdir = "/vast/palmer/scratch/jetz/ss4224")
} else {
    root <- "~/clim_risk_phylosdm"
    epath <- "~/env"
    message("Running locally")
}

if (interactive()) {
    EXP_ROOT <- "v0"
    EXP_ID <- "sub1000"
    DATASET <- "amphibians"
    CLUSTER <- "Rani1"
    FSP <- "ALL"
    REPNO <- 1
    NREP <- 1
    MODEL_TYPE <- "TMB"
    MODEL_NAME <- "05_lgcp_corrected.cpp"
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
message(sprintf("Processing cluster: %s, rep: %d", CLUSTER, REPNO))

# ---- Paths ----
dpath <- file.path(root, "raw_data", DATASET)
out_dir <- file.path(dpath, CLUSTER)
soft_clip_dir <- file.path(root, "analysis", "soft_clips")

# ---- Configuration ----
THIN <- TRUE # Convert counts to presence/absence (1/0)
USE_SOFT_CLIPS <- TRUE # Whether to use soft clip values as offset

# ---- Load cluster data ----
run_file <- file.path(out_dir, paste0(CLUSTER, "_run_files.Rdata"))
if (!file.exists(run_file)) {
    stop(sprintf("Run file not found: %s. Run 03-gen_data.R first.", run_file))
}
load(run_file) # loads 'store' object
message(sprintf("Loaded data: %d sites, %d species", nrow(store$x), length(store$sps)))

# ---- Load indices ----
indices_file <- file.path(out_dir, paste0(CLUSTER, "_indices.Rdata"))
if (!file.exists(indices_file)) {
    stop(sprintf("Indices file not found: %s. Run 04-data_indices.R first.", indices_file))
}
load(indices_file) # loads 'list_of_indices'

# ---- Extract store components ----
store_y <- store$y
store_x <- store$x
store_cood <- store$cood
tree <- store$tree
sps <- store$sps
species_index_map <- store$species_index_map

# ---- Create phylogenetic distance matrix ----
if (inherits(tree, "phylo")) {
    distance_matrix <- ape::cophenetic.phylo(tree)
    # Reorder to match sps order (species indices in y map to positions in sps)
    distance_matrix <- distance_matrix[sps, sps, drop = FALSE]
} else {
    distance_matrix <- tree
    # Ensure order matches sps
    if (!is.null(rownames(distance_matrix))) {
        distance_matrix <- distance_matrix[sps, sps, drop = FALSE]
    }
}

# ---- Process specified replicate ----
message(sprintf("Processing replicate %d...", REPNO))

# Get train/test indices for this rep
idx <- list_of_indices[[paste0("rep_", REPNO)]]
idx_train <- sort(idx$train)

# Subset y to training sites only
y <- store_y[store_y$site %in% idx_train, ]

# ---- Synchronize species: drop any with no observations in training data ----
observed_species_idx <- sort(unique(y$species))
n_observed <- length(observed_species_idx)
n_total <- length(sps)

if (n_observed < n_total) {
    # Use species_index_map to map original indices to names
    # species_index_map: named vector, names = original index, values = species name
    observed_orig_idx <- as.character(observed_species_idx)
    observed_species_names <- as.character(species_index_map[observed_orig_idx])
    missing_idx <- setdiff(names(species_index_map), observed_orig_idx)
    missing_species_names <- as.character(species_index_map[missing_idx])
    message(sprintf(
        "  WARNING: %d of %d species have no training observations",
        n_total - n_observed, n_total
    ))
    message(sprintf("  Dropping species: %s", paste(missing_species_names, collapse = ", ")))

    # Subset sps to only observed species (by name)
    sps <- observed_species_names

    # Subset distance_matrix to only observed species (by name)
    distance_matrix <- distance_matrix[sps, sps, drop = FALSE]

    # Remap species indices in y to contiguous 1:n_observed
    name_to_new_idx <- setNames(seq_along(sps), sps)
    # Map y$species (original index) to species name, then to new index
    y_species_names <- as.character(species_index_map[as.character(y$species)])
    y$species <- as.integer(name_to_new_idx[y_species_names])

    message(sprintf(
        "  After sync: %d species, D_phylo=%dx%d, species range=[%d,%d]",
        length(sps), nrow(distance_matrix), ncol(distance_matrix),
        min(y$species), max(y$species)
    ))
}

# Get unique sites in y (in order they appear)
unique_sites <- unique(y$site)

# Subset x and cood to matching sites
x <- store_x[unique_sites, , drop = FALSE]
cood <- store_cood[unique_sites, , drop = FALSE]

# ---- Add quadratic terms for bio1 and bio13 ----
if ("CHELSA_bio_1" %in% colnames(x)) {
    x <- cbind(x, CHELSA_bio_1_sq = x[, "CHELSA_bio_1"]^2)
    message("  Added CHELSA_bio_1^2")
}
if ("CHELSA_bio_13" %in% colnames(x)) {
    x <- cbind(x, CHELSA_bio_13_sq = x[, "CHELSA_bio_13"]^2)
    message("  Added CHELSA_bio_13^2")
}

# ---- Extract soft clip values for each observation (optional) ----
if (USE_SOFT_CLIPS) {
    # soft_clip is species-specific, so we need to get it per observation in y
    message("  Extracting soft clip values...")
    soft_clip_values <- numeric(nrow(y))
    soft_clip_found <- logical(nrow(y))

    for (obs_idx in seq_len(nrow(y))) {
        sp_idx <- y$species[obs_idx]
        sp_name <- sps[sp_idx]
        site_idx <- y$site[obs_idx]

        # Get coordinates for this site
        site_cood <- cood[as.character(site_idx), , drop = FALSE]

        # Load soft clip raster for this species
        soft_clip_file <- file.path(soft_clip_dir, paste0(sp_name, "_soft_clip.tif"))
        if (file.exists(soft_clip_file)) {
            soft_clip_rast <- terra::rast(soft_clip_file)
            # Extract value at site location
            val <- terra::extract(soft_clip_rast, cbind(site_cood$lon, site_cood$lat))[1, 1]
            soft_clip_values[obs_idx] <- ifelse(is.na(val), 1, val)
            soft_clip_found[obs_idx] <- TRUE
        } else {
            soft_clip_values[obs_idx] <- 1 # Default to 1 if no soft clip
            soft_clip_found[obs_idx] <- FALSE
        }
    }

    # Log transform soft clip values (for use as offset)
    # Add small value to avoid log(0)
    # If soft_clip file wasn't found, set log to 0 (equivalent to no offset)
    soft_clip_log <- log(pmax(soft_clip_values, 1e-6))
    soft_clip_log[!soft_clip_found] <- 0
    message(sprintf(
        "  Soft clip range: [%.4f, %.4f], log range: [%.2f, %.2f]",
        min(soft_clip_values), max(soft_clip_values),
        min(soft_clip_log), max(soft_clip_log)
    ))
    message(sprintf("  Soft clip files found for %d/%d observations", sum(soft_clip_found), length(soft_clip_found)))
} else {
    message("  Skipping soft clip extraction (USE_SOFT_CLIPS = FALSE)")
    soft_clip_log <- rep(0, nrow(y))
}

# Create new site index mapping (1:n_sites)
site_map <- setNames(seq_along(unique_sites), unique_sites)
y$site_new <- site_map[as.character(y$site)]

# Extract y components
y_species_idx <- y$species
y_values <- y$count
y_sites <- y$site_new

# Convert to presence/absence if THIN = TRUE
if (THIN) {
    y_values <- ifelse(y_values > 0, 1L, 0L)
}

# Validate y values
if (min(y_values) != 0) {
    stop("y_values should include zeros for presence-absence data.")
}

# ---- Create stan data ----
rmv <- which(colnames(x) == "site")
if (length(rmv) > 0) x <- x[, -c(rmv)]
colnames(x) <- c(
    "meanTemp", "tempSeason", "precipWet", "precipSeason", "cloudCover",
    "EVI", "TRI", "elevation", "Intercept", "meanTemp2", "precipWet2"
)

order <- c(
    "Intercept", "meanTemp", "meanTemp2", "tempSeason",
    "precipWet", "precipWet2", "precipSeason", "cloudCover", "EVI",
    "TRI", "elevation"
)
x <- x[, order]
x_matrix <- as.matrix(x)

model_data <- list(
    N = nrow(x_matrix),
    J = length(unique(y_species_idx)),
    K = ncol(x_matrix),
    N_obs = length(y_values),
    species = as.integer(y_species_idx),
    site = as.integer(y_sites),
    X = x_matrix,
    y = as.integer(y_values),
    D_phylo = distance_matrix,
    tree = tree,
    cood = cood,
    species_names = sps,
    offset = soft_clip_log,
    species_index_map = store$species_index_map
)

# ---- Save stan_data ----
out_file <- file.path(out_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_model_data.Rdata"))
save(model_data, file = out_file)
message(sprintf("  Saved: %s", basename(out_file)))

message(sprintf("Done! Processed rep %d for cluster %s", REPNO, CLUSTER))
