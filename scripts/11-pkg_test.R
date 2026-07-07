# Purpose: Package TEST data for model evaluation
# This script processes a single cluster for a specified replicate
# Input: cluster_name and rep number passed as command line argument or set manually
# Output: test_data .Rdata files saved to raw_data/<DATASET>/<CLUSTER>/

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
    CLUSTER <- "Hyla1"
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
original_sps <- store$sps
species_index_map <- store$species_index_map

# ---- Diagnostic: Check store structure ----
message(sprintf(
    "  store_x: %d rows, rownames[1:5]=%s",
    nrow(store_x), paste(head(rownames(store_x), 5), collapse = ",")
))
message(sprintf("  store_x has 'site' column: %s", "site" %in% colnames(store_x)))
message(sprintf(
    "  store_cood: %d rows, rownames[1:5]=%s",
    nrow(store_cood), paste(head(rownames(store_cood), 5), collapse = ",")
))
message(sprintf("  store_y$site range: [%d, %d]", min(store_y$site), max(store_y$site)))

# ---- Load training model_data to get trained species ----
model_data_file <- file.path(out_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_model_data.Rdata"))
if (!file.exists(model_data_file)) {
    stop(sprintf("Training model_data not found: %s. Run training data script first.", model_data_file))
}
load(model_data_file) # loads model_data
trained_species_names <- model_data$species_names
trained_distance_matrix <- model_data$D_phylo
trained_J <- model_data$J
message(sprintf("Loaded training model_data: %d trained species", length(trained_species_names)))

# ---- Process specified replicate ----
message(sprintf("Processing replicate %d...", REPNO))

# Get train/test indices for this rep
idx <- list_of_indices[[paste0("rep_", REPNO)]]
idx_test <- sort(idx$test)

# Subset y to testing sites only
y <- store_y[store_y$site %in% idx_test, ]
message(sprintf("  Test observations before species sync: %d", nrow(y)))

# ---- Synchronize species: keep only species that were in training data ----
# Map y$species (original index) to species names using species_index_map
y_orig_species_idx <- y$species
y_species_names <- as.character(species_index_map[as.character(y_orig_species_idx)])

# Filter to only observations of trained species
trained_mask <- y_species_names %in% trained_species_names
n_dropped_obs <- sum(!trained_mask)
if (n_dropped_obs > 0) {
    dropped_species <- unique(y_species_names[!trained_mask])
    message(sprintf(
        "  WARNING: Dropping %d observations of %d species not in training data",
        n_dropped_obs, length(dropped_species)
    ))
    message(sprintf("  Dropped species: %s", paste(dropped_species, collapse = ", ")))
    y <- y[trained_mask, ]
    y_species_names <- y_species_names[trained_mask]
}

# Use trained species list and distance matrix
sps <- trained_species_names
distance_matrix <- trained_distance_matrix

# Remap species indices in y to match training data indices (1:n_trained_species)
name_to_new_idx <- setNames(seq_along(sps), sps)
y$species <- as.integer(name_to_new_idx[y_species_names])

# Validate species remapping
stopifnot(all(y$species >= 1 & y$species <= length(sps)))

message(sprintf(
    "  After sync: %d trained species, %d test observations, species range=[%d,%d]",
    length(sps), nrow(y), min(y$species), max(y$species)
))

# Get unique sites in y (preserving original site IDs)
unique_sites <- unique(y$site)
message(sprintf("  Unique test sites: %d", length(unique_sites)))

# ---- Subset x and cood using the site column for matching ----
# store_x and store_cood have a 'site' column that matches y$site values
# Use this for reliable matching instead of rownames

if ("site" %in% colnames(store_x)) {
    # Match by site column
    x_site_match <- match(unique_sites, store_x$site)
    if (any(is.na(x_site_match))) {
        stop(sprintf(
            "Some test sites not found in store_x: %s",
            paste(unique_sites[is.na(x_site_match)], collapse = ", ")
        ))
    }
    x <- store_x[x_site_match, , drop = FALSE]

    cood_site_match <- match(unique_sites, store_cood$site)
    if (any(is.na(cood_site_match))) {
        stop(sprintf(
            "Some test sites not found in store_cood: %s",
            paste(unique_sites[is.na(cood_site_match)], collapse = ", ")
        ))
    }
    cood <- store_cood[cood_site_match, , drop = FALSE]
} else {
    # Fallback: try rowname matching
    message("  Note: store_x has no 'site' column, using rowname matching")
    x <- store_x[as.character(unique_sites), , drop = FALSE]
    cood <- store_cood[as.character(unique_sites), , drop = FALSE]
}

# Validate subsetting worked
stopifnot(nrow(x) == length(unique_sites))
stopifnot(nrow(cood) == length(unique_sites))

message(sprintf("  Subsetted x: %d rows, cood: %d rows", nrow(x), nrow(cood)))

# ---- Add quadratic terms for bio1 and bio13 ----
if ("CHELSA_bio_1" %in% colnames(x)) {
    x <- cbind(x, CHELSA_bio_1_sq = x[, "CHELSA_bio_1"]^2)
    message("  Added CHELSA_bio_1^2")
}
if ("CHELSA_bio_13" %in% colnames(x)) {
    x <- cbind(x, CHELSA_bio_13_sq = x[, "CHELSA_bio_13"]^2)
    message("  Added CHELSA_bio_13^2")
}

# ---- Create site mapping: original site ID -> new contiguous index ----
# This maps unique_sites[i] -> i
site_id_to_new_idx <- setNames(seq_along(unique_sites), as.character(unique_sites))

# ---- Extract soft clip values for each observation (optional) ----
if (USE_SOFT_CLIPS) {
    message("  Extracting soft clip values...")
    soft_clip_values <- numeric(nrow(y))
    soft_clip_found <- logical(nrow(y))

    for (obs_idx in seq_len(nrow(y))) {
        sp_idx <- y$species[obs_idx]
        sp_name <- sps[sp_idx]
        site_id <- y$site[obs_idx] # Original site ID

        # Map original site ID to row index in subsetted cood
        cood_row_idx <- site_id_to_new_idx[as.character(site_id)]

        if (is.na(cood_row_idx)) {
            warning(sprintf("Site %d not found in mapping, using default soft_clip=1", site_id))
            soft_clip_values[obs_idx] <- 1
            soft_clip_found[obs_idx] <- FALSE
            next
        }

        # Get coordinates for this site
        site_cood <- cood[cood_row_idx, , drop = FALSE]

        # Load soft clip raster for this species
        soft_clip_file <- file.path(soft_clip_dir, paste0(sp_name, "_soft_clip.tif"))
        if (file.exists(soft_clip_file)) {
            soft_clip_rast <- terra::rast(soft_clip_file)
            val <- terra::extract(soft_clip_rast, cbind(site_cood$lon, site_cood$lat))[1, 1]
            soft_clip_values[obs_idx] <- ifelse(is.na(val), 1, val)
            soft_clip_found[obs_idx] <- TRUE
        } else {
            soft_clip_values[obs_idx] <- 1
            soft_clip_found[obs_idx] <- FALSE
        }
    }

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

# ---- Map site IDs in y to new contiguous indices ----
y$site_new <- site_id_to_new_idx[as.character(y$site)]

# Validate site mapping
stopifnot(!any(is.na(y$site_new)))
stopifnot(all(y$site_new >= 1 & y$site_new <= length(unique_sites)))

# Extract y components
y_species_idx <- y$species
y_values <- y$count
y_sites <- y$site_new

# Convert to presence/absence if THIN = TRUE
if (THIN) {
    y_values <- ifelse(y_values > 0, 1L, 0L)
}

# Validate y values include zeros (background points)
if (min(y_values) != 0) {
    warning("y_values does not include zeros - no background points in test set?")
}

# ---- Create test data matrix ----
rmv <- which(colnames(x) %in% c("site"))
if (length(rmv) > 0) x <- x[, -rmv, drop = FALSE]

colnames(x) <- c(
    "meanTemp", "tempSeason", "precipWet", "precipSeason", "cloudCover",
    "EVI", "TRI", "elevation", "Intercept", "meanTemp2", "precipWet2"
)

col_order <- c(
    "Intercept", "meanTemp", "meanTemp2", "tempSeason",
    "precipWet", "precipWet2", "precipSeason", "cloudCover", "EVI",
    "TRI", "elevation"
)
x <- x[, col_order]
x_matrix <- as.matrix(x)

# Remove site column from cood if present
if ("site" %in% colnames(cood)) {
    cood <- cood[, !colnames(cood) %in% "site", drop = FALSE]
}

# Reset rownames to 1:N to match site_new indexing
rownames(x_matrix) <- seq_len(nrow(x_matrix))
rownames(cood) <- seq_len(nrow(cood))

# ---- Create test_data list ----
# J should match training J to ensure species indices are compatible
test_data <- list(
    N = nrow(x_matrix),
    J = length(sps), # Use training J, not length(unique(y_species_idx))
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

# ---- Final validation ----
message("\n  --- Final Validation ---")
message(sprintf("  X: %d rows x %d cols", nrow(test_data$X), ncol(test_data$X)))
message(sprintf("  cood: %d rows x %d cols", nrow(test_data$cood), ncol(test_data$cood)))
message(sprintf("  N_obs: %d", test_data$N_obs))
message(sprintf("  N (unique sites): %d", test_data$N))
message(sprintf("  J (trained species): %d", test_data$J))
message(sprintf("  Unique species in test: %d", length(unique(test_data$species))))
message(sprintf("  Site range: [%d, %d]", min(test_data$site), max(test_data$site)))
message(sprintf("  Species range: [%d, %d]", min(test_data$species), max(test_data$species)))

# Validate dimensions match
stopifnot(test_data$N == nrow(test_data$X))
stopifnot(test_data$N == nrow(test_data$cood))
stopifnot(test_data$N_obs == length(test_data$y))
stopifnot(test_data$N_obs == length(test_data$species))
stopifnot(test_data$N_obs == length(test_data$site))
stopifnot(test_data$N_obs == length(test_data$offset))
stopifnot(max(test_data$site) <= test_data$N)
stopifnot(max(test_data$species) <= test_data$J)

# ---- Save test_data ----
out_file <- file.path(out_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_test_data.Rdata"))
save(test_data, file = out_file)
message(sprintf("\n  Saved: %s", basename(out_file)))

message(sprintf("Done! Processed rep %d for cluster %s", REPNO, CLUSTER))


# # ---- Diagnostic plots for test data ----
# library(ggplot2)

# # Load world map
# world_file <- file.path(root, "raw_data", "world.Rdata")
# if (file.exists(world_file)) {
#     load(world_file) # loads 'world' as sf object
# } else {
#     world <- rnaturalearth::ne_countries(returnclass = "sf")
# }

# # Plot all test sites (all species combined)
# all_plot_df <- data.frame(
#     lon = test_data$cood$lon[test_data$site],
#     lat = test_data$cood$lat[test_data$site],
#     presence = ifelse(test_data$y > 0, "Presence", "Background")
# )

# p_all <- ggplot() +
#     geom_sf(data = world, fill = "lightgray", color = "white") +
#     geom_point(data = all_plot_df, aes(x = lon, y = lat, color = presence), alpha = 0.5, size = 1) +
#     coord_sf(
#         xlim = range(test_data$cood$lon) + c(-2, 2),
#         ylim = range(test_data$cood$lat) + c(-2, 2)
#     ) +
#     scale_color_manual(values = c("Presence" = "blue", "Background" = "red")) +
#     theme_minimal() +
#     ggtitle(sprintf("All TEST sites for cluster %s (rep %d)", CLUSTER, REPNO))
# print(p_all)

# # Plot per-species test sites
# p <- list()
# for (i in seq_along(test_data$species_names)) {
#     sp_name <- test_data$species_names[i]
#     sp_idx <- i  # Species index matches position in species_names

#     # Get observations for this species
#     sp_mask <- test_data$species == sp_idx

#     if (sum(sp_mask) == 0) {
#         message(sprintf("Species %d (%s): No test observations, skipping plot", i, sp_name))
#         next
#     }

#     # Get site indices for this species' observations
#     sp_sites <- test_data$site[sp_mask]
#     sp_y <- test_data$y[sp_mask]

#     plot_df <- data.frame(
#         lon = test_data$cood$lon[sp_sites],
#         lat = test_data$cood$lat[sp_sites],
#         presence = ifelse(sp_y > 0, "Presence", "Background")
#     )

#     n_pres <- sum(sp_y > 0)
#     n_bg <- sum(sp_y == 0)

#     p[[i]] <- ggplot() +
#         geom_sf(data = world, fill = "lightgray", color = "white") +
#         geom_point(data = plot_df, aes(x = lon, y = lat, color = presence), alpha = 0.6, size = 2) +
#         coord_sf(
#             xlim = range(plot_df$lon) + c(-3, 3),
#             ylim = range(plot_df$lat) + c(-3, 3)
#         ) +
#         scale_color_manual(values = c("Presence" = "blue", "Background" = "red")) +
#         theme_minimal() +
#         theme(panel.grid = element_blank()) +
#         ggtitle(sprintf("TEST: %s (%d pres, %d bg) - Cluster %s Rep %d",
#                         sp_name, n_pres, n_bg, CLUSTER, REPNO))

#     print(p[[i]])
# }

# # Optional: Compare train vs test for a single species
# compare_train_test <- function(species_idx, model_data, test_data, world) {
#     sp_name <- test_data$species_names[species_idx]

#     # Training data
#     train_mask <- model_data$species == species_idx
#     train_sites <- model_data$site[train_mask]

#     # Handle both matrix and data.frame for cood
#     if (is.matrix(model_data$cood) || is.data.frame(model_data$cood)) {
#         train_lon <- model_data$cood[train_sites, "lon"]
#         train_lat <- model_data$cood[train_sites, "lat"]
#     } else {
#         train_lon <- model_data$cood$lon[train_sites]
#         train_lat <- model_data$cood$lat[train_sites]
#     }

#     train_df <- data.frame(
#         lon = train_lon,
#         lat = train_lat,
#         presence = ifelse(model_data$y[train_mask] > 0, "Presence", "Background"),
#         split = "Train"
#     )

#     # Test data
#     test_mask <- test_data$species == species_idx
#     test_sites <- test_data$site[test_mask]

#     # Handle both matrix and data.frame for cood
#     if (is.matrix(test_data$cood) || is.data.frame(test_data$cood)) {
#         test_lon <- test_data$cood[test_sites, "lon"]
#         test_lat <- test_data$cood[test_sites, "lat"]
#     } else {
#         test_lon <- test_data$cood$lon[test_sites]
#         test_lat <- test_data$cood$lat[test_sites]
#     }

#     test_df <- data.frame(
#         lon = test_lon,
#         lat = test_lat,
#         presence = ifelse(test_data$y[test_mask] > 0, "Presence", "Background"),
#         split = "Test"
#     )

#     combined_df <- rbind(train_df, test_df)

#     ggplot() +
#         geom_sf(data = world, fill = "lightgray", color = "white") +
#         geom_point(data = combined_df,
#                    aes(x = lon, y = lat, color = presence, shape = split),
#                    alpha = 0.6, size = 2) +
#         coord_sf(
#             xlim = range(combined_df$lon, na.rm = TRUE) + c(-3, 3),
#             ylim = range(combined_df$lat, na.rm = TRUE) + c(-3, 3)
#         ) +
#         scale_color_manual(values = c("Presence" = "blue", "Background" = "red")) +
#         scale_shape_manual(values = c("Train" = 16, "Test" = 17)) +
#         theme_minimal() +
#         theme(panel.grid = element_blank()) +
#         ggtitle(sprintf("Train vs Test: %s - Cluster %s Rep %d", sp_name, CLUSTER, REPNO))
# }

# # Example: compare train/test for first species
# print(compare_train_test(6, model_data, test_data, world))
