# Purpose: Generate train/test split indices for all clusters
# Input: DATASET passed as command line argument or set manually
# Output: Index list saved to raw_data/<DATASET>/<CLUSTER>/<CLUSTER>_indices.Rdata
# This script is pretty lightweight so it's actually not designed for DSQ, just run it interactively or as one job
# Clear workspace
rm(list = ls())


# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    message("Running on HPC")
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
} else {
    message("Running locally")
    root <- "~/clim_risk_phylosdm"
}

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
checks_dir <- file.path(root, "checks")
dir.create(checks_dir, showWarnings = FALSE, recursive = TRUE)


# ---- Configuration ----
NREPS <- 5
TRAIN_SPLIT <- 0.70 # 70% train, 30% test

# ---- Load cluster list ----
load(file.path(dpath, "spList.Rdata")) # loads spList
if (CLUSTER == "ALL") {
    clusters <- names(spList)
} else {
    clusters <- CLUSTER
}

message(sprintf("Found %d clusters to process", length(clusters)))

# Track failed clusters
failed_clusters <- data.frame(
    cluster = character(),
    reason = character(),
    stringsAsFactors = FALSE
)

# ---- Loop over clusters ----
for (i in seq_along(clusters)) {
    CLUSTER <- clusters[i]
    message(sprintf("\n[%d/%d] Processing cluster: %s", i, length(clusters), CLUSTER))

    tryCatch(
        {
            out_dir <- file.path(dpath, CLUSTER)

            # Load cluster data
            data_file <- file.path(out_dir, paste0(CLUSTER, "_run_files.Rdata"))
            if (!file.exists(data_file)) {
                warning(sprintf("Data file not found: %s. Skipping.", data_file))
                failed_clusters <- rbind(failed_clusters, data.frame(
                    cluster = CLUSTER,
                    reason = "Data file not found",
                    stringsAsFactors = FALSE
                ))
                next
            }
            load(data_file) # loads 'store' object
            message(sprintf("  Loaded data: %d sites, %d species", nrow(store$x), length(store$sps)))

            # Create 70/30 train/test split indices
            n <- nrow(store$x)
            n_train <- floor(n * TRAIN_SPLIT)

            set.seed(42) # For reproducibility
            list_of_indices <- list()

            for (j in seq_len(NREPS)) {
                # Random sample of row indices for training
                train_idx <- sample(seq_len(n), size = n_train, replace = FALSE)
                test_idx <- setdiff(seq_len(n), train_idx)

                list_of_indices[[j]] <- list(
                    train = train_idx,
                    test = test_idx
                )
                message(sprintf("    Rep %d: %d train, %d test", j, length(train_idx), length(test_idx)))
            }
            names(list_of_indices) <- paste0("rep_", seq_len(NREPS))

            # Save indices
            out_file <- file.path(out_dir, paste0(CLUSTER, "_indices.Rdata"))
            save(list_of_indices, file = out_file)
            message(sprintf("  Saved %d reps to: %s", NREPS, basename(out_file)))

            # Clean up
            rm(store, list_of_indices)
            gc()
        },
        error = function(e) {
            warning(sprintf("Error processing cluster %s: %s", CLUSTER, e$message))
            failed_clusters <<- rbind(failed_clusters, data.frame(
                cluster = CLUSTER,
                reason = e$message,
                stringsAsFactors = FALSE
            ))
        }
    )
}

# ---- Write failed clusters to checks/ ----
if (nrow(failed_clusters) > 0) {
    failed_file <- file.path(checks_dir, paste0("04-data_indices_failed_clusters_", DATASET, ".csv"))
    write.csv(failed_clusters, failed_file, row.names = FALSE)
    message(sprintf("\nWrote %d failed clusters to: %s", nrow(failed_clusters), failed_file))
} else {
    message("\nNo failed clusters!")
}

message(sprintf("\nDone! Processed %d clusters.", length(clusters)))


# Diagnostic plots
# Uncomment to visualize the first split
# library(ggplot2)

# yc <- merge(store$cood, store$y, by = "site")
# yc_train <- yc[which(yc$site %in% list_of_indices$rep_1$train), ]
# world <- rnaturalearth::ne_countries(returnclass = "sf")

# ggplot(data = yc_train) +
#     geom_point(aes(x = lon, y = lat, color = as.factor(ifelse
#     (count > 1, 1, 0))), alpha = 0.5, size = .5) +
#     geom_sf(data = world, fill = NA, color = "grey50") +
#     coord_sf(xlim = range(store$cood$lon), ylim = range(store$cood$lat)) +
#     labs(color = "Presence") +
#     theme_minimal() +
#     ggtitle(sprintf("Cluster: %s - Training Sites (Rep 1)", CLUSTER))

# yc_test <- yc[which(yc$site %in% list_of_indices$rep_1$test), ]

# # And test --
# ggplot(data = yc_test) +
#     geom_point(aes(x = lon, y = lat, color = as.factor(ifelse
#     (count > 1, 1, 0))), alpha = 0.5, size = .5) +
#     geom_sf(data = world, fill = NA, color = "grey50") +
#     coord_sf(xlim = range(store$cood$lon), ylim = range(store$cood$lat)) +
#     labs(color = "Presence") +
#     theme_minimal() +
#     ggtitle(sprintf("Cluster: %s - Test Sites (Rep 1)", CLUSTER))


# # For a given species
# index <- 8
# yc_train_sp <- yc_train[which(yc_train$species == index), ]

# ggplot(data = yc_train_sp) +
#     geom_point(aes(x = lon, y = lat, color = as.factor(ifelse
#     (count > 1, 1, 0))), alpha = 0.5, size = 2) +
#     geom_sf(data = world, fill = "grey50", color = "grey50", alpha = .1) +
#     coord_sf(xlim = range(yc_train_sp$lon) + c(-4, 4), ylim = range(yc_train_sp$lat) + c(-4, 4)) +
#     labs(color = "Presence") +
#     theme_minimal() +
#     theme(panel.grid = element_blank()) +
#     ggtitle(sprintf("Cluster: %s - Training Sites (Rep 1), for a given species %d", CLUSTER, index))

# yc_test_sp <- yc_test[which(yc_test$species == index), ]

# ggplot(data = yc_test_sp) +
#     geom_point(aes(x = lon, y = lat, color = as.factor(ifelse
#     (count > 1, 1, 0))), alpha = 0.5, size = 2) +
#     geom_sf(data = world, fill = "grey50", color = "grey50", alpha = .1) +
#     coord_sf(xlim = range(yc_train_sp$lon) + c(-4, 4), ylim = range(yc_train_sp$lat) + c(-4, 4)) +
#     labs(color = "Presence") +
#     theme_minimal() +
#     theme(panel.grid = element_blank()) +
#     ggtitle(sprintf("Cluster: %s - Test Sites (Rep 1), for a given species %d", CLUSTER, index))
