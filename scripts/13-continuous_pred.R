# ====================== CLIM_RISK_PHYLOSDM | SPATIAL PREDICTION ======================
# Runs spatial prediction for all species in a cluster using best-performing model
# Uses PER-SPECIES extents to avoid memory/disk issues with disjunct distributions
# ============================================================================

rm(list = ls())
suppressPackageStartupMessages({
    library(terra)
    library(rstan)
})

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    # Use SLURM_JOB_ID if available, otherwise fallback to PID
    job_id <- Sys.getenv("SLURM_JOB_ID")
    if (job_id == "") job_id <- as.character(Sys.getpid())
    temp_dir <- file.path("/vast/palmer/scratch/jetz/ss4224/clim_risk", paste0("spatial_dir_", job_id))
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
    epath <- "/vast/palmer/pi/jetz/ss4224/env"
    message("Running on HPC")

    # Create temp directory if needed
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)


    # Set environment and terra options
    Sys.setenv(TMPDIR = temp_dir)
    terra::terraOptions(
        tempdir = temp_dir,
        memfrac = 0.1,
        todisk = TRUE
    )

    # ---- Clean up old terra temp files at start ----
    message("Cleaning up old terra temp files...")
    old_temps <- list.files(temp_dir, pattern = "^spat_", full.names = TRUE)
    if (length(old_temps) > 0) {
        unlink(old_temps)
        message(sprintf("  Removed %d old temp files", length(old_temps)))
    }
} else {
    root <- "~/clim_risk_phylosdm"
    epath <- "~/env"
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
message(sprintf("Spatial prediction for cluster: %s, rep: %d", CLUSTER, REPNO))

# ---- Configuration ----
BUFFER_DEG <- 5
MAX_EXTENT_DEG <- 80 # Maximum extent in any direction (degrees)
MAX_CELLS <- 200e6 # Maximum cells (~200 million) before warning
WRITE_RELATIVE <- TRUE
USE_SOFT_CLIPS <- TRUE

# ---- Paths ----
scripts_directory <- file.path(root, "scripts")
raw_dir <- file.path(root, "raw_data")
res_dir <- file.path(root, "res")
analysis_dir <- file.path(root, "analysis")
dpath <- file.path(raw_dir, DATASET)
data_dir <- file.path(dpath, CLUSTER)

# Output directories - use PI storage for outputs
if (HPC != "FALSE") {
    pred_dir <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/v0/spatial_pred"
} else {
    pred_dir <- file.path(analysis_dir, EXP_ROOT, "spatial_pred")
}
soft_clip_dir <- file.path(analysis_dir, "soft_clips")
eval_dir <- file.path(analysis_dir, EXP_ROOT, "eval")

dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

# ---- File naming ----
base_prefix <- paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO)

# ---- Environment files ----
env_files <- c(
    "CHELSA_bio_1.tif",
    "CHELSA_bio_4.tif",
    "CHELSA_bio_13.tif",
    "CHELSA_bio_15.tif",
    "cloudCover.tif",
    "Annual_EVI.tif",
    "TRI.tif",
    "elevation_1KMmean_SRTM.tif"
)

# Model variable order (must match model_data$X column order exactly)
model_var_order <- c(
    "Intercept", "meanTemp", "meanTemp2", "tempSeason",
    "precipWet", "precipWet2", "precipSeason", "cloudCover", "EVI",
    "TRI", "elevation"
)

# ---- Helper functions ----
exp_marginal_over_f <- function(sigma_vec) {
    exp(0.5 * (sigma_vec^2))
}

write_rast <- function(r, path) {
    terra::writeRaster(r, path, overwrite = TRUE)
    message("  Wrote: ", basename(path))
}

get_species_extent <- function(cood, buffer_deg = 5, max_extent_deg = 80) {
    # Compute extent from coordinates
    lon_range <- range(cood[, "lon"], na.rm = TRUE)
    lat_range <- range(cood[, "lat"], na.rm = TRUE)

    # Add buffer
    lon_min <- lon_range[1] - buffer_deg
    lon_max <- lon_range[2] + buffer_deg
    lat_min <- lat_range[1] - buffer_deg
    lat_max <- lat_range[2] + buffer_deg

    # Check extent size
    lon_span <- lon_max - lon_min
    lat_span <- lat_max - lat_min

    # Cap if too large
    if (lon_span > max_extent_deg) {
        lon_center <- mean(c(lon_min, lon_max))
        lon_min <- lon_center - max_extent_deg / 2
        lon_max <- lon_center + max_extent_deg / 2
        message(sprintf("  Capped longitude span from %.1f to %.1f degrees", lon_span, max_extent_deg))
    }

    if (lat_span > max_extent_deg) {
        lat_center <- mean(c(lat_min, lat_max))
        lat_min <- lat_center - max_extent_deg / 2
        lat_max <- lat_center + max_extent_deg / 2
        message(sprintf("  Capped latitude span from %.1f to %.1f degrees", lat_span, max_extent_deg))
    }

    terra::ext(lon_min, lon_max, lat_min, lat_max)
}

load_and_prepare_env_rasters <- function(extent, env_files, epath, scales_df, model_var_order) {
    # Load and crop all rasters for this extent
    env_rast_list <- list()
    for (i in seq_along(env_files)) {
        f <- file.path(epath, env_files[i])
        if (!file.exists(f)) {
            stop(sprintf("Environment file not found: %s", f))
        }
        r <- terra::rast(f)
        r <- terra::crop(r, extent, snap = "out")
        env_rast_list[[i]] <- r
    }

    # Use first raster as template and resample others to match
    template_rast <- env_rast_list[[1]]
    for (i in 2:length(env_rast_list)) {
        if (!terra::compareGeom(template_rast, env_rast_list[[i]], stopOnError = FALSE)) {
            env_rast_list[[i]] <- terra::resample(env_rast_list[[i]], template_rast, method = "bilinear")
        }
    }

    # Stack
    r_env_raw <- terra::rast(env_rast_list)
    names(r_env_raw) <- c(
        "meanTemp", "tempSeason", "precipWet", "precipSeason",
        "cloudCover", "EVI", "TRI", "elevation"
    )

    # Scale
    r_env_scaled <- r_env_raw
    for (v in names(r_env_scaled)) {
        scale_row <- which(scales_df$variable == v)
        if (length(scale_row) == 0) next
        mu <- scales_df$mean[scale_row]
        sd_val <- scales_df$sd[scale_row]
        if (sd_val == 0) sd_val <- 1
        r_env_scaled[[v]] <- (r_env_scaled[[v]] - mu) / sd_val
    }

    # Add quadratic terms
    r_env_scaled$meanTemp2 <- r_env_scaled$meanTemp^2
    r_env_scaled$precipWet2 <- r_env_scaled$precipWet^2

    # Add intercept
    r_intercept <- terra::rast(r_env_scaled[[1]])
    terra::values(r_intercept) <- 1
    names(r_intercept) <- "Intercept"

    # Combine and reorder
    r_env_all <- c(r_intercept, r_env_scaled)
    r_env_final <- r_env_all[[model_var_order]]

    return(r_env_final)
}

cleanup_terra_temps <- function(temp_dir = NULL) {
    # Get temp directory
    if (is.null(temp_dir)) {
        temp_dir <- terra::terraOptions()$tempdir
    }

    # Remove terra temp files manually
    if (!is.null(temp_dir) && dir.exists(temp_dir)) {
        temp_files <- list.files(temp_dir, pattern = "^spat_", full.names = TRUE)
        if (length(temp_files) > 0) {
            unlink(temp_files)
        }
    }
    # Also try terra's built-in cleanup (wrapped in tryCatch)
    tryCatch(
        {
            terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
        },
        error = function(e) {
            # Silently ignore
        }
    )

    gc()
}

# ---- Load cluster data ----
run_file <- file.path(data_dir, paste0(CLUSTER, "_run_files.Rdata"))
if (!file.exists(run_file)) {
    stop(sprintf("Run file not found: %s", run_file))
}
load(run_file) # loads 'store'
message(sprintf("Loaded store: %d sites, %d species", nrow(store$x), length(store$sps)))

species_names <- store$sps
species_index_map <- store$species_index_map
J <- length(species_names)

# ---- Load model_data to get trained species and verify variable order ----
model_data_file <- file.path(data_dir, paste0(base_prefix, "_model_data.Rdata"))
if (!file.exists(model_data_file)) {
    stop(sprintf("Model data file not found: %s", model_data_file))
}
load(model_data_file) # loads 'model_data'
trained_species <- model_data$species_names
K <- model_data$K
message(sprintf("Loaded model_data: %d trained species, K=%d", length(trained_species), K))
message(sprintf("Model variable order: %s", paste(colnames(model_data$X), collapse = ", ")))

# Verify our assumed order matches
stopifnot(all(colnames(model_data$X) == model_var_order))

# ---- Load scaling stats ----
scales_file <- file.path(data_dir, paste0(CLUSTER, "_env_scales.csv"))
if (!file.exists(scales_file)) {
    stop(sprintf("Scales file not found: %s", scales_file))
}
scales_df <- read.csv(scales_file, stringsAsFactors = FALSE)

# Rename scales_df variables to match model names
scales_df$variable <- c(
    "meanTemp", "tempSeason", "precipWet", "precipSeason",
    "cloudCover", "EVI", "TRI", "elevation"
)
message("Loaded scaling statistics")

# ---- Load model fit ----
result_file <- file.path(res_dir, paste0(base_prefix, "_", MODEL_TYPE, "_fit.Rdata"))
if (!file.exists(result_file)) {
    stop(sprintf("Result file not found: %s", result_file))
}
load(result_file) # loads 'result'
stan_fit <- result$fit
message("Loaded model fit")

# ---- Load conditional predictions ----
cond_pred_dir <- file.path(analysis_dir, EXP_ROOT, "cond_pred")
cond_pred_file <- file.path(cond_pred_dir, paste0(base_prefix, "_cond_pred.Rdata"))
has_cond_pred <- file.exists(cond_pred_file)
if (has_cond_pred) {
    load(cond_pred_file) # loads 'B_cond_list'
    message(sprintf("Loaded conditional predictions: %d species", length(B_cond_list)))
} else {
    message("Conditional predictions not found, will use model coefficients only")
    B_cond_list <- NULL
}

# ---- Load best model selection ----
best_model_file <- file.path(eval_dir, paste0(base_prefix, "_best_model.csv"))
if (file.exists(best_model_file)) {
    best_model_df <- read.csv(best_model_file, stringsAsFactors = FALSE)
    message(sprintf("Loaded best model selection: %d species", nrow(best_model_df)))
} else {
    message("Best model file not found, defaulting to 'model' for all species")
    best_model_df <- data.frame(
        species_name = trained_species,
        best_model = rep("model", length(trained_species)),
        stringsAsFactors = FALSE
    )
}

# ---- Extract posterior samples ----
posterior_B <- rstan::extract(stan_fit, pars = "B")$B # (S, J, K)
posterior_sigma_f <- rstan::extract(stan_fit, pars = "sigma_f")$sigma_f # (S)
S <- dim(posterior_B)[1]
message(sprintf(
    "Posterior: S=%d samples, J=%d species, K=%d predictors",
    S, dim(posterior_B)[2], dim(posterior_B)[3]
))

# ---- Correction factor for random effect ----
cf_bar <- mean(exp_marginal_over_f(posterior_sigma_f))
message(sprintf("Mean correction factor (exp(0.5*sigma_f^2)): %.4f", cf_bar))

# ---- Main prediction loop (PER-SPECIES EXTENT) ----
message(sprintf("\n========== Starting predictions for %d trained species ==========", length(trained_species)))

for (j in seq_along(trained_species)) {
    sp_name <- trained_species[j]
    message(sprintf("\n[%d/%d] Processing species: %s", j, length(trained_species), sp_name))

    # Check if output already exists
    out_file <- file.path(pred_dir, paste0(base_prefix, "_", sp_name, "_relprob.tif"))
    if (file.exists(out_file)) {
        message("  Output already exists, skipping")
        next
    }

    tryCatch(
        {
            # ---- Get species-specific extent ----
            # Find sites where this species occurs
            sp_idx_in_store <- which(species_names == sp_name)
            if (length(sp_idx_in_store) == 0) {
                message("  Species not found in store, skipping")
                next
            }

            sp_obs_mask <- store$y$species == sp_idx_in_store
            sp_sites <- unique(store$y$site[sp_obs_mask])
            sp_cood <- store$cood[sp_sites, , drop = FALSE]

            message(sprintf("  Species has %d unique sites", length(sp_sites)))

            # Compute species-specific extent
            sp_extent <- get_species_extent(sp_cood, buffer_deg = BUFFER_DEG, max_extent_deg = MAX_EXTENT_DEG)

            lon_span <- sp_extent$xmax - sp_extent$xmin
            lat_span <- sp_extent$ymax - sp_extent$ymin
            est_cells <- (lon_span * 120) * (lat_span * 120) # ~1km resolution

            message(sprintf(
                "  Extent: lon [%.2f, %.2f], lat [%.2f, %.2f]",
                sp_extent$xmin, sp_extent$xmax, sp_extent$ymin, sp_extent$ymax
            ))
            message(sprintf(
                "  Extent size: %.1f x %.1f deg (~%.0f million cells)",
                lon_span, lat_span, est_cells / 1e6
            ))

            if (est_cells > MAX_CELLS) {
                warning(sprintf("  Large extent: %.0f million cells. Processing anyway but may be slow.", est_cells / 1e6))
            }

            # ---- Load and prepare environment rasters for this extent ----
            message("  Loading environment rasters...")
            r_env_final <- load_and_prepare_env_rasters(
                sp_extent, env_files, epath, scales_df, model_var_order
            )

            message(sprintf(
                "  Env raster: %d x %d cells, %d layers",
                nrow(r_env_final), ncol(r_env_final), terra::nlyr(r_env_final)
            ))

            # Verify dimensions
            stopifnot(terra::nlyr(r_env_final) == K)

            # ---- Determine which model to use ----
            best_row <- which(best_model_df$species_name == sp_name)
            if (length(best_row) == 0) {
                best_model <- "model"
                message("  Species not in best_model_df, defaulting to 'model'")
            } else {
                best_model <- best_model_df$best_model[best_row]
            }
            message(sprintf("  Using: %s", best_model))

            # ---- Get coefficients based on best model ----
            if (best_model == "cond_pred" && has_cond_pred && sp_name %in% names(B_cond_list)) {
                B_cond_sp <- as.matrix(B_cond_list[[sp_name]])
                beta_bar <- as.vector(colMeans(B_cond_sp))
                message("  Using conditional prediction coefficients")
            } else {
                beta_bar <- as.vector(colMeans(posterior_B[, j, , drop = FALSE]))
                message("  Using model coefficients")
            }

            message(sprintf(
                "  beta_bar: length=%d, range=[%.3f, %.3f]",
                length(beta_bar), min(beta_bar), max(beta_bar)
            ))

            if (length(beta_bar) != K) {
                stop(sprintf("Coefficient length mismatch: got %d, expected %d", length(beta_bar), K))
            }

            # ---- Compute linear predictor using weighted sum ----
            message("  Computing linear predictor...")
            linpred <- r_env_final[[1]] * beta_bar[1]
            for (k in 2:K) {
                linpred <- linpred + r_env_final[[k]] * beta_bar[k]
            }
            names(linpred) <- "linpred"

            # Free up env raster memory
            rm(r_env_final)
            gc()

            # ---- Add soft clip if available ----
            if (USE_SOFT_CLIPS) {
                soft_clip_file <- file.path(soft_clip_dir, paste0(sp_name, "_soft_clip.tif"))
                if (file.exists(soft_clip_file)) {
                    soft_clip <- terra::rast(soft_clip_file)
                    soft_clip <- terra::crop(soft_clip, sp_extent, snap = "out")
                    soft_clip <- terra::resample(soft_clip, linpred, method = "bilinear")
                    soft_clip_log <- log(soft_clip + 1e-6)
                    linpred <- linpred + soft_clip_log
                    rm(soft_clip, soft_clip_log)
                    message("  Added soft clip offset")
                } else {
                    message("  No soft clip file found")
                }
            }

            # ---- Expected intensity: exp(linpred) * correction_factor ----
            lambda_mean <- exp(linpred) * cf_bar
            rm(linpred)

            # ---- Relative probability ----
            if (WRITE_RELATIVE) {
                maxv <- terra::global(lambda_mean, "max", na.rm = TRUE)[1, 1]
                if (!is.na(maxv) && maxv > 0) {
                    rel_prob <- lambda_mean / maxv
                } else {
                    warning("  Max lambda is NA or 0, skipping relative scaling")
                    rel_prob <- lambda_mean
                }
            } else {
                rel_prob <- lambda_mean
            }
            rm(lambda_mean)

            # ---- Write output ----
            write_rast(rel_prob, out_file)

            # ---- Print summary stats ----
            stats <- terra::global(rel_prob, c("min", "mean", "max"), na.rm = TRUE)
            message(sprintf(
                "  Stats: min=%.4f, mean=%.4f, max=%.4f",
                stats$min, stats$mean, stats$max
            ))

            # ---- Clean up ----
            rm(rel_prob)
            cleanup_terra_temps()

            message("  Done!")
        },
        error = function(e) {
            message(sprintf("  ERROR: %s", e$message))
            # Clean up even on error
            cleanup_terra_temps()
        }
    )
}

# ---- Final cleanup ----
message("\nFinal cleanup...")
cleanup_terra_temps()

message(sprintf("\n========== Done! Predictions saved to: %s ==========", pred_dir))
