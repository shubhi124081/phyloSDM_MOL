# ====================== CLIM_RISK_PHYLOSDM | THRESHOLD COMPUTATION ======================
# Computes per-species thresholds (and AUC, sens, spec) from spatial prediction rasters
# Uses test data to determine optimal thresholds for binary classification
# ========================================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(terra)
    library(pROC)
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
    root <- "~/phyloSDM_MOL"
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
message(sprintf("Threshold computation for cluster: %s, rep: %d", CLUSTER, REPNO))

# ---- Paths ----
raw_dir <- file.path(root, "raw_data")
analysis_dir <- file.path(root, "analysis")
dpath <- file.path(raw_dir, DATASET)
data_dir <- file.path(dpath, CLUSTER)

# Input/output directories
pred_dir <- file.path(analysis_dir, EXP_ROOT, "spatial_pred")
eval_dir <- file.path(analysis_dir, EXP_ROOT, "eval")
threshold_dir <- file.path(eval_dir, "thresholds")

dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)

# ---- File naming ----
base_prefix <- paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO)

# ---- Helper functions ----

# Compute ROC/AUC and best-threshold metrics
roc_metrics <- function(obs_binary, pred_numeric) {
    # Remove NAs
    valid_idx <- !is.na(obs_binary) & !is.na(pred_numeric)
    obs_binary <- obs_binary[valid_idx]
    pred_numeric <- pred_numeric[valid_idx]

    if (length(obs_binary) < 2 || length(unique(obs_binary)) < 2) {
        return(list(auc = NA, sens = NA, spec = NA, thr = NA))
    }

    roc_obj <- pROC::roc(response = obs_binary, predictor = pred_numeric, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    # Get best threshold using Youden's J
    best_coords <- pROC::coords(roc_obj, "best",
        best.method = "youden",
        ret = c("threshold", "sensitivity", "specificity")
    )

    # Handle case where multiple "best" thresholds are returned
    if (is.data.frame(best_coords) || is.matrix(best_coords)) {
        thr_val <- best_coords$threshold[1]
        sens_val <- best_coords$sensitivity[1]
        spec_val <- best_coords$specificity[1]
    } else {
        thr_val <- best_coords["threshold"]
        sens_val <- best_coords["sensitivity"]
        spec_val <- best_coords["specificity"]
    }

    list(
        auc = auc_val,
        sens = as.numeric(sens_val),
        spec = as.numeric(spec_val),
        thr = as.numeric(thr_val)
    )
}

# ---- Load test data ----
test_data_file <- file.path(data_dir, paste0(base_prefix, "_test_data.Rdata"))
if (!file.exists(test_data_file)) {
    stop(sprintf("Test data file not found: %s", test_data_file))
}
load(test_data_file) # loads 'test_data'
message(sprintf(
    "Loaded test data: N=%d sites, J=%d species, N_obs=%d observations",
    test_data$N, test_data$J, test_data$N_obs
))

# Extract test data components
species_names <- test_data$species_names
cood_test <- test_data$cood
y_test <- test_data$y
y_species <- test_data$species
y_sites <- test_data$site

J <- length(species_names)

# ---- Initialize result containers ----
auc_values <- rep(NA_real_, J)
sensitivity_values <- rep(NA_real_, J)
specificity_values <- rep(NA_real_, J)
threshold_values <- rep(NA_real_, J)
n_test_obs <- rep(0L, J)
n_presences <- rep(0L, J)

# ---- Main loop over species ----
message(sprintf("\n========== Computing thresholds for %d species ==========", J))

for (j in seq_len(J)) {
    sp_name <- species_names[j]
    message(sprintf("\n[%d/%d] Processing species: %s", j, J, sp_name))

    tryCatch(
        {
            # Get test observations for this species
            sp_mask <- y_species == j
            n_obs <- sum(sp_mask)
            n_test_obs[j] <- n_obs

            if (n_obs == 0) {
                message("  No test observations, skipping")
                next
            }

            # Get site indices and observed values
            sp_sites <- y_sites[sp_mask]
            sp_y <- y_test[sp_mask]
            n_presences[j] <- sum(sp_y > 0)

            message(sprintf(
                "  Test obs: %d (%d presences, %d absences)",
                n_obs, n_presences[j], n_obs - n_presences[j]
            ))

            # Check we have both presences and absences
            if (n_presences[j] == 0 || n_presences[j] == n_obs) {
                message("  Need both presences and absences for ROC, skipping")
                next
            }

            # Get coordinates for test sites
            sp_cood <- cood_test[sp_sites, , drop = FALSE]

            # Load prediction raster
            pred_file <- file.path(pred_dir, paste0(base_prefix, "_", sp_name, "_relprob.tif"))
            if (!file.exists(pred_file)) {
                message(sprintf("  Prediction raster not found: %s", basename(pred_file)))
                next
            }

            pred_rast <- terra::rast(pred_file)

            # Extract predicted values at test locations
            pred_vals <- terra::extract(pred_rast,
                cbind(sp_cood[, "lon"], sp_cood[, "lat"]),
                method = "simple"
            )
            pred_numeric <- as.numeric(pred_vals[[1]])

            # Check for valid predictions
            n_valid <- sum(!is.na(pred_numeric))
            if (n_valid < n_obs) {
                message(sprintf("  Warning: %d/%d predictions are NA", n_obs - n_valid, n_obs))
            }
            if (n_valid < 2) {
                message("  Not enough valid predictions, skipping")
                next
            }

            # Compute ROC metrics
            obs_binary <- as.integer(sp_y > 0)
            metrics <- roc_metrics(obs_binary, pred_numeric)

            auc_values[j] <- metrics$auc
            sensitivity_values[j] <- metrics$sens
            specificity_values[j] <- metrics$spec
            if (is.na(metrics$thr) | !is.finite(metrics$thr)) {
                threshold_values[j] <- 0
            } else {
                threshold_values[j] <- metrics$thr
            }


            message(sprintf(
                "  AUC=%.3f, Threshold=%.4f, Sens=%.3f, Spec=%.3f",
                metrics$auc, metrics$thr, metrics$sens, metrics$spec
            ))

            # Clean up
            rm(pred_rast, pred_vals)
            gc()
        },
        error = function(e) {
            message(sprintf("  ERROR: %s", e$message))
        }
    )
}

# ---- Compile results ----
threshold_df <- data.frame(
    EXP_ROOT = rep(EXP_ROOT, J),
    EXP_ID = rep(EXP_ID, J),
    CLUSTER = rep(CLUSTER, J),
    REPNO = rep(REPNO, J),
    species_idx = seq_len(J),
    species_name = species_names,
    n_test_obs = n_test_obs,
    n_presences = n_presences,
    AUC = auc_values,
    sensitivity = sensitivity_values,
    specificity = specificity_values,
    threshold = threshold_values
)
if (any(is.na(threshold_df$threshold))) {
    message("Found NA thresholds: replacing with 0")
    which_na <- which(is.na(threshold_df$threshold))
    threshold_df[which_na, "threshold"] <- 0
}

# ---- Save results ----
# Save as CSV
out_csv <- file.path(threshold_dir, paste0(base_prefix, "_thresholds.csv"))
write.csv(threshold_df, out_csv, row.names = FALSE)
message(sprintf("\nSaved thresholds CSV: %s", out_csv))

# Save as Rdata
out_rdata <- file.path(threshold_dir, paste0(base_prefix, "_thresholds.Rdata"))
save(threshold_df, file = out_rdata)
message(sprintf("Saved thresholds Rdata: %s", out_rdata))

# ---- Print summary ----
message("\n========== Threshold Summary ==========")
message(sprintf("Species processed: %d", J))
message(sprintf("Species with valid thresholds: %d", sum(!is.na(threshold_values))))
message(sprintf(
    "Mean AUC: %.3f (SD: %.3f)",
    mean(auc_values, na.rm = TRUE), sd(auc_values, na.rm = TRUE)
))
message(sprintf(
    "AUC range: [%.3f, %.3f]",
    min(auc_values, na.rm = TRUE), max(auc_values, na.rm = TRUE)
))
message(sprintf(
    "Mean threshold: %.4f (SD: %.4f)",
    mean(threshold_values, na.rm = TRUE), sd(threshold_values, na.rm = TRUE)
))

# Print per-species summary
message("\nPer-species results:")
print(threshold_df[, c("species_name", "n_test_obs", "n_presences", "AUC", "threshold")])

message(sprintf("\n========== Done! Thresholds saved to: %s ==========", threshold_dir))
