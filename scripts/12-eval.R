rm(list = ls())

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
    message("Running on HPC")
} else {
    root <- "~/phyloSDM_MOL"
    message("Running locally")
}

scripts_directory <- file.path(root, "scripts")
res_dir <- file.path(root, "res")
raw_dir <- file.path(root, "raw_data")
analysis_dir <- file.path(root, "analysis")

# ---- Libraries ----
library(mvtnorm)
library(pROC)


# ---- Functions ----
predict_intensity_rf <- function(posterior_B, posterior_sigma_f, x_new, y_new, offset_new) {
    # Extract dimensions
    num_samples <- dim(posterior_B)[1]
    N_new <- nrow(x_new)

    log_lambda_samples <- matrix(NA, nrow = num_samples, ncol = N_new)
    pb <- txtProgressBar(min = 0, max = num_samples, style = 3)

    for (i in 1:num_samples) {
        for (n in 1:N_new) {
            species_n <- y_new$species[n] # species index
            f_new <- rnorm(1, mean = 0, sd = posterior_sigma_f[i])
            log_lambda_samples[i, n] <- sum(x_new[n, ] * posterior_B[i, species_n, ]) +
                f_new + offset_new[n]
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)

    # Return mean predicted intensity (lambda)
    lambda_samples <- exp(log_lambda_samples)
    predicted_intensity <- colMeans(lambda_samples)

    return(predicted_intensity)
}

predict_intensity_cond <- function(B_cond_df, posterior_sigma_f, x_new, offset_new) {
    # Uses conditional B samples from cp_LGCP_simple output
    # B_cond_df: data.frame with S rows (posterior draws) and K columns (predictors)
    B_cond_mat <- as.matrix(B_cond_df) # (S, K)
    num_samples <- nrow(B_cond_mat)
    N_new <- nrow(x_new)

    # Ensure we have matching number of posterior samples
    if (num_samples != length(posterior_sigma_f)) {
        # Subsample or recycle to match
        sample_idx <- seq_len(min(num_samples, length(posterior_sigma_f)))
        B_cond_mat <- B_cond_mat[sample_idx, , drop = FALSE]
        posterior_sigma_f <- posterior_sigma_f[sample_idx]
        num_samples <- length(sample_idx)
    }

    log_lambda_samples <- matrix(NA, nrow = num_samples, ncol = N_new)
    pb <- txtProgressBar(min = 0, max = num_samples, style = 3)

    for (i in 1:num_samples) {
        for (n in 1:N_new) {
            f_new <- rnorm(1, mean = 0, sd = posterior_sigma_f[i])
            log_lambda_samples[i, n] <- sum(x_new[n, ] * B_cond_mat[i, ]) +
                f_new + offset_new[n]
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)

    # Return mean predicted intensity (lambda)
    lambda_samples <- exp(log_lambda_samples)
    predicted_intensity <- colMeans(lambda_samples)

    return(predicted_intensity)
}

evaluate_species <- function(predicted_intensity, y_test_species_values, cood_test_species,
                             pred_dir, filename_prefix, species_name, ii, n_species) {
    # Scale predictions
    predicted_intensity_scaled <- predicted_intensity / sum(predicted_intensity, na.rm = TRUE)

    # Binary response
    test_binary <- ifelse(y_test_species_values > 0, 1, 0)

    # Calculate AUC + associated metrics
    roc_obj <- pROC::roc(response = test_binary, predictor = predicted_intensity_scaled, quiet = TRUE)
    auc_val <- pROC::auc(roc_obj)

    best_thresh <- pROC::coords(roc_obj, "best", best.method = "youden", ret = "threshold")
    thresh_val <- as.numeric(best_thresh[1])

    sens <- pROC::coords(roc_obj, x = "best", ret = "sensitivity")
    sens_val <- as.numeric(sens[1])

    specs <- pROC::coords(roc_obj, x = "best", ret = "specificity")
    spec_val <- as.numeric(specs[1])

    # Boyce Index (placeholder - set to NA)
    boyce_val <- NA

    # Confusion matrix
    pred_binary <- ifelse(predicted_intensity_scaled >= thresh_val, 1, 0)
    conf_mat <- table(
        Predicted = factor(pred_binary, levels = c(0, 1)),
        Observed = factor(test_binary, levels = c(0, 1))
    )

    tp_val <- conf_mat["1", "1"]
    tn_val <- conf_mat["0", "0"]
    fp_val <- conf_mat["1", "0"]
    fn_val <- conf_mat["0", "1"]

    # Save predictions to CSV
    pred_df <- data.frame(
        cood_test_species,
        true_binary = test_binary,
        pred_binary = pred_binary,
        pred_intensity = predicted_intensity_scaled
    )

    pred_filename <- paste0(filename_prefix, "_", species_name, "_predictions.csv")
    write.csv(pred_df, file = file.path(pred_dir, pred_filename), row.names = FALSE)

    message(sprintf("  Species %d/%d (%s): AUC=%.3f", ii, n_species, species_name, auc_val))

    return(list(
        auc = auc_val,
        boyce = boyce_val,
        threshold = thresh_val,
        sensitivity = sens_val,
        specificity = spec_val,
        TP = tp_val,
        TN = tn_val,
        FP = fp_val,
        FN = fn_val
    ))
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
message(sprintf("Evaluation for cluster: %s, rep: %d", CLUSTER, REPNO))

# ---- Paths ----
dpath <- file.path(raw_dir, DATASET)
data_dir <- file.path(dpath, CLUSTER)
eval_dir <- file.path(analysis_dir, EXP_ROOT, "eval")
pred_dir_model <- file.path(eval_dir, "predictions_model")
pred_dir_cond <- file.path(eval_dir, "predictions_cond")

if (!dir.exists(eval_dir)) {
    dir.create(eval_dir, recursive = TRUE)
    message(sprintf("Created eval directory: %s", eval_dir))
}
if (!dir.exists(pred_dir_model)) {
    dir.create(pred_dir_model, recursive = TRUE)
    message(sprintf("Created model predictions directory: %s", pred_dir_model))
}
if (!dir.exists(pred_dir_cond)) {
    dir.create(pred_dir_cond, recursive = TRUE)
    message(sprintf("Created conditional predictions directory: %s", pred_dir_cond))
}

# ---- Load test data ----
test_data_file <- file.path(data_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_test_data.Rdata"))
if (!file.exists(test_data_file)) {
    stop(sprintf("Test data file not found: %s", test_data_file))
}
load(test_data_file) # loads 'test_data'
message(sprintf(
    "Loaded test data: N=%d, J=%d, K=%d, N_obs=%d",
    test_data$N, test_data$J, test_data$K, test_data$N_obs
))

# ---- Load model fit ----
result_file <- file.path(res_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_", MODEL_TYPE, "_fit.Rdata"))
if (!file.exists(result_file)) {
    stop(sprintf("Result file not found: %s", result_file))
}
load(result_file) # loads 'result'
message("Loaded model fit")

# ---- Load conditional predictions ----
# Conditional predictions are saved to: analysis/<EXP_ROOT>/cond_pred/
cond_pred_dir <- file.path(analysis_dir, EXP_ROOT, "cond_pred")
cond_pred_file <- file.path(cond_pred_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_cond_pred.Rdata"))
has_cond_pred <- file.exists(cond_pred_file)
if (has_cond_pred) {
    load(cond_pred_file) # loads 'B_cond_list'
    message(sprintf("Loaded conditional predictions: %d species", length(B_cond_list)))
} else {
    message(sprintf("Conditional predictions file not found: %s", cond_pred_file))
    message("Will skip conditional evaluation")
}

# ---- Extract test data components ----
x_test <- test_data$X
y_test_values <- test_data$y
y_test_species <- test_data$species
y_test_sites <- test_data$site
cood_test <- test_data$cood
offset_test <- test_data$offset
species_names <- test_data$species_names
species_index_map <- test_data$species_index_map

# Dimensions
N_test <- test_data$N
K <- test_data$K
J <- test_data$J

# ---- Extract posterior samples ----
# result$posterior (local/cmdstanr fits) is already extracted plain arrays;
# result$fit (HPC/rstan fits) is a live stanfit object needing rstan::extract().
if (!is.null(result$posterior)) {
    posterior_B <- result$posterior$B # (num_samples, J, K)
    posterior_sigma_f <- result$posterior$sigma_f # (num_samples)
} else {
    posterior_B <- rstan::extract(result$fit, pars = "B")$B # (num_samples, J, K)
    posterior_sigma_f <- rstan::extract(result$fit, pars = "sigma_f")$sigma_f # (num_samples)
}

# ---- Get unique species indices present in test data ----
species_indices <- sort(unique(y_test_species))
n_species <- length(species_indices)
message(sprintf("Evaluating %d species", n_species))

# ---- Initialize result vectors for MODEL evaluation ----
model_auc <- numeric(n_species)
model_boyce <- numeric(n_species)
model_threshold <- numeric(n_species)
model_sensitivity <- numeric(n_species)
model_specificity <- numeric(n_species)
model_TN <- numeric(n_species)
model_TP <- numeric(n_species)
model_FP <- numeric(n_species)
model_FN <- numeric(n_species)

# ---- Initialize result vectors for CONDITIONAL evaluation ----
cond_auc <- numeric(n_species)
cond_boyce <- numeric(n_species)
cond_threshold <- numeric(n_species)
cond_sensitivity <- numeric(n_species)
cond_specificity <- numeric(n_species)
cond_TN <- numeric(n_species)
cond_TP <- numeric(n_species)
cond_FP <- numeric(n_species)
cond_FN <- numeric(n_species)

# ---- Filename prefixes ----
base_prefix <- paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO)
model_prefix <- paste0(base_prefix, "_model")
cond_prefix <- paste0(base_prefix, "_cond")

# ============================================================
# LOOP 1: MODEL-BASED EVALUATION
# ============================================================
message("\n========== MODEL-BASED EVALUATION ==========")

for (ii in seq_len(n_species)) {
    species_idx <- species_indices[ii]
    species_name <- species_names[species_idx]

    tryCatch(
        {
            # Subset to this species (observation-level mask)
            sp_mask <- y_test_species == species_idx

            # Get the site indices for this species' observations
            sites_for_species <- y_test_sites[sp_mask]

            # Subset observation-level data
            y_test_species_values <- y_test_values[sp_mask]
            offset_test_species <- offset_test[sp_mask]

            # Subset site-level data using site indices
            x_test_species <- x_test[sites_for_species, , drop = FALSE]
            cood_test_species <- cood_test[sites_for_species, , drop = FALSE]

            # Create y_test_species dataframe for predict function compatibility
            y_test_sp_df <- data.frame(
                species = y_test_species[sp_mask],
                count = y_test_species_values
            )

            # Predict intensity using model parameters
            predicted_intensity <- predict_intensity_rf(
                posterior_B, posterior_sigma_f,
                x_test_species, y_test_sp_df, offset_test_species
            )

            # Evaluate
            eval_results <- evaluate_species(
                predicted_intensity, y_test_species_values, cood_test_species,
                pred_dir_model, model_prefix, species_name, ii, n_species
            )

            model_auc[ii] <- eval_results$auc
            model_boyce[ii] <- eval_results$boyce
            model_threshold[ii] <- eval_results$threshold
            model_sensitivity[ii] <- eval_results$sensitivity
            model_specificity[ii] <- eval_results$specificity
            model_TP[ii] <- eval_results$TP
            model_TN[ii] <- eval_results$TN
            model_FP[ii] <- eval_results$FP
            model_FN[ii] <- eval_results$FN
        },
        error = function(e) {
            message(sprintf("  Error processing species %d (%s): %s", species_idx, species_name, e$message))
            model_auc[ii] <<- NA
            model_threshold[ii] <<- NA
            model_boyce[ii] <<- NA
            model_sensitivity[ii] <<- NA
            model_specificity[ii] <<- NA
            model_TP[ii] <<- NA
            model_TN[ii] <<- NA
            model_FP[ii] <<- NA
            model_FN[ii] <<- NA
        }
    )
}

# ============================================================
# LOOP 2: CONDITIONAL PREDICTION EVALUATION
# ============================================================
if (has_cond_pred) {
    message("\n========== CONDITIONAL PREDICTION EVALUATION ==========")

    for (ii in seq_len(n_species)) {
        species_idx <- species_indices[ii]
        species_name <- species_names[species_idx]

        tryCatch(
            {
                # Check if this species has conditional predictions
                if (!species_name %in% names(B_cond_list)) {
                    message(sprintf("  Species %s not found in B_cond_list, skipping", species_name))
                    cond_auc[ii] <- NA
                    next
                }

                # Subset to this species (observation-level mask)
                sp_mask <- y_test_species == species_idx

                # Get the site indices for this species' observations
                sites_for_species <- y_test_sites[sp_mask]

                # Subset observation-level data
                y_test_species_values <- y_test_values[sp_mask]
                offset_test_species <- offset_test[sp_mask]

                # Subset site-level data using site indices
                x_test_species <- x_test[sites_for_species, , drop = FALSE]
                cood_test_species <- cood_test[sites_for_species, , drop = FALSE]

                # Get conditional B for this species
                # B_cond_list is a named list, each element is a data.frame (S rows, K cols)
                B_cond_species <- B_cond_list[[species_name]]

                # Predict intensity using conditional parameters
                predicted_intensity <- predict_intensity_cond(
                    B_cond_species, posterior_sigma_f,
                    x_test_species, offset_test_species
                )

                # Evaluate
                eval_results <- evaluate_species(
                    predicted_intensity, y_test_species_values, cood_test_species,
                    pred_dir_cond, cond_prefix, species_name, ii, n_species
                )

                cond_auc[ii] <- eval_results$auc
                cond_boyce[ii] <- eval_results$boyce
                cond_threshold[ii] <- eval_results$threshold
                cond_sensitivity[ii] <- eval_results$sensitivity
                cond_specificity[ii] <- eval_results$specificity
                cond_TP[ii] <- eval_results$TP
                cond_TN[ii] <- eval_results$TN
                cond_FP[ii] <- eval_results$FP
                cond_FN[ii] <- eval_results$FN
            },
            error = function(e) {
                message(sprintf("  Error processing species %d (%s): %s", species_idx, species_name, e$message))
                cond_auc[ii] <<- NA
                cond_threshold[ii] <<- NA
                cond_boyce[ii] <<- NA
                cond_sensitivity[ii] <<- NA
                cond_specificity[ii] <<- NA
                cond_TP[ii] <<- NA
                cond_TN[ii] <<- NA
                cond_FP[ii] <<- NA
                cond_FN[ii] <<- NA
            }
        )
    }
} else {
    # Fill with NA if no conditional predictions
    cond_auc[] <- NA
    cond_boyce[] <- NA
    cond_threshold[] <- NA
    cond_sensitivity[] <- NA
    cond_specificity[] <- NA
    cond_TP[] <- NA
    cond_TN[] <- NA
    cond_FP[] <- NA
    cond_FN[] <- NA
}

# ============================================================
# COMPILE AND SAVE RESULTS
# ============================================================

# ---- Model results ----
results_model_df <- data.frame(
    EXP_ROOT = rep(EXP_ROOT, n_species),
    EXP_ID = rep(EXP_ID, n_species),
    CLUSTER = rep(CLUSTER, n_species),
    REPNO = rep(REPNO, n_species),
    eval_type = rep("model", n_species),
    species_idx = species_indices,
    species_name = species_names[species_indices],
    AUC = model_auc,
    Boyce = model_boyce,
    threshold = model_threshold,
    sensitivity = as.numeric(model_sensitivity),
    specificity = as.numeric(model_specificity),
    TP = as.numeric(model_TP),
    TN = as.numeric(model_TN),
    FP = as.numeric(model_FP),
    FN = as.numeric(model_FN)
)

# ---- Conditional results ----
results_cond_df <- data.frame(
    EXP_ROOT = rep(EXP_ROOT, n_species),
    EXP_ID = rep(EXP_ID, n_species),
    CLUSTER = rep(CLUSTER, n_species),
    REPNO = rep(REPNO, n_species),
    eval_type = rep("conditional", n_species),
    species_idx = species_indices,
    species_name = species_names[species_indices],
    AUC = cond_auc,
    Boyce = cond_boyce,
    threshold = cond_threshold,
    sensitivity = as.numeric(cond_sensitivity),
    specificity = as.numeric(cond_specificity),
    TP = as.numeric(cond_TP),
    TN = as.numeric(cond_TN),
    FP = as.numeric(cond_FP),
    FN = as.numeric(cond_FN)
)

# ---- Combined results ----
results_df <- rbind(results_model_df, results_cond_df)

# ---- Save separate files ----
eval_model_filename <- paste0(base_prefix, "_eval_model.Rdata")
eval_cond_filename <- paste0(base_prefix, "_eval_cond.Rdata")
eval_combined_filename <- paste0(base_prefix, "_eval_combined.Rdata")

save(results_model_df, file = file.path(eval_dir, eval_model_filename))
message(sprintf("Saved model evaluation results: %s", eval_model_filename))

if (has_cond_pred) {
    save(results_cond_df, file = file.path(eval_dir, eval_cond_filename))
    message(sprintf("Saved conditional evaluation results: %s", eval_cond_filename))
}

save(results_df, file = file.path(eval_dir, eval_combined_filename))
message(sprintf("Saved combined evaluation results: %s", eval_combined_filename))

# ---- Print summary ----
message("\n========== EVALUATION SUMMARY ==========")

message("\n--- Model-based Evaluation ---")
message(sprintf("Species evaluated: %d", n_species))
message(sprintf("Mean AUC: %.3f (SD: %.3f)", mean(results_model_df$AUC, na.rm = TRUE), sd(results_model_df$AUC, na.rm = TRUE)))
message(sprintf("AUC range: [%.3f, %.3f]", min(results_model_df$AUC, na.rm = TRUE), max(results_model_df$AUC, na.rm = TRUE)))

if (has_cond_pred) {
    message("\n--- Conditional Prediction Evaluation ---")
    message(sprintf("Species evaluated: %d", n_species))
    message(sprintf("Mean AUC: %.3f (SD: %.3f)", mean(results_cond_df$AUC, na.rm = TRUE), sd(results_cond_df$AUC, na.rm = TRUE)))
    message(sprintf("AUC range: [%.3f, %.3f]", min(results_cond_df$AUC, na.rm = TRUE), max(results_cond_df$AUC, na.rm = TRUE)))

    message("\n--- Comparison ---")
    auc_diff <- results_cond_df$AUC - results_model_df$AUC
    message(sprintf("Mean AUC improvement (cond - model): %.3f (SD: %.3f)", mean(auc_diff, na.rm = TRUE), sd(auc_diff, na.rm = TRUE)))
    message(sprintf("Species improved: %d / %d", sum(auc_diff > 0, na.rm = TRUE), sum(!is.na(auc_diff))))
}

# ---- Create best model selection CSV ----
message("\n--- Creating best model selection file ---")

# Get all species from training data (use species_names which has all trained species)
all_species <- test_data$species_names

# Create a lookup for model AUC by species name
model_auc_lookup <- setNames(results_model_df$AUC, results_model_df$species_name)
cond_auc_lookup <- setNames(results_cond_df$AUC, results_cond_df$species_name)

# Determine best model for each species
best_model_df <- data.frame(
    species_name = all_species,
    best_model = character(length(all_species)),
    model_AUC = numeric(length(all_species)),
    cond_AUC = numeric(length(all_species)),
    stringsAsFactors = FALSE
)

for (i in seq_along(all_species)) {
    sp <- all_species[i]

    model_auc <- model_auc_lookup[sp]
    cond_auc <- cond_auc_lookup[sp]

    best_model_df$model_AUC[i] <- ifelse(is.null(model_auc), NA, model_auc)
    best_model_df$cond_AUC[i] <- ifelse(is.null(cond_auc), NA, cond_auc)

    # Determine best model
    # Default to "model" if either AUC is NA/NULL or if they're equal
    if (is.na(model_auc) || is.na(cond_auc)) {
        best_model_df$best_model[i] <- "model"
    } else if (cond_auc > model_auc) {
        best_model_df$best_model[i] <- "cond_pred"
    } else {
        best_model_df$best_model[i] <- "model"
    }
}

# Summary stats
n_model_best <- sum(best_model_df$best_model == "model")
n_cond_best <- sum(best_model_df$best_model == "cond_pred")
n_na <- sum(is.na(best_model_df$model_AUC) | is.na(best_model_df$cond_AUC))

message(sprintf(
    "Best model selection: %d model, %d cond_pred, %d with NA (defaulted to model)",
    n_model_best - n_na, n_cond_best, n_na
))

# Save to CSV
best_model_filename <- paste0(base_prefix, "_best_model.csv")
write.csv(best_model_df, file = file.path(eval_dir, best_model_filename), row.names = FALSE)
message(sprintf("Saved best model selection: %s", best_model_filename))

message(sprintf("\nDone! Evaluation complete for cluster %s, rep %d", CLUSTER, REPNO))
