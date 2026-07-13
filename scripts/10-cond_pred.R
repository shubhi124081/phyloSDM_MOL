# Purpose: Generate conditional predictions for species using fitted model
# This script loads fitted Stan models and generates predictions for held-out species
# Input: Command line args or interactive settings
# Output: Conditional predictions saved to analysis/<EXP_ROOT>/cond_pred/

# Clear workspace
rm(list = ls())

# ---- Libraries ----
library(rstan)
library(mvtnorm)

# ============================================================
# HELPER FUNCTION: Conditional predictor
# ============================================================

# Simple, stable per-draw conditional predictor (species-level GP over betas)
# Returns: data.frame with S rows (posterior draws) and K columns (predictors)
# Notes:
# - Assumes new_index is a single focal species.
# - Variance is species-level (same sd for all predictors), which matches your model.
# - Optional nugget on K_RR if you have species-specific st_devs draws.

cp_LGCP_simple <- function(
    B_draws, # [S, J, K] posterior draws of species-environment coefficients
    alpha_draws, # [S] posterior draws of GP amplitude
    rho_draws, # [S] posterior draws of GP length scale
    D_phylo, # J x J distance matrix among species
    observed_index, # integer vector of observed-species indices (R)
    new_index, # single integer index of focal species (D)
    stan_data, # to pick K and column names
    use_draws = TRUE, # TRUE: propagate uncertainty; FALSE: plug-in means (S=1)
    st_devs_draws = NULL, # optional S x J matrix of species SDs; if provided, adds nugget
    jitter = 1e-6 # numeric jitter for numerical stability
    ) {
    stopifnot(length(new_index) == 1)

    S <- dim(B_draws)[1]
    J <- dim(B_draws)[2]
    K <- dim(B_draws)[3]

    # Optional plug-in (fast, less correct): collapse to means
    if (!use_draws) {
        B_draws <- array(apply(B_draws, c(2, 3), mean), dim = c(1, J, K))
        alpha_draws <- mean(alpha_draws)
        rho_draws <- mean(rho_draws)
        if (!is.null(st_devs_draws)) {
            # collapse st_devs if provided as draws
            st_devs_draws <- matrix(colMeans(st_devs_draws), nrow = 1) # [1, J]
        }
        S <- 1
    }

    # Kernel function (squared exponential on phylo distances)
    ker <- function(a, r, D) a^2 * exp(-(D^2) / (2 * r^2))

    out <- matrix(NA_real_, nrow = S, ncol = K)

    for (s in seq_len(S)) {
        # Build kernel for this draw
        Kfull <- ker(alpha_draws[s], rho_draws[s], D_phylo)

        # Partition by species
        K_RR <- Kfull[observed_index, observed_index, drop = FALSE]
        K_DR <- Kfull[new_index, observed_index, drop = FALSE] # 1 x |R|
        K_DD <- Kfull[new_index, new_index, drop = FALSE] # scalar

        # Optional nugget on observed block (species-specific residual sd)
        if (!is.null(st_devs_draws)) {
            stopifnot(ncol(st_devs_draws) == J)
            K_RR <- K_RR + diag(st_devs_draws[s, observed_index]^2, nrow = length(observed_index))
        }

        # Stable solve via Cholesky + jitter
        L_RR <- chol(K_RR + diag(jitter, nrow(K_RR)))
        solve_KRR <- function(B) backsolve(L_RR, forwardsolve(t(L_RR), B))

        # Collect observed betas for ALL predictors at once: (|R| x K)
        B_R_mat <- do.call(cbind, lapply(seq_len(K), function(k) B_draws[s, observed_index, k]))

        # Conditional mean for ALL predictors: (1 x |R|) %*% (|R| x K) = (1 x K)
        mu_vec <- as.vector(K_DR %*% solve_KRR(B_R_mat))

        # Conditional variance is species-level (same for each predictor)
        var_D <- as.numeric(K_DD - K_DR %*% solve_KRR(t(K_DR)))
        sd_D <- sqrt(max(var_D, jitter))

        # One posterior-predictive draw per predictor for this draw s
        out[s, ] <- rnorm(K, mean = mu_vec, sd = sd_D)
    }

    df <- as.data.frame(out)
    colnames(df) <- colnames(stan_data$X)
    df
}

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

# ---- Command line arguments or manual set ----
# ---- Command line arguments or manual set ----
if (interactive()) {
    EXP_ROOT <- "v0"
    EXP_ID <- "full_test"
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
message(sprintf("Conditional prediction for cluster: %s, rep: %d", CLUSTER, REPNO))

# ---- Paths ----
dpath <- file.path(raw_dir, DATASET)
data_dir <- file.path(dpath, CLUSTER)

# Create output directory
cond_pred_dir <- file.path(analysis_dir, EXP_ROOT, "cond_pred")
if (!dir.exists(cond_pred_dir)) {
    dir.create(cond_pred_dir, recursive = TRUE)
    message(sprintf("Created output directory: %s", cond_pred_dir))
}

# ---- Load result file ----
# Result file naming: {EXP_ROOT}_{EXP_ID}_{CLUSTER}_{FSP}_rep_{REPNO}_{MODEL_TYPE}_fit.Rdata
result_file <- file.path(res_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_", MODEL_TYPE, "_fit.Rdata"))
if (!file.exists(result_file)) {
    stop(sprintf("Result file not found: %s", result_file))
}
load(result_file) # loads 'result'
message(sprintf("Loaded result file: %s", basename(result_file)))

# result$posterior (local/cmdstanr fits) is already extracted plain arrays;
# result$fit (HPC/rstan fits) is a live stanfit object needing rstan::extract().
if (!is.null(result$posterior)) {
    B_draws <- result$posterior$B
    alpha_draws <- result$posterior$alpha
    rho_draws <- result$posterior$rho
} else {
    stan_fit <- result$fit
    B_draws <- rstan::extract(stan_fit, pars = "B")$B
    alpha_draws <- rstan::extract(stan_fit, pars = "alpha")$alpha
    rho_draws <- rstan::extract(stan_fit, pars = "rho")$rho
}

# ---- Load model data file ----
# Model data file naming: {EXP_ROOT}_{EXP_ID}_{CLUSTER}_{FSP}_rep_{REPNO}_model_data.Rdata
# Okay - accidentally saved everything to raw_data instead of data so loading from there right now
model_data_file <- file.path(data_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_model_data.Rdata"))
if (!file.exists(model_data_file)) {
    stop(sprintf("Model data file not found: %s", model_data_file))
}
load(model_data_file) # loads 'model_data'
message(sprintf("Loaded model data: N=%d, J=%d, K=%d", model_data$N, model_data$J, model_data$K))

# ---- Prepare data for conditional prediction ----
# Some dimensions we may need
K <- model_data$K
J <- model_data$J
N_obs <- model_data$N_obs
N <- model_data$N

# Data we need
D_phylo <- model_data$D_phylo # phylogenetic distance matrix
sps <- model_data$species_names # species names

# Prepare stan_data structure for cp_LGCP_simple function
stan_data <- list(
    X = model_data$X,
    K = K
)

# ---- Run conditional prediction ----
message("Starting conditional prediction...")
B_cond_list <- list()

for (new_index in seq_len(length(sps))) {
    observed_index <- seq_len(J)[-new_index] # Indices of species with training data
    sp_name <- sps[new_index]
    message(sprintf("  Predicting for species %d/%d: %s", new_index, length(sps), sp_name))

    # Call the function
    B_cond_list[[new_index]] <- cp_LGCP_simple(B_draws, alpha_draws, rho_draws, D_phylo, observed_index, new_index, stan_data)
}
names(B_cond_list) <- sps

# ---- Save the results ----
out_file <- file.path(cond_pred_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_cond_pred.Rdata"))
save(B_cond_list, file = out_file)
message(sprintf("Done! Conditional predictions saved to: %s", out_file))
