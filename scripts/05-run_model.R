# Purpose: Run species distribution model (Stan or TMB)
# This script loads packaged data and fits the specified model
# Input: Command line args or interactive settings
# Output: Model fit saved to res/<EXP_ID>/

# Clear workspace
rm(list = ls())

# ---- Libraries ----
library(yaml)

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
    message("Running on HPC")
} else {
    root <- "~/clim_risk_phylosdm"
    message("Running locally")
}

scripts_directory <- file.path(root, "scripts")
res_dir <- file.path(root, "res")

# ---- Command line arguments or manual set ----
if (interactive()) {
    EXP_ROOT <- "v0"
    EXP_ID <- "sub1000"
    DATASET <- "amphibians"
    CLUSTER <- "Rani1"
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
message(sprintf("Model type: %s, Model name: %s", MODEL_TYPE, MODEL_NAME))

# ---- Paths ----
dpath <- file.path(root, "raw_data", DATASET)
data_dir <- file.path(dpath, CLUSTER)

# ---- Load model_data ----

model_data_file <- file.path(data_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_model_data.Rdata"))
if (!file.exists(model_data_file)) {
    stop(sprintf("Model data file not found: %s. Run 04-pkg_data.R first.
    \n Note: Exp. root and exp. ID need to be the same between packaging data and model runs", model_data_file))
}
load(model_data_file) # loads 'model_data'
message(sprintf(
    "Loaded data: N=%d, J=%d, K=%d, N_obs=%d",
    model_data$N, model_data$J, model_data$K, model_data$N_obs
))

# ---- Synchronize species indices with phylogeny ----
# Some species may have been dropped during data processing (e.g., missing env data)
# Need to ensure D_phylo, species_names, and species indices are consistent
observed_species_idx <- sort(unique(model_data$species))
D_phylo_species <- seq_len(nrow(model_data$D_phylo))
n_observed <- length(observed_species_idx)


if (!identical(as.integer(observed_species_idx), D_phylo_species)) {
    message(sprintf("  WARNING: Only %d of %d species have observations. Synchronizing...", n_observed, model_data$J))

    # Subset species_names to only observed species
    model_data$species_names <- model_data$species_names[observed_species_idx]

    # Subset D_phylo to only observed species
    model_data$D_phylo <- model_data$D_phylo[observed_species_idx, observed_species_idx, drop = FALSE]

    # Create mapping from old indices to new contiguous indices (1:n_observed)
    old_to_new <- setNames(seq_along(observed_species_idx), observed_species_idx)
    model_data$species <- as.integer(old_to_new[as.character(model_data$species)])

    # Update J
    model_data$J <- n_observed

    message(sprintf(
        "  Synchronized: J=%d, D_phylo=%dx%d, species range=[%d,%d]",
        model_data$J, nrow(model_data$D_phylo), ncol(model_data$D_phylo),
        min(model_data$species), max(model_data$species)
    ))
}

# ---- Prepare data ----
# Ensure offset doesn't have zeros (for log)
if ("offset" %in% names(model_data)) {
    model_data$offset <- model_data$offset + 1e-10
}
model_specs <- list(
    "iter" = 8000,
    "warmup" = 3000,
    "chains" = 1,
    "thin" = 10,
    "cores" = 1
)

# ============================================================
# RUN MODEL
# ============================================================

if (MODEL_TYPE == "STAN") {
    # ---- Stan Model ----
    message("Running Stan model...")
    library(rstan)
    options(mc.cores = model_specs$cores)

    # Source stan model code
    source(file.path(scripts_directory, "05-stan_model.R"))

    # Get model code (model name without extension)
    model_code_name <- MODEL_NAME
    if (!exists(model_code_name)) {
        stop(sprintf("Stan model '%s' not found in 05-stan_model.R", model_code_name))
    }
    code <- get(model_code_name)

    # X is site-level (N sites x K), but TMB expects observation-level (N_obs x K)
    # Expand X using site indices
    X_obs <- as.matrix(model_data$X)[model_data$site, , drop = FALSE]


    # Scale the phylogenetic distance matrix
    D_phylo_max <- max(model_data$D_phylo)
    if (D_phylo_max > 0) {
        model_data$D_phylo <- model_data$D_phylo / D_phylo_max
        message("  Scaled phylogenetic distance matrix by max value")
    }

    # Prepare data for Stan (remove non-data elements)
    standata <- model_data[c("N", "J", "K", "species", "X", "y", "D_phylo", "species")]
    standata$N <- model_data$N_obs # Stan expects N = number of observations
    standata$X <- X_obs

    # Add offset if available
    if (!is.null(model_data$offset)) {
        standata$offset <- model_data$offset
        standata$use_offset <- 1L
        message("  Using offset in model")
    } else {
        standata$offset <- rep(0, model_data$N_obs)
        standata$use_offset <- 0L
        message("  No offset (disabled)")
    }

    # Fit model
    fit <- rstan::stan(
        model_code = code,
        data = standata,
        iter = model_specs$iter,
        thin = model_specs$thin,
        warmup = model_specs$warmup,
        chains = model_specs$chains,
        cores = model_specs$cores
    )

    # Package results
    result <- list(
        fit = fit,
        model_type = "STAN",
        model_name = MODEL_NAME,
        cluster = CLUSTER,
        repno = REPNO
        # config = config
    )
} else if (MODEL_TYPE == "TMB") {
    # ---- TMB Model ----
    message("Running TMB model...")
    library(TMB)

    # Compile TMB model if needed
    cpp_file <- file.path(scripts_directory, MODEL_NAME)
    if (!file.exists(cpp_file)) {
        stop(sprintf("TMB model file not found: %s", cpp_file))
    }

    # Get model name without .cpp extension
    model_base <- gsub("\\.cpp$", "", MODEL_NAME)
    dll_path <- file.path(scripts_directory, model_base)
    dll_file <- paste0(dll_path, if (.Platform$OS.type == "windows") ".dll" else ".so")

    # Compile if DLL doesn't exist or cpp is newer than DLL
    needs_compile <- !file.exists(dll_file) ||
        file.mtime(cpp_file) > file.mtime(dll_file)

    if (needs_compile) {
        message("Compiling TMB model...")
        # Remove old compiled files first
        suppressWarnings({
            try(file.remove(paste0(dll_path, ".o")), silent = TRUE)
            try(file.remove(dll_file), silent = TRUE)
        })
        old_wd <- getwd()
        setwd(scripts_directory)
        TMB::compile(MODEL_NAME)
        setwd(old_wd)
    }
    dyn.load(TMB::dynlib(dll_path))

    # Prepare data for TMB
    # X is site-level (N sites x K), but TMB expects observation-level (N_obs x K)
    # Expand X using site indices
    X_obs <- as.matrix(model_data$X)[model_data$site, , drop = FALSE]

    # Scale the phylogenetic distance matrix
    D_phylo_max <- max(model_data$D_phylo)
    if (D_phylo_max > 0) {
        model_data$D_phylo <- model_data$D_phylo / D_phylo_max
        message("  Scaled phylogenetic distance matrix by max value")
    }

    tmb_data <- list(
        N = model_data$N_obs,
        J = model_data$J,
        K = model_data$K,
        species = model_data$species, # Already 1-based from R
        X = X_obs,
        y = as.numeric(model_data$y),
        D_phylo = as.matrix(model_data$D_phylo),
        offset = as.numeric(model_data$offset)
        # offset = as.numeric(rep(0, length(model_data$offset)))
    )

    # Validate data bounds
    if (min(tmb_data$species) < 1 || max(tmb_data$species) > tmb_data$J) {
        stop(sprintf(
            "Species indices out of bounds! Range: [%d, %d], expected [1, %d]",
            min(tmb_data$species), max(tmb_data$species), tmb_data$J
        ))
    }
    if (nrow(tmb_data$X) != tmb_data$N) {
        stop(sprintf("X has %d rows but N=%d", nrow(tmb_data$X), tmb_data$N))
    }
    if (nrow(tmb_data$D_phylo) != tmb_data$J || ncol(tmb_data$D_phylo) != tmb_data$J) {
        stop(sprintf("D_phylo is %dx%d but J=%d", nrow(tmb_data$D_phylo), ncol(tmb_data$D_phylo), tmb_data$J))
    }
    message(sprintf(
        "  Data validated: species range [%d, %d], X: %dx%d, D_phylo: %dx%d",
        min(tmb_data$species), max(tmb_data$species),
        nrow(tmb_data$X), ncol(tmb_data$X),
        nrow(tmb_data$D_phylo), ncol(tmb_data$D_phylo)
    ))

    # Initial parameters (must match PARAMETER declarations in TMB model)
    tmb_par <- list(
        B = matrix(0, tmb_data$J, tmb_data$K),
        log_alpha = 0.5,
        log_rho = 1.0,
        log_sigma_f = 0.5,
        f = rep(0, tmb_data$N)
    )

    # Build objective function
    message("Building TMB objective function...")
    obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = tmb_par,
        DLL = model_base,
        random = c("f"),
        silent = TRUE
    )

    # Optimize
    message("Optimizing...")
    t_start <- Sys.time()
    opt <- nlminb(obj$par, obj$fn, obj$gr,
        control = list(iter.max = model_specs$iter, eval.max = model_specs$iter * 2)
    )
    t_end <- Sys.time()
    run_time <- as.numeric(difftime(t_end, t_start, units = "mins"))
    message(sprintf("Optimization complete in %.2f minutes", run_time))

    # Get standard errors
    message("Computing standard errors...")
    rep <- TMB::sdreport(obj)

    fix <- summary(rep, "fixed")
    b_rows <- grep("^B", rownames(fix))
    stopifnot(length(b_rows) == tmb_data$J * tmb_data$K)

    B_hat <- matrix(fix[b_rows, 1], nrow = tmb_data$J, ncol = tmb_data$K)
    B_se <- matrix(fix[b_rows, 2], nrow = tmb_data$J, ncol = tmb_data$K)


    # Extract hyperparameters
    hyper_summary <- summary(rep, "report")

    # Package results
    result <- list(
        opt = opt,
        rep = rep,
        B_hat = B_hat,
        B_se = B_se,
        hyperparameters = hyper_summary,
        run_time_mins = run_time,
        convergence = opt$convergence,
        model_type = "TMB",
        model_name = MODEL_NAME,
        cluster = CLUSTER,
        repno = REPNO,
        # config = config,
        species_names = model_data$species_names
    )
} else {
    stop(sprintf("Unknown model type: %s. Use 'STAN' or 'TMB'.", MODEL_TYPE))
}

# ---- Save results ----
out_file <- file.path(res_dir, paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO, "_", MODEL_TYPE, "_fit.Rdata"))
save(result, file = out_file)
message(sprintf("Done! Results saved to: %s", out_file))
