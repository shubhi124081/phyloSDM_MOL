LGCP_background <- "
  data {
    int<lower=1> N;               // Total number of observations (presences + background)
    int<lower=1> J;               // Number of species
    int<lower=1> K;               // Number of environmental covariates
    int<lower=1> species[N];      // Species index for each observation
    matrix[N, K] X;               // Environmental covariate matrix
    int<lower=0, upper=1> y[N];   // Presence (1) or background (0)
    matrix[J, J] D_phylo;         // Phylogenetic distance matrix
    vector[N] offset;             // Log effort or area per observation (e.g., log(duration × observers × distance))
}

parameters {
    matrix[J, K] B;               // Species-environment regression coefficients
    real<lower=0> alpha;          // GP amplitude (variance parameter)
    real<lower=0> rho;            // GP length scale (phylogenetic distance)
    real<lower=0> sigma_f;        // SD for spatial random effects
    vector[N] f;                  // Spatial random effects (site-specific)
}

transformed parameters {
    matrix[J, J] K_phylo;         // Phylogenetic covariance matrix
    matrix[J, J] L_B;             // Cholesky decomposition of K_phylo

    // Squared exponential kernel
    for (i in 1:J) {
        for (j in 1:J) {
            K_phylo[i, j] = alpha^2 * exp(-D_phylo[i, j]^2 / (2 * rho^2));
        }
    }

    // Add jitter for numerical stability
    K_phylo = K_phylo + diag_matrix(rep_vector(1e-6, J));
    L_B = cholesky_decompose(K_phylo);
}

model {
    // Priors
    alpha ~ normal(0, 1);
    rho ~ normal(0, 1);
    sigma_f ~ normal(0, 1);

    // GP prior on each column of B (K environmental covariates)
    for (k in 1:K) {
        B[, k] ~ multi_normal_cholesky(rep_vector(0, J), L_B);
    }

    // Spatial random effect
    f ~ normal(0, sigma_f);

    // Poisson likelihood with offset and log-link (presence = 1, background = 0)
    for (n in 1:N) {
        int s = species[n];
        real log_lambda = dot_product(X[n], B[s, ]) + f[n] + offset[n];
        y[n] ~ poisson_log(log_lambda);
    }
}
"
