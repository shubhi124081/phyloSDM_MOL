#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
    using namespace density;

    // ------------------ DATA ------------------
    DATA_INTEGER(N);
    DATA_INTEGER(J);
    DATA_INTEGER(K);
    DATA_IVECTOR(species); // 1-based indices from R
    DATA_MATRIX(X);        // N x K
    DATA_IVECTOR(y);       // length N (0/1 or counts)
    DATA_MATRIX(D_phylo);  // J x J
    DATA_VECTOR(offset);   // length N (log-offset)

    // ---------------- PARAMETERS ---------------
    PARAMETER_MATRIX(B); // J x K
    PARAMETER(log_alpha);
    PARAMETER(log_rho);
    PARAMETER(log_sigma_f);
    PARAMETER_VECTOR(f); // length N

    // Positive transforms
    Type alpha = exp(log_alpha);
    Type rho = exp(log_rho);
    Type sigma_f = exp(log_sigma_f);

    Type nll = 0.0;

    // ------------------ PRIORS -----------------
    // Match Stan: alpha, rho, sigma_f ~ Normal(0,1) constrained positive (half-normal)
    // Because we parameterize with logs, include Jacobian terms.
    nll -= dnorm(alpha, Type(0.0), Type(1.0), true) + log_alpha;
    nll -= dnorm(rho, Type(0.0), Type(1.0), true) + log_rho;
    nll -= dnorm(sigma_f, Type(0.0), Type(1.0), true) + log_sigma_f;

    // ------------------ GP prior on B ------------------
    matrix<Type> Kphy(J, J);
    for (int i = 0; i < J; ++i)
    {
        for (int j = 0; j < J; ++j)
        {
            Type d = D_phylo(i, j);
            Kphy(i, j) = alpha * alpha * exp(-(d * d) / (Type(2.0) * rho * rho));
        }
    }
    for (int j = 0; j < J; ++j)
        Kphy(j, j) += Type(1e-6);

    MVNORM_t<Type> mvn(Kphy);
    for (int k = 0; k < K; ++k)
    {
        nll += mvn(B.col(k));
    }

    // ------------------ Random effect prior ------------------
    for (int n = 0; n < N; ++n)
    {
        nll -= dnorm(f(n), Type(0.0), sigma_f, true);
    }

    // ------------------ Likelihood ------------------
    for (int n = 0; n < N; ++n)
    {
        int s = species(n) - 1; // convert to 0-based for C++ indexing

        Type xb = Type(0.0);
        for (int k = 0; k < K; ++k)
        {
            xb += X(n, k) * B(s, k);
        }

        Type log_lambda = xb + f(n) + offset(n);
        nll -= dpois(Type(y(n)), exp(log_lambda), true);
    }

    ADREPORT(alpha);
    ADREPORT(rho);
    ADREPORT(sigma_f);

    return nll;
}
