functions {
  // Your existing helper on include path
  #include component_product.stan
}

data {
  int<lower=1> T;                       // time points
  int<lower=1> C;                       // ALR dimension (K_comps - 1)
  array[T] vector[C] Y;                 // responses on ALR scale

  int<lower=0> N;                       // # covariates across all components
  array[C] int<lower=0> K;              // covariate counts by component
  array[T] vector[N] X;                 // covariates by time

  int<lower=1> P;                       // VAR order (>=1)
  int<lower=0> Q;                       // VMA order (can be 0)

  int<lower=0> T_new;                   // forecast steps
  array[T_new] vector[N] X_new;         // future covariates

  int prior_only;                       // 1 = ignore likelihood
  real<lower=0> beta_sd;                // prior sd for beta
  real<lower=0> rw_scale_sd;            // half-N prior sd for innovation scales (A,B)
}

transformed data {
  int M = max(P, Q);                    // warm-up length for recursion
}

parameters {
  // Exogenous effects (kept time-invariant)
  vector[N] beta;

  // Long-run (mean-reversion targets) for VAR and VMA
  array[P] matrix[C, C] A_bar;
  array[Q] matrix[C, C] B_bar;

  // Initial states at t = 1
  array[P] matrix[C, C] A1;
  array[Q] matrix[C, C] B1;

  // Non-centered innovations for AR(1) evolutions:
  // z_A[t-1,p] ~ N(0, I), z_B[t-1,q] ~ N(0, I)
  array[T - 1, P] matrix[C, C] z_A;
  array[T - 1, Q] matrix[C, C] z_B;

  // Mean-reversion (AR(1)) parameters
  vector<lower=0, upper=1>[P] rho_A;    // persistence for A
  vector<lower=0>[P]          tau_A;    // innovation scale for A
  vector<lower=0, upper=1>[Q] rho_B;    // persistence for B
  vector<lower=0>[Q]          tau_B;    // innovation scale for B

  // Homoscedastic observation covariance on ALR scale
  vector<lower=0>[C] sigma;
  cholesky_factor_corr[C] L_Omega;
}

transformed parameters {
  // Covariate contribution per component (unchanged)
  array[T] vector[C] Xbeta = component_product(X, beta, K);

  // Cholesky of Sigma
  matrix[C, C] L_Sigma = diag_pre_multiply(sigma, L_Omega);

  // Full time paths with mean reversion to A_bar/B_bar
  array[T, P] matrix[C, C] A_t;
  for (p in 1:P) A_t[1, p] = A1[p];
  for (t in 2:T) {
    for (p in 1:P) {
      A_t[t, p] = A_bar[p]
                + rho_A[p] * (A_t[t - 1, p] - A_bar[p])
                + tau_A[p] * z_A[t - 1, p];
    }
  }

  array[T, Q] matrix[C, C] B_t;
  if (Q > 0) {
    for (q in 1:Q) B_t[1, q] = B1[q];
    for (t in 2:T) {
      for (q in 1:Q) {
        B_t[t, q] = B_bar[q]
                  + rho_B[q] * (B_t[t - 1, q] - B_bar[q])
                  + tau_B[q] * z_B[t - 1, q];
      }
    }
  }

  // Conditional means via VARMA recursion
  array[T] vector[C] eta;
  for (t in 1:M) eta[t] = Y[t];  // warm start
  for (t in (M + 1):T) {
    vector[C] ar = rep_vector(0, C);
    vector[C] ma = rep_vector(0, C);

    for (p in 1:P)
      ar += A_t[t, p] * (Y[t - p] - Xbeta[t - p]);

    if (Q > 0)
      for (q in 1:Q)
        ma += B_t[t, q] * (Y[t - q] - eta[t - q]);

    eta[t] = ar + ma + Xbeta[t];
  }
}

model {
  // Priors
  beta ~ normal(0, beta_sd);

  // Long-run matrices: diagonals near modest persistence, off-diagonals shrink to 0
  for (p in 1:P) {
    diagonal(A_bar[p]) ~ normal(0.5, 0.3);
    for (i in 1:C) for (j in 1:C) if (i != j) A_bar[p][i, j] ~ normal(0, 0.2);
  }
  if (Q > 0) {
    for (q in 1:Q) {
      diagonal(B_bar[q]) ~ normal(0.5, 0.3);
      for (i in 1:C) for (j in 1:C) if (i != j) B_bar[q][i, j] ~ normal(0, 0.2);
    }
  }

  // Initial states shrink toward long-run means for stability
  for (p in 1:P) to_vector(A1[p] - A_bar[p]) ~ normal(0, 0.2);
  if (Q > 0) for (q in 1:Q) to_vector(B1[q] - B_bar[q]) ~ normal(0, 0.2);

  // Mean-reversion strength and innovation scales
  rho_A ~ beta(9, 1);                 // mean ≈ 0.90
  tau_A ~ normal(0, rw_scale_sd);     // half-N via <lower=0>

  if (Q > 0) {
    rho_B ~ beta(9, 1);
    tau_B ~ normal(0, rw_scale_sd);
  }

  // Non-centered innovations
  for (t in 1:(T - 1)) for (p in 1:P) to_vector(z_A[t, p]) ~ normal(0, 1);
  if (Q > 0) for (t in 1:(T - 1)) for (q in 1:Q) to_vector(z_B[t, q]) ~ normal(0, 1);

  // Homoscedastic Sigma
  L_Omega ~ lkj_corr_cholesky(3);
  sigma   ~ normal(0, 0.5);

  // Likelihood
  if (!prior_only) {
    for (t in (M + 1):T)
      Y[t] ~ multi_normal_cholesky(eta[t], L_Sigma);
  }
}

generated quantities {
  // Pointwise log-lik
  vector[T - M] log_lik;
  {
    matrix[C, C] L_S = diag_pre_multiply(sigma, L_Omega);
    for (t in (M + 1):T)
      log_lik[t - M] = multi_normal_cholesky_lpdf(Y[t] | eta[t], L_S);
  }

  // In-sample one-step predictions (training window)
  array[T] vector[C] Y_hat_train;
  {
    matrix[C, C] L_S = diag_pre_multiply(sigma, L_Omega);
    for (t in 1:T) {
      vector[C] mu;
      if (t <= M) {
        mu = Y[t];
      } else {
        vector[C] ar = rep_vector(0, C);
        vector[C] ma = rep_vector(0, C);

        for (p in 1:P)
          ar += A_t[t, p] * (Y[t - p] - Xbeta[t - p]);

        if (Q > 0)
          for (q in 1:Q)
            ma += B_t[t, q] * (Y[t - q] - eta[t - q]);

        mu = ar + ma + Xbeta[t];
      }
      Y_hat_train[t] = multi_normal_cholesky_rng(mu, L_S);
    }
  }

  // Out-of-sample forecasts with deterministic mean reversion (no lock-in)
  array[T_new] vector[C] Y_hat;
  if (T_new > 0) {
    array[T_new] vector[C] Xbeta_new = component_product(X_new, beta, K);
    matrix[C, C] L_S = diag_pre_multiply(sigma, L_Omega);

    // Start from last in-sample coefficient states and revert each step
    array[P] matrix[C, C] A_for;
    for (p in 1:P) A_for[p] = A_t[T, p];

    array[Q] matrix[C, C] B_for;
    if (Q > 0) for (q in 1:Q) B_for[q] = B_t[T, q];

    array[T_new] vector[C] eta_new;

    for (h in 1:T_new) {
      // Deterministic reversion of A and B toward their long-run means
      for (p in 1:P)
        A_for[p] = A_bar[p] + rho_A[p] * (A_for[p] - A_bar[p]);
      if (Q > 0)
        for (q in 1:Q)
          B_for[q] = B_bar[q] + rho_B[q] * (B_for[q] - B_bar[q]);

      vector[C] ar = rep_vector(0, C);
      vector[C] ma = rep_vector(0, C);

      for (p in 1:P) {
        vector[C] Y_lag;
        vector[C] Xbeta_lag;
        if (h < p + 1) {
          Y_lag     = Y[T + h - p];
          Xbeta_lag = Xbeta[T + h - p];
        } else {
          Y_lag     = Y_hat[h - p];
          Xbeta_lag = Xbeta_new[h - p];
        }
        ar += A_for[p] * (Y_lag - Xbeta_lag);
      }

      if (Q > 0) {
        for (q in 1:Q) {
          vector[C] Y_lag;
          vector[C] eta_lag;
          if (h < q + 1) {
            Y_lag  = Y[T + h - q];
            eta_lag = eta[T + h - q];
          } else {
            Y_lag  = Y_hat[h - q];
            eta_lag = eta_new[h - q];
          }
          ma += B_for[q] * (Y_lag - eta_lag);
        }
      }

      eta_new[h] = ar + ma + Xbeta_new[h];
      Y_hat[h]   = multi_normal_cholesky_rng(eta_new[h], L_S);
    }
  }
}
