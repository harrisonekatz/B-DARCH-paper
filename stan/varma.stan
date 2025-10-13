functions {
  #include component_product.stan
}

data {
  int<lower=1> T;                       // number of time periods
  int<lower=1> C;                       // number of categories
  array[T] vector[C] Y;                 // response array

  int<lower=0> N;                       // total number of covariates
  array[C] int<lower=0> K;              // number of covariates by component
  array[T] vector[N] X;                 // covariates across all components

  int<lower=0> P;                       // number of auto-regressive lags
  int<lower=0> Q;                       // number of moving average lags

  int<lower=0> T_new;                   // number of time periods to forecast
  array[T_new] vector[N] X_new;         // new covariates across all components

  int prior_only;                       // whether to ignore the likelihood (1)
  real<lower=0> beta_sd;                // SD of the prior for beta
}

transformed data {
  int M = max(P, Q);
}

parameters {
  vector[N] beta;                   // coefficients for covariates
  array[P] matrix[C, C] A;          // VAR coefficients
  array[Q] matrix[C, C] B;          // VMA coefficients

  vector<lower=0>[C] sigma;         // SDs of white noise process
  cholesky_factor_corr[C] L_Omega;  // white noise correlation matrix
}

transformed parameters {
  array[T] vector[C] Xbeta = component_product(X, beta, K);

  array[T] vector[C] eta;

  matrix[C, C] L_Sigma = diag_pre_multiply(sigma, L_Omega);

  eta[1:M] = Y[1:M];

  for (t in (M + 1):T) {
    vector[C] ar = rep_vector(0, C);
    vector[C] ma = rep_vector(0, C);

    for (p in 1:P) {
      ar += A[p] * (Y[t - p] - Xbeta[t - p]);
    }

    for (q in 1:Q) {
      ma += B[q] * (Y[t - q] - eta[t - q]);
    }

    eta[t] = ar + ma + Xbeta[t];
  }
}

model {
  // Prior model
  beta ~ normal(0, beta_sd);

  for (p in 1:P) {
    diagonal(A[p]) ~ normal(0, 1);

    for (i in 1:C) {
      for (j in 1:C) {
        if (i != j) {
          A[p, i, j] ~ normal(0, 1);
        }
      }
    }
  }

  for (q in 1:Q) {
    diagonal(B[q]) ~ normal(0, 1);

    for (i in 1:C) {
      for (j in 1:C) {
        if (i != j) {
          B[q, i, j] ~ normal(0, 1);
        }
      }
    }
  }

  L_Omega ~ lkj_corr_cholesky(3);
  sigma ~ normal(0, 0.5);

  // Observational model
  if (!prior_only) {
    Y[(M + 1):T] ~ multi_normal_cholesky(eta[(M + 1):T], L_Sigma);
  }
}

generated quantities {
  array[T] vector[C] Y_hat_train;  // Fitted values for training data
  array[T_new] vector[C] Y_hat;    // Predictions for new data
  vector[T - M] log_lik;

  matrix[C, C] Omega = multiply_lower_tri_self_transpose(L_Omega);
  matrix[C, C] Sigma = quad_form_diag(Omega, sigma);

  // Log-likelihood contributions
  for (t in (M + 1):T) {
    log_lik[t - M] = multi_normal_cholesky_lpdf(Y[t] | eta[t], L_Sigma);
  }

  // Generate fitted values for training data
  for (t in 1:T) {
    Y_hat_train[t] = multi_normal_cholesky_rng(eta[t], L_Sigma);
  }

  // Predictions for out-of-sample data
  if (T_new > 0) {
    array[T_new] vector[C] Xbeta_new = component_product(X_new, beta, K);
    array[T_new] vector[C] eta_new;

    for (t in 1:T_new) {
      vector[C] ar = rep_vector(0, C);
      vector[C] ma = rep_vector(0, C);

      for (p in 1:P) {
        vector[C] Y_lag;
        vector[C] Xbeta_lag;

        if (t < p + 1) {
          Y_lag = Y[T + t - p];
          Xbeta_lag = Xbeta[T + t - p];
        } else {
          Y_lag = Y_hat[t - p];
          Xbeta_lag = Xbeta_new[t - p];
        }

        ar += A[p] * (Y_lag - Xbeta_lag);
      }

      for (q in 1:Q) {
        vector[C] Y_lag;
        vector[C] eta_lag;

        if (t < q + 1) {
          Y_lag = Y[T + t - q];
          eta_lag = eta[T + t - q];
        } else {
          Y_lag = Y_hat[t - q];
          eta_lag = eta_new[t - q];
        }

        ma += B[q] * (Y_lag - eta_lag);
      }

      eta_new[t] = ar + ma + Xbeta_new[t];
      Y_hat[t] = multi_normal_cholesky_rng(eta_new[t], L_Sigma);
    }
  }
}
