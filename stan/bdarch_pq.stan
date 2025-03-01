functions {
  #include alr.stan
  #include component_product.stan
}

data {
  int<lower=1> T;                       // number of time periods
  int<lower=2> C;                       // number of categories
  array[T] simplex[C] Y;                // response array
  int<lower=1, upper=C> ref;            // ALR reference element of simplex

  int<lower=0> N;                       // total number of covariates
  array[C - 1] int<lower=0> K;          // number of covariates by component
  array[T] vector[N] X;                 // covariates across all components

  int<lower=1> K_phi;                   // number of population effects for phi
  matrix[T, K_phi] X_phi;               // design matrix for phi

  int<lower=0> P;                       // number of auto-regressive lags
  int<lower=0> Q;                       // number of moving average lags

  int<lower=0> T_new;                   // number of time periods to forecast
  array[T_new] vector[N] X_new;         // new covariates across all components
  matrix[T_new, K_phi] X_phi_new;       // design matrix for forecasting phi

  int prior_only;                       // whether to ignore the likelihood (1)
  real<lower=0> beta_sd;                // SD of the prior for beta
}

transformed data {
  int M = max(P, Q);
  // length T so we can use it in generated quantities
  array[T] vector[C - 1] alr_Y;

  for (t in 1:T) {
    alr_Y[t] = alr(Y[t], ref);

  }
  // Array to hold initial values for phi
  vector[M] phi_init;

  // Generate uniform random initial values for phi

}


parameters {
  vector[N] beta;                   // coefficients for covariates
  array[P] matrix[C - 1, C - 1] A;  // VAR coefficients
  array[Q] matrix[C - 1, C - 1] B;  // VMA coefficients
  vector[K_phi] beta_phi;           // coefficients for phi covariates
  vector[P] ar_phi;                 // AR coefficients for phi
  vector[Q] ma_phi;                 // AR coefficients for phi
}

transformed parameters {
  vector[T] phi;
  array[T] vector[C - 1] Xbeta = component_product(X, beta, K);
  array[T] vector[C - 1] eta;
  array[T - M] vector[C] alpha;
  vector[T] Xbeta_phi = X_phi * beta_phi;

  // Compute phi with AR term

 for (m in 1:M) {
    eta[m] = alr_Y[m];
    phi[m] = Xbeta_phi[m];
  }


  for (t in (M + 1):T) {
    vector[C - 1] ar = rep_vector(0, C - 1);
    vector[C - 1] ma = rep_vector(0, C - 1);
    vector[1] phi_ar = rep_vector(0,1);
    vector[1] phi_ma = rep_vector(0,1);

    for (p in 1:P) {
      ar += A[p] * (alr_Y[t - p] - Xbeta[t - p]);
      phi_ar += ar_phi[p] * (phi[t - p]-Xbeta_phi[t-p]);
    }

    for (q in 1:Q) {
      ma += B[q] * (alr_Y[t - q] - eta[t - q]);
      phi_ma += ma_phi[q] * sum(square(alr_Y[t - q] - eta[t - q]));
    }

    eta[t] = ar + ma + Xbeta[t];
    phi[t] = phi_ar[1] + phi_ma[1] + Xbeta_phi[t];
  }

  for (t in (M + 1):T) {
    alpha[t - M] = exp(phi[t]) * alrinv(eta[t], ref);
  }
}

model {
  // prior model
  beta ~ normal(0, beta_sd);

  for (p in 1:P) {
    // different prior for diagonal elements
    diagonal(A[p]) ~ normal(0, 1);

    for (i in 1:(C - 1)) {
      for (j in 1:(C - 1)) {
        if (i != j) {
          A[p, i, j] ~ normal(0, 1);
        }
      }
    }
  }

  for (q in 1:Q) {
    diagonal(B[q]) ~ normal(0, 1);

    for (i in 1:(C - 1)) {
      for (j in 1:(C - 1)) {
        if (i != j) {
          B[q, i, j] ~ normal(0, 1);
        }
      }
    }
  }

  // different prior for the intercept
  beta_phi[1] ~ normal(7, 1.5);
  if (num_elements(beta_phi) > 1) {
    beta_phi[2:] ~ normal(0, .1);
  }

  ar_phi ~ normal(.5, .25); // prior for AR coefficients for phi
  ma_phi ~ normal(-.75, .5);
  // observational model
  if (!prior_only) {
    Y[(M + 1):T] ~ dirichlet(alpha);
  }
}

generated quantities {
  array[T_new] simplex[C] Y_hat;
  vector[T - M] log_lik;

  for (t in (M + 1):T) {
    log_lik[t - M] = dirichlet_lpdf(Y[t] | alpha[t - M]);
  }

  if (T_new > 0) {
    {
      array[T_new] vector[C - 1] Xbeta_new = component_product(X_new, beta, K);
      vector[T_new] Xbeta_phi_new = X_phi_new * beta_phi;
      array[T_new] vector[C - 1] alr_Y_hat;
      array[T_new] vector[C - 1] eta_new;
      vector[T_new] phi_new;
      array[T_new] vector[C] alpha_new;

      for (t in 1:T_new) {
        vector[C - 1] alr_Y_lag;
        vector[C - 1] Xbeta_lag;
        vector[1] Xbeta_phi_lag;
        real phi_lag;
        vector[C - 1] eta_lag;
        vector[C - 1] ar = rep_vector(0, C - 1);
        vector[C - 1] ma = rep_vector(0, C - 1);
        vector[1] phi_ar = rep_vector(0,1);
        vector[1] phi_ma = rep_vector(0,1);



        for (p in 1:P) {
          if (t < p + 1) {
            phi_lag = phi[T + t - p];
            alr_Y_lag = alr_Y[T + t - p];
            Xbeta_lag = Xbeta[T + t - p];
            Xbeta_phi_lag[1] = Xbeta_phi[T + t - p];
          } else {
            phi_lag = phi_new[t - p];
            alr_Y_lag = alr_Y_hat[t - p];
            Xbeta_lag = Xbeta_new[t - p];
            Xbeta_phi_lag[1] = Xbeta_phi_new[t - p];
          }

          ar += A[p] * (alr_Y_lag - Xbeta_lag);
          phi_ar += ar_phi[p]*(phi_lag-Xbeta_phi_lag);
        }

        for (q in 1:Q) {
          if (t < q + 1) {
            phi_lag = phi[T + t - q];
            alr_Y_lag = alr_Y[T + t - q];
            eta_lag = eta[T + t - q];
          } else {
             phi_lag = phi_new[t - q];
            alr_Y_lag = alr_Y_hat[t - q];
            eta_lag = eta_new[t - q];
          }

          ma += B[q] * (alr_Y_lag - eta_lag);
          phi_ma += ma_phi[q]* sum(square(alr_Y_lag - eta_lag));
        }
        phi_new[t] = phi_ar[1]+phi_ma[1]+Xbeta_phi_new[t];
        eta_new[t] = ar + ma + Xbeta_new[t];
        alpha_new[t] = exp(phi_new[t]) * alrinv(eta_new[t], ref);

        Y_hat[t] = dirichlet_rng(alpha_new[t]);
                alr_Y_hat[t] = alr(Y_hat[t], ref);
      }
    }
  }
}
