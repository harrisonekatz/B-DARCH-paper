functions {
  #include alr.stan
  #include component_product.stan
}

data {
  int<lower=1> T;                       // number of time periods
  int<lower=2> C;                       // number of categories
  array[T] simplex[C] Y;                // response array
  int<lower=1, upper=C> ref;            // ALR reference element of simplex
  int<lower=0> P;                       // number of auto-regressive lags
  int<lower=0> Q;                       // number of moving average lags
  int<lower=0> T_new;                   // number of time periods to forecast
}

transformed data {
  array[T] vector[C - 1] alr_Y;

  for (t in 1:T) {
    alr_Y[t] = alr(Y[t], ref);


  // Generate uniform random initial values for phi

     // Generate uniform random values in [3, 5]
  }

}

parameters {
  vector[C - 1] beta;                   // coefficients for covariates
  matrix[C - 1, C - 1] A;               // VAR coefficients
  real phi_int;                         // coefficients for phi covariates
  real alpha1;                          // AR coefficients for phi
  real gamma;                           // AR coefficients for phi
}

transformed parameters {
  vector[T] phi;
  array[T] vector[C - 1] eta;
  array[T] vector[C] alpha;

  // Initialize eta and phi
  eta[1] = alr_Y[1];
  phi[1] = phi_int;
  alpha[1] = exp(phi[1]) * alrinv(eta[1], ref);

  // Model the known DGP
  for (t in 2:T) {
    eta[t] = beta + A * (alr_Y[t-1] - beta);
    phi[t] =  phi_int + alpha1 * (phi[t-1] - phi_int) + gamma * sum(square((alr_Y[t-1]) - eta[t-1]));
    alpha[t] = exp(phi[t]) * alrinv(eta[t], ref);
  }
}

model {
  beta ~ normal(0, .5);
  to_vector(A) ~ normal(0, 1);
  alpha1 ~ normal(.5, .15);
  gamma ~ normal(-.75, .15);
  phi_int~normal(7,1.5);
  // Likelihood
  for (t in 2:T) {
    Y[t] ~ dirichlet(alpha[t]);
  }
}


generated quantities {
  array[T] simplex[C] Y_hat_train;      // Fitted values for training data
  array[T] simplex[C] E_Y_hat_train;    // Expected values for training data
  array[T] vector[C] alpha_hat;    // Expected values for training data

  // Generate fitted values for training data
  for (t in 1:T) {
    // Simulate Y_hat_train[t] from the Dirichlet distribution
    Y_hat_train[t] = dirichlet_rng(alpha[t]);
    // Compute the expected value of Y[t]
    E_Y_hat_train[t] = alpha[t] / sum(alpha[t]);
    alpha_hat[t] = alpha[t];
  }



  array[T_new] simplex[C] Y_hat;  // Predictions for new data
        vector[T_new] phi_new;
        array[T_new] vector[C - 1] eta_new;

  // Posterior predictive checks and predictions

  // Predictions for out-of-sample data
  if (T_new > 0) {
    {
      array[T_new] vector[C - 1] alr_Y_hat;
      array[T_new] vector[C] alpha_new;

      // Initialize first values for predictions
     for (t in 1:T_new) {
  vector[C - 1] alr_Y_lag;
  real phi_lag;
  real trend_lag;
  vector[C - 1] eta_lag;
  vector[C - 1] ar = rep_vector(0, C - 1);
  vector[C - 1] ma = rep_vector(0, C - 1);
  real phi_ar = 0;    // Initialize OUTSIDE the loop
  real phi_ma = 0;

  // Autoregressive component
  for (p in 1:P) {
    if (t >= p + 1) {
      phi_lag = phi_new[t - p];
      alr_Y_lag = alr_Y_hat[t - p];
      trend_lag=t-p;

    } else {
      // Handle out-of-bounds: use values from the training data
      phi_lag = phi[T + t - p];
      alr_Y_lag = alr_Y[T + t - p];
      trend_lag=T + t - p;
    }

    ar += A * (alr_Y_lag - beta);
    phi_ar =  alpha1 * (phi_lag - phi_int);
  }

  // Moving average component
  for (q in 1:Q) {
    if (t >= q + 1) {
      phi_lag = phi_new[t - q];
      alr_Y_lag = alr_Y_hat[t - q];
      eta_lag = eta_new[t - q];
    } else {
      // Handle out-of-bounds
      phi_lag = phi[T + t - q];
      alr_Y_lag = alr_Y[T + t - q];
      eta_lag = eta[T + t - q];
    }

    phi_ma = gamma * sum(square(alr_Y_lag - eta_lag));
  }

  // Compute new phi, eta, alpha, and predictions
  phi_new[t] = phi_ar + phi_ma + phi_int;
  eta_new[t] = ar  + beta;
  alpha_new[t] = exp(phi_new[t]) * alrinv(eta_new[t], ref);

  Y_hat[t] = dirichlet_rng(alpha_new[t]);
  alr_Y_hat[t] = alr(Y_hat[t],ref);
}
}
}}
