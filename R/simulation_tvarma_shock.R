## simulation script
library(dplyr)
library(MASS)
library(compositions)
library(brms)
library(tidyr)
library(ggplot2)
library(patchwork)
source(file="R/compositional.R")
library(purrr)

library(cmdstanr)

library(MASS)  # For multivariate normal simulation
# For rdirichlet


observations_test <- 40

observations_training <- 60




n_sim <- 50


pacf_results_btvarma_list <- vector("list", n_sim)
pacf_results_bdarch_list <- vector("list", n_sim)
pacf_results_bdarma_list <- vector("list", n_sim)

results_btvarma <- vector("list", n_sim)
results_bdarma <- vector("list", n_sim)
results_bdarch <- vector("list", n_sim)

pacf_results_tvpvar_list <- vector("list", n_sim)
results_tvpvar <- vector("list", n_sim)

A_sim <- vector("list", n_sim)
beta_sim <- vector("list", n_sim)
y_sim <- vector("list", n_sim)
phi_beta_sim <-  vector("list", n_sim)
set.seed(1007)

total_obs <- observations_test+observations_training


for(i in 1:n_sim){


  valid_simulation <- FALSE

  while(!valid_simulation) {
    # Try to run the simulation in this loop
    tryCatch({
      # Function to generate a random positive-definite covariance matrix
      generate_cov_matrix <- function(dimension) {
        random_matrix <- matrix(runif(dimension^2, -0.3, 0.3), nrow = dimension)
        cov_matrix <- t(random_matrix) %*% random_matrix  # Ensure positive-definite
        return(cov_matrix)
      }

      # Initial setup
      A  <- matrix(runif(16, -.75, .75), nrow = 4)
      beta_y <- rnorm(5, mean = 0.2, sd = 0.01)
      beta_y <- beta_y / sum(beta_y)
      beta <- alr(beta_y)
      phi_trend <- 0
      alpha <- 0.8
      gamma <- -.95
      phi <- matrix(data = NA, nrow = total_obs)
      y <- matrix(data = NA, nrow = total_obs, ncol = 5)
      y[1,] <- c(0.2, 0.2, 0.2, 0.2, 0.2) + rnorm(5, 0, sd = 0.01)
      y[1,] <- y[1,] / sum(y[1,])
      eta <- matrix(data = NA, nrow = total_obs, ncol = 4)
      eta[1,] <- alr(y[1,])
      missreported <- rpois(10, lambda = 6)
      missreported <- cumsum(missreported)
      missreported <- missreported[which(missreported <= observations_training)]

      # Generate a single sigma for the whole simulation, with a complex covariance structure
      sigma <- runif(1, 0.05, 0.5)
      cov_matrix <- generate_cov_matrix(4)  # Generate a complex covariance matrix

      # Start simulation loop for tVARMA
      for (k in 2:total_obs) {
        # Use the single sigma with the more complex covariance matrix
        eta[k,] <- mvrnorm(1, mu = beta + (A %*% (alr(y[k-1,]) - beta)), Sigma = sigma^2 * cov_matrix)

        # Calculate the transformation back to the simplex
        # Update y
        y[k,] = alrinv(eta[k,])

        # Check if any value in y[k,] is less than or equal to 0.001
        if (min(y[k,]) <= 0.001) {
          valid_simulation <- FALSE
          stop("Simulation invalid: min(y[k,]) <= 0.001")  # Force the tryCatch to restart
        }
      }

      # Handling misreported values
      for (l in 1:length(missreported)) {
        y_miss <- runif(5, min = 0, max = 1)
        y_miss <- y_miss / sum(y_miss)
        y[missreported[l],] <- y_miss
      }

      # If the simulation reaches this point without errors, mark it as valid
      valid_simulation <- TRUE

    }, error = function(e) {
      # Handle the error without stopping the script
      message("An error occurred: ", e$message)
      # valid_simulation remains FALSE, so the loop continues
    })
  }





  # Plot results
  #plot(1:total_obs, exp(phi), type = "l")
  matplot(y, type = "l", lty = 1, pch = 19, col = 1:5, xlab = "Index", ylab = "Value", main = "Plot of Matrix Columns")

  y_train <- y[1:observations_training,]
  y_test <- y[(observations_training+1):total_obs,]

  # If the simulation reaches this point without errors, mark it as valid








  y_train_df <- y_train |>
    as.data.frame() |>
    dplyr::mutate(idx = row_number()) |>
    tidyr::pivot_longer(-idx) |>
    dplyr::mutate(across(name, factor)) %>%
    group_by(idx) %>%
    summarise(proportions = list(value)) %>%
    pull(proportions)






  # Identify which elements are not the same size


  ds_var <- "t"
  compositional_var <- "currency"


  T_new <- observations_test
  forecast_window_days <- observations_test

  max_actuals_time <- length(y_train_df)




  C <- 5

  # Check that each vector is a simplex (sums to 1)



  results_list_bdama <- list()
  results_list_bdarch <- list()
  results_list_tvaram <- list()

  results_list_tvpvar <- list()





  # Print the final combined data frame



  P_x <- 1
  Q_x <- 0

  P_phi <- 1
  Q_phi <- 1




  X_component <- as.list(rep(1,observations_training))



  X_component <- purrr::map(
    seq_len(observations_training),
    ~ rep(1, 4)  # Each element is a vector of length 3 filled with 1's
  )


  X_new_component <- purrr::map(
    seq_len(observations_test),
    ~ rep(1, 4)  # Each element is a vector of length 3 filled with 1's
  )


  X_phi <- as.list(rep(1,observations_training))





  X_phi <- do.call(rbind, X_phi)





  X_phi_new <- as.list(rep(1,observations_test))


  X_phi_new <- do.call(rbind, X_phi_new)






  stan_data <- list(
    `T` = length(y_train_df),
    C = C,
    Y = y_train_df,
    ref = 1,
    N = length(X_component[[1]]),
    K = purrr::map_int(X_component[[1]], length),
    X = X_component,
    K_phi = ncol(X_phi),
    X_phi = X_phi,
    P = 1,
    Q = 0,
    T_new = T_new,
    X_new = X_new_component,
    X_phi_new = X_phi_new,
    prior_only = 0,
    beta_sd = 0.3
  )

  #cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.31.0/")

  stan_dir <- "stan"

  stan_file <- file.path(stan_dir, "bdarma.stan")

  mod <- cmdstanr::cmdstan_model(stan_file, include_paths = stan_dir)



  fit_bdarma <- mod$sample(
    data = stan_data,
    seed = 1234,
    init = .25,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 250,
    iter_sampling = 250
  )





  Y_hat_bdarma <- fit_bdarma$draws("Y_hat", format = "draws_df")



  Y_hat_long_raw_bdarma <-
    Y_hat_bdarma |>
    dplyr::select(.draw, starts_with("Y_hat")) |>
    tidyr::pivot_longer(starts_with("Y_hat"))


  Y_hat_long_bdarma <-
    Y_hat_long_raw_bdarma |>
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      qlo = quantile(value, probs = 0.05, na.rm = TRUE, names = FALSE),
      qhi = quantile(value, probs = 0.95, na.rm = TRUE, names = FALSE),
      .by = name
    ) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",") |>
    dplyr::select(all_of(c(ds_var, compositional_var)), mean, qlo, qhi)



  Y_hat_long_bdarma$forecast <- "bdarma"



  Y_hat_long_bdarma$t <- as.numeric(Y_hat_long_bdarma$t)

  colnames(y_test) <- c("1","2","3","4","5")

  colnames(y_train) <- c("1","2","3","4","5")




  y_train_long <- as_tibble(y_train) %>%
    mutate(t = 1:nrow(y_train)) %>%
    pivot_longer(cols = -t, names_to = "currency", values_to = "mean") %>%
    mutate(currency = as.character(as.numeric(currency)), # Ensure currency is character type
           qlo = NA,  # Placeholder for qlo
           qhi = NA,  # Placeholder for qhi
           forecast = "actuals") %>%
    dplyr::select(t, currency, mean, forecast) # Ensure correct column order




  y_test_long <- as_tibble(y_test) %>%
    mutate(t = 1:nrow(y_test)) %>%
    pivot_longer(cols = -t, names_to = "currency", values_to = "mean") %>%
    mutate(currency = as.character(as.numeric(currency)), # Ensure currency is character type
           qlo = NA,  # Placeholder for qlo
           qhi = NA,  # Placeholder for qhi
           forecast = "actuals") %>%
    dplyr::select(t, currency, mean, qlo, qhi, forecast) # Ensure correct column order

  y_test_long$t <- (y_test_long$t)

  validation_r1_bdarma <- full_join(Y_hat_long_bdarma,y_test_long)


  bdarma_valid <-  ggplot(validation_r1_bdarma, aes(x = t, y = mean, color = forecast)) +
    geom_line() +
    geom_ribbon(aes(ymin = qlo, ymax = qhi), alpha = 0.1) +
    facet_wrap(~currency, scales = 'free_y') +
    labs(title = "B-DARMA: forecast vs actuals",
         x = "date",
         y = "proportion",
         color = "forecast") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))








  ##### y train




  alpha_train <- fit_bdarma$draws("alpha_hat", format = "draws_df")



  alpha_hat_long <-
    alpha_train |>
    dplyr::select(.draw, starts_with("alpha_hat")) |>
    tidyr::pivot_longer(starts_with("alpha_hat"))



  alpha_hat_long_bdarma <-
    alpha_hat_long |>
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      .by = name
    ) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",") |>
    dplyr::select(all_of(c(ds_var, compositional_var)), mean)

  alpha_hat_long_bdarma <- alpha_hat_long_bdarma %>%
    group_by(t) %>%                      # Group by time point
    mutate(alpha_0 = sum(mean),           # Compute the sum of means (alpha_0) for each time point
           sd = sqrt((mean * (alpha_0 - mean)) / (alpha_0^2 * (alpha_0 + 1)))) %>%
    ungroup() %>% dplyr::select(t,currency,sd)


  Y_hat_bdarma_train <- fit_bdarma$draws("Y_hat_train", format = "draws_df")



  Y_hat_long_raw_bdarma_train <-
    Y_hat_bdarma_train |>
    dplyr::select(.draw, starts_with("Y_hat")) |>
    tidyr::pivot_longer(starts_with("Y_hat"))


  Y_hat_long_bdarma_train <-
    Y_hat_long_raw_bdarma_train |>
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      .by = name
    ) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",") |>
    dplyr::select(all_of(c(ds_var, compositional_var)), mean)



  Y_hat_long_bdarma_train$forecast <- "bdarma"



  Y_hat_long_bdarma_train$t <- as.numeric(Y_hat_long_bdarma_train$t)


  validation_r1_bdarma_train <- full_join(Y_hat_long_bdarma_train,y_train_long)


  bdarma_valid <-  ggplot(validation_r1_bdarma, aes(x = t, y = mean, color = forecast)) +
    geom_line() +
    geom_ribbon(aes(ymin = qlo, ymax = qhi), alpha = 0.1) +
    facet_wrap(~currency, scales = 'free_y') +
    labs(title = "B-DARMA: forecast vs actuals",
         x = "date",
         y = "proportion",
         color = "forecast") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))






  posterior_samples <- as_draws_df(fit_bdarma)
  # Extract and reshape `beta`, `beta_phi`, `ar_phi`, `ma_phi`
  coefficients <- posterior_samples %>%
    dplyr::select(starts_with("beta"), starts_with("ar_phi"), starts_with("ma_phi")) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")



  unique(coefficients$parameter)



  # Extract and reshape matrices A and B
  # Extract A
  A_samples <- posterior_samples %>%
    dplyr::select(matches("^A\\[.*\\]")) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")


  unique(A_samples$parameter)
  # Extract B
  B_samples <- posterior_samples %>%
    dplyr::select(starts_with("B")) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

  # Combine all the data into a single data frame for plotting
  all_parameters <- bind_rows(
    coefficients %>% mutate(type = "Coefficients"),
    A_samples %>% mutate(type = "A (VAR)"),
    B_samples %>% mutate(type = "B (VMA)")
  )

  # Plot the densities for all parameters
  ggplot(all_parameters, aes(x = value, fill = parameter)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ parameter, scales = "free") +
    theme_minimal() +
    labs(title = "Density of Model Parameters", x = "Value", y = "Density")






  forecasts_bdarma <- validation_r1_bdarma %>%
    filter(forecast == "bdarma") %>%
    dplyr::select(t, currency, mean) %>%
    rename(forecast_mean = mean)

  actuals_bdarma <- validation_r1_bdarma %>%
    filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    rename(actual_mean = mean)


  forecasts_bdarma_train <- validation_r1_bdarma_train %>%
    filter(forecast == "bdarma") %>%
    dplyr::select(t, currency, mean) %>%
    rename(forecast_mean = mean)

  actuals_bdarma_train <- validation_r1_bdarma_train %>%
    filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    rename(actual_mean = mean)






  alr_transformed_forecast <- forecasts_bdarma %>%
    group_by(t) %>%
    mutate(alr_value_forecast = if_else(currency != "1",
                                        log(forecast_mean / forecast_mean[currency == "1"]),
                                        NA_real_)) %>%
    ungroup() %>% na.omit()


  alr_transformed_actuals <- actuals_bdarma %>%
    group_by(t) %>%
    mutate(alr_value_actual = if_else(currency != "1",
                                      log(actual_mean / actual_mean[currency == "1"]),
                                      NA_real_)) %>%
    ungroup() %>% na.omit()



  merged_df_bdarma <- forecasts_bdarma %>%
    inner_join(actuals_bdarma, by = c("t", "currency"))


  merged_df_bdarma_train <- forecasts_bdarma_train %>%
    inner_join(actuals_bdarma_train, by = c("t", "currency"))


  merged_df_alr_bdarma <- alr_transformed_forecast %>%
    inner_join(alr_transformed_actuals, by = c("t", "currency"))


  colnames(merged_df_alr_bdarma)


  squared_residual_bdarma <-  merged_df_alr_bdarma %>%
    mutate(squared_error = (alr_value_forecast - alr_value_actual)^2) %>%
    group_by(t) %>%
    summarise(value=sum(squared_error))




  metrics_bdarma <-  merged_df_bdarma %>%
    mutate(squared_error = (actual_mean - forecast_mean)^2) %>%
    mutate(absolute_error=abs((actual_mean - forecast_mean))) %>%
    ungroup() %>%
    summarise(SSR=sum(squared_error),MAE=mean(absolute_error),mSSR=sqrt(mean(squared_error)),TAE=sum(absolute_error))


  alpha_hat_long_bdarma$t <- as.numeric(alpha_hat_long_bdarma$t)


  head(alpha_hat_long_bdarma)
  head(merged_df_bdarma_train)
  squared_residual_bdarma_series <- merged_df_bdarma_train %>%
    left_join(alpha_hat_long_bdarma)%>%
    arrange(t,currency) %>%
    mutate(standard_error=(actual_mean - forecast_mean)/sd) %>%
    group_by(t) %>% summarise(squared_error=sum(standard_error^2))







  pacf_results_bdarma <- squared_residual_bdarma_series %>%
    do({
      # Extract squared errors for the current group
      squared_errors <- .$squared_error

      # Calculate PACF
      pacf_result <- pacf(squared_errors, plot = FALSE)

      # Create a tibble with PACF values and lags
      tibble(
        lag = 1:length(pacf_result$acf),
        pacf_value = pacf_result$acf
      )
    })


  g1 <- ggplot(pacf_results_bdarma, aes(x = lag, y = pacf_value)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, fill = "skyblue") +
    theme_minimal()+
    labs(title = "PACF of B-DARMA ",subtitle = "DGP: DARCH(1,1)", x = "Lag", y = "PACF Value") +
    geom_hline(yintercept = 0.2, linetype = "dotted", color = "blue") +
    geom_hline(yintercept = -0.2, linetype = "dotted", color = "blue") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylim(-.75,.75)




  stan_data <- list(
    `T` = length(y_train_df),
    C = C,
    Y = y_train_df,
    ref = 1,
    N = length(X_component[[1]]),
    K = purrr::map_int(X_component[[1]], length),
    X = X_component,
    K_phi = ncol(X_phi),
    X_phi = X_phi,
    P = 1,
    Q = 0,
    T_new = T_new,
    X_new = X_new_component,
    X_phi_new = X_phi_new,
    prior_only = 0,
    beta_sd = 0.3
  )


  #cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.31.0/")

  stan_file_darch <- file.path(stan_dir, "bdarch_pq_sim_trend.stan")

  mod_darch <- cmdstanr::cmdstan_model(stan_file_darch, include_paths = stan_dir)


  init_function <- function() {
    list(               # Initial value for beta_phi_ma
      phi = rep(800, `T`)                  # Initial values for phi
    )
  }


  fit_darch <- mod_darch$sample(
    data = stan_data,
    seed = 1234,
    init = .25,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 250,
    iter_sampling = 250,
    max_treedepth = 10,
    adapt_delta = .8
  )








  posterior_samples <- as_draws_df(fit_darch)
  # Extract and reshape `beta`, `beta_phi`, `ar_phi`, `ma_phi`
  coefficients <- posterior_samples %>%
    dplyr::select(starts_with("beta"), starts_with("alpha1"), starts_with("gamma"),starts_with("phi_int")) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")




  summary_stats <- coefficients %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      `2.5%` = quantile(value, 0.025),
      `50%` = quantile(value, 0.5),  # Median
      `97.5%` = quantile(value, 0.975),
      .groups = 'drop'
    )

  # Display the summary statistics
  print(summary_stats)


  unique(coefficients$parameter)



  # Extract and reshape matrices A and B
  # Extract A
  A_samples <- posterior_samples %>%
    dplyr::select(matches("^A\\[.*\\]")) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")



  summary_stats_A <- A_samples %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      `2.5%` = quantile(value, 0.025),
      `50%` = quantile(value, 0.5),  # Median
      `97.5%` = quantile(value, 0.975),
      .groups = 'drop'
    )

  print(summary_stats_A)




  # Extract B
  B_samples <- posterior_samples %>%
    dplyr::select(starts_with("B")) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

  summary_stats_B <- B_samples %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      `2.5%` = quantile(value, 0.025),
      `50%` = quantile(value, 0.5),  # Median
      `97.5%` = quantile(value, 0.975),
      .groups = 'drop'
    )

  print(summary_stats_B)


  # Combine all the data into a single data frame for plotting
  all_parameters <- bind_rows(
    coefficients %>% mutate(type = "Coefficients"),
    A_samples %>% mutate(type = "A (VAR)"),
    B_samples %>% mutate(type = "B (VMA)")
  )

  # Plot the densities for all parameters
  ggplot(all_parameters, aes(x = value, fill = parameter)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ parameter, scales = "free") +
    theme_minimal() +
    labs(title = "Density of Model Parameters", x = "Value", y = "Density")

  Y_hat <- fit_darch$draws("Y_hat", format = "draws_df")


  Y_hat_train <- fit_darch$draws("Y_hat_train", format = "draws_df")



  alpha_train <- fit_darch$draws("alpha_hat", format = "draws_df")



  Y_hat_long_raw_bdarch <-
    Y_hat |>
    dplyr::select(.draw, starts_with("Y_hat")) |>
    tidyr::pivot_longer(starts_with("Y_hat"))


  Y_hat_long_train_raw_bdarch <-
    Y_hat_train |>
    dplyr::select(.draw, starts_with("Y_hat_train")) |>
    tidyr::pivot_longer(starts_with("Y_hat_train"))


  alpha_hat_long <-
    alpha_train |>
    dplyr::select(.draw, starts_with("alpha_hat")) |>
    tidyr::pivot_longer(starts_with("alpha_hat"))



  alpha_hat_long_bdarch <-
    alpha_hat_long |>
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      .by = name
    ) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",") |>
    dplyr::select(all_of(c(ds_var, compositional_var)), mean)

  alpha_hat_long_bdarch <- alpha_hat_long_bdarch %>%
    group_by(t) %>%                      # Group by time point
    mutate(alpha_0 = sum(mean),           # Compute the sum of means (alpha_0) for each time point
           sd = sqrt((mean * (alpha_0 - mean)) / (alpha_0^2 * (alpha_0 + 1)))) %>%
    ungroup() %>% dplyr::select(t,currency,sd)


  Y_hat_long_train_bdarch <-
    Y_hat_long_train_raw_bdarch |>
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      .by = name
    ) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",") |>
    dplyr::select(all_of(c(ds_var, compositional_var)), mean)

  Y_hat_long_bdarch <-
    Y_hat_long_raw_bdarch |>
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      qlo = quantile(value, probs = 0.05, na.rm = TRUE, names = FALSE),
      qhi = quantile(value, probs = 0.95, na.rm = TRUE, names = FALSE),
      .by = name
    ) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",") |>
    dplyr::select(all_of(c(ds_var, compositional_var)), mean, qlo, qhi)




  Y_hat_long_bdarch$forecast <- "bdarch"


  Y_hat_long_train_bdarch$forecast <- "bdarch"


  Y_hat_long_train_bdarch$t <- as.numeric(Y_hat_long_train_bdarch$t)



  Y_hat_long_bdarch$t <- as.numeric(Y_hat_long_bdarch$t)



  validation_r1_bdarch <- full_join(Y_hat_long_bdarch,y_test_long)

  validation_r1_bdarch_train <- full_join(Y_hat_long_train_bdarch,y_train_long)


  validation_r1_bdarch_bdarma <- full_join(Y_hat_long_bdarma,validation_r1_bdarch)

  bdarch_valid <-  ggplot(validation_r1_bdarch_bdarma, aes(x = t, y = mean, color = forecast)) +
    geom_line() +
    geom_ribbon(aes(ymin = qlo, ymax = qhi), alpha = 0.1) +
    facet_wrap(~currency, scales = 'free_y') +
    labs(title = "DGP: DARMA(1,1)-DARCH(1,1)",
         subtitle="forecast vs actuals",
         x = "t",
         y = "proportion",
         color = "forecast") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  bdarch_train <-  ggplot(validation_r1_bdarch_train, aes(x = t, y = mean, color = forecast)) +
    geom_line() +
    facet_wrap(~currency, scales = 'free_y') +
    labs(title = "DGP: DARMA(1,1)-DARCH(1,1)",
         subtitle="forecast vs actuals",
         x = "t",
         y = "proportion",
         color = "forecast") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  bdarch_train

  bdarch_valid





  forecasts_bdarch_train <- validation_r1_bdarch_train %>%
    filter(forecast == "bdarch") %>%
    dplyr::select(t, currency, mean) %>%
    rename(forecast_mean = mean)

  actuals_bdarch_train <- validation_r1_bdarch_train %>%
    filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    rename(actual_mean = mean)


  # Merge the forecasts and actuals
  merged_df_bdarch_train <- forecasts_bdarch_train %>%
    inner_join(actuals_bdarch_train, by = c("t", "currency"))



  alpha_hat_long_bdarch$t <- as.numeric(alpha_hat_long_bdarch$t)


  squared_residual_bdarch_series_train <- merged_df_bdarch_train %>%
    left_join(alpha_hat_long_bdarch)%>%
    arrange(t,currency) %>%
    mutate(standard_error=(actual_mean - forecast_mean)/sd) %>%
    group_by(t) %>% summarise(squared_error=sum(standard_error^2))





  forecasts_bdarch <- validation_r1_bdarch %>%
    filter(forecast == "bdarch") %>%
    dplyr::select(t, currency, mean) %>%
    rename(forecast_mean = mean)

  actuals_bdarch <- validation_r1_bdarch %>%
    filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    rename(actual_mean = mean)




  #
  alr_transformed_forecast <- forecasts_bdarch %>%
    group_by(t) %>%
    mutate(alr_value_forecast = if_else(currency != "1",
                                        log(forecast_mean / forecast_mean[currency == "1"]),
                                        NA_real_)) %>%
    ungroup() %>% na.omit()


  alr_transformed_actuals <- actuals_bdarch %>%
    group_by(t) %>%
    mutate(alr_value_actual = if_else(currency != "1",
                                      log(actual_mean / actual_mean[currency == "1"]),
                                      NA_real_)) %>%
    ungroup() %>% na.omit()




  # Merge the forecasts and actuals
  merged_df_bdarch <- forecasts_bdarch %>%
    inner_join(actuals_bdarch, by = c("t", "currency"))



  merged_df_alr_btvarma <- alr_transformed_forecast %>%
    inner_join(alr_transformed_actuals, by = c("t", "currency"))


  head(merged_df_alr_btvarma)


  squared_residual_btvarma <-  merged_df_alr_btvarma %>%
    mutate(squared_error = (alr_value_actual - alr_value_forecast)^2) %>%
    group_by(t) %>%
    summarise(value=sum(squared_error))



  squared_residual_bdarch_series <- merged_df_alr_btvarma %>%
    arrange(t,currency) %>%
    mutate(squared_error = (actual_mean - forecast_mean)^2)%>%
    group_by(t) %>% summarise(squared_error=sum(squared_error))



  head(merged_df_alr_btvarma)

  squared_residual_bdarch_series <- merged_df_alr_btvarma %>%
    left_join(alpha_hat_long_bdarch)%>%
    arrange(t,currency) %>%
    mutate(standard_error=(actual_mean - forecast_mean)/sd) %>%
    group_by(t) %>% summarise(squared_error=sum(standard_error^2))





  metrics_bdarch <-  merged_df_bdarch %>%
    mutate(squared_error = (actual_mean - forecast_mean)^2) %>%
    mutate(absolute_error=abs((actual_mean - forecast_mean))) %>%
    ungroup() %>%
    summarise(SSR=sum(squared_error),MAE=mean(absolute_error),mSSR=sqrt(mean(squared_error)),TAE=sum(absolute_error))

  pacf_results_bdarch_series <- squared_residual_bdarch_series_train %>%
    do({
      # Extract squared errors for the current group
      squared_errors <- .$squared_error

      # Calculate PACF
      pacf_result <- pacf(squared_errors, plot = FALSE)

      # Create a tibble with PACF values and lags
      tibble(
        lag = 1:length(pacf_result$acf),
        pacf_value = pacf_result$acf
      )
    })


  g2 <- ggplot(pacf_results_bdarch_series, aes(x = lag, y = pacf_value)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, fill = "skyblue") +
    labs(title = "PACF of B-DARCH ",subtitle = "DGP: DARCH(1,1)", x = "Lag", y = "PACF Value") +
    theme_minimal() +
    geom_hline(yintercept = 0.2, linetype = "dotted", color = "blue") +
    geom_hline(yintercept = -0.2, linetype = "dotted", color = "blue") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylim(-.75,.75)




  pacf_results_btvarma <- squared_residual_btvarma %>%
    do({
      # Extract squared errors for the current group
      squared_errors <- .$value

      # Calculate PACF
      pacf_result <- pacf(squared_errors, plot = FALSE)

      # Create a tibble with PACF values and lags
      tibble(
        lag = 1:length(pacf_result$acf),
        pacf_value = pacf_result$acf
      )
    })



  pacf_bdarch <- ggplot(pacf_results_btvarma, aes(x = lag, y = pacf_value)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, fill = "skyblue") +
    labs(title = "PACF of B-DARCH ", x = "Lag", y = "PACF Value") +
    theme_minimal() + ylim(-.5,.5)+
    geom_hline(yintercept = 0.2, linetype = "dotted", color = "blue") +
    geom_hline(yintercept = -0.2, linetype = "dotted", color = "blue") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))



  ### tvarma









  y_train_df_alr <- as.data.frame(y_train)
  y_train_df_alr$t <- 1:observations_training
  y_train_df_alr <- matrix(data=NA,nrow=observations_training,ncol=4)

  for(j in 1:observations_training){
    y_train_df_alr[j,]<- alr(y_train[j,1:5])
  }

  y_train_df_alr <- as.data.frame(y_train_df_alr)
  y_train_df_alr$t <- 1:observations_training



  y_train_df_long <- y_train_df_alr %>%
    pivot_longer(cols = starts_with("V"), names_to = "variable", values_to = "alr_value_actual")

  # Grouping by 't' and summarizing as a list
  proportions_list <- y_train_df_long %>%
    group_by(t) %>%
    summarise(proportions = list(alr_value_actual)) %>%
    pull(proportions)





  stan_data <- list(
    `T` = length(proportions_list),
    C = C-1,
    Y = proportions_list,
    N = length(X_component[[1]]),
    K = purrr::map_int(X_component[[1]], length),
    X = X_component,
    K_phi = ncol(X_phi),
    X_phi = X_phi,
    P = 1,
    Q = 0,
    T_new = T_new,
    X_new = X_new_component,
    prior_only = 0,
    beta_sd = 1
  )

  cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.31.0/")

  stan_dir <- "stan"

  stan_file <- file.path(stan_dir, "varma.stan")

  mod_varma <- cmdstanr::cmdstan_model(stan_file, include_paths = stan_dir)



  fit_varma <- mod_varma$sample(
    data = stan_data,
    seed = 1234,
    init = .25,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 250,
    iter_sampling = 250
  )




  Y_hat_varma_sigma <- fit_varma$draws("sigma", format = "draws_df")



  Y_hat_varma_sigma<-
    Y_hat_varma_sigma |>
    dplyr::select(.draw, starts_with("sigma")) |>
    tidyr::pivot_longer(starts_with("sigma")) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])"))%>%
    group_by(name) %>% summarise(value=mean(value))


  Y_hat_varma_sigma$currency <- 1:4


  Y_hat_varma_sigma <- Y_hat_varma_sigma %>% dplyr::select(sd=value,currency)

  Y_hat_varma_train <- fit_varma$draws("Y_hat_train", format = "draws_df")

  head(Y_hat_varma_sigma)

  Y_hat_long_raw_varma_train <-
    Y_hat_varma_train |>
    dplyr::select(.draw, starts_with("Y_hat")) |>
    tidyr::pivot_longer(starts_with("Y_hat")) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",")

  Y_hat_long_raw_varma_exp <- Y_hat_long_raw_varma_train
  Y_hat_long_raw_varma_exp$value <- exp(Y_hat_long_raw_varma_exp$value)


  Y_hat_long_raw_varma_exp <- Y_hat_long_raw_varma_exp[,c(1,3,4,5)]

  Y_hat_wide_train <- Y_hat_long_raw_varma_exp %>%
    pivot_wider(
      names_from = currency,
      values_from = value,
      names_prefix = "currency_",
      values_fill = list(value = 0)  # Fill missing values with 0 if necessary
    ) %>%
    rename_with(~ paste0(.x, "_value"), starts_with("currency_"))

  Y_hat_transformed_train <- Y_hat_wide_train %>%
    group_by(.draw, t) %>%
    mutate(across(starts_with("currency_"), ~ .x / (1 + sum(c_across(starts_with("currency_")))))) %>%
    ungroup()



  Y_hat_transformed_train$currency_0_value <- 1-(Y_hat_transformed_train$currency_1_value+Y_hat_transformed_train$currency_2_value+Y_hat_transformed_train$currency_3_value+Y_hat_transformed_train$currency_4_value)


  Y_hat_long_train <- Y_hat_transformed_train %>%
    pivot_longer(
      cols = starts_with("currency_"),
      names_to = "currency",
      values_to = "value"
    ) %>%
    mutate(currency = sub("currency_", "", currency),   # Remove "currency_" prefix
           currency = sub("_value", "", currency))



  Y_hat_long_raw_varm_train <-
    Y_hat_long_train |>
    group_by(t,currency) %>%
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      qlo = quantile(value, probs = 0.05, na.rm = TRUE, names = FALSE),
      qhi = quantile(value, probs = 0.95, na.rm = TRUE, names = FALSE),
    ) |>
    dplyr::select(all_of(c(ds_var, compositional_var)), mean, qlo, qhi)


  Y_hat_long_raw_varm_train$currency <- as.numeric(Y_hat_long_raw_varm_train$currency)+1

  Y_hat_long_raw_varm_train$currency <- paste0(Y_hat_long_raw_varm_train$currency)






  Y_hat_varma <- fit_varma$draws("Y_hat", format = "draws_df")

  Y_hat_long_raw_varma <-
    Y_hat_varma |>
    dplyr::select(.draw, starts_with("Y_hat")) |>
    tidyr::pivot_longer(starts_with("Y_hat")) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",")

  Y_hat_long_raw_varma_exp <- Y_hat_long_raw_varma
  Y_hat_long_raw_varma_exp$value <- exp(Y_hat_long_raw_varma_exp$value)


  Y_hat_long_raw_varma_exp <- Y_hat_long_raw_varma_exp[,c(1,3,4,5)]

  Y_hat_wide <- Y_hat_long_raw_varma_exp %>%
    pivot_wider(
      names_from = currency,
      values_from = value,
      names_prefix = "currency_",
      values_fill = list(value = 0)  # Fill missing values with 0 if necessary
    ) %>%
    rename_with(~ paste0(.x, "_value"), starts_with("currency_"))

  Y_hat_transformed <- Y_hat_wide %>%
    group_by(.draw, t) %>%
    mutate(across(starts_with("currency_"), ~ .x / (1 + sum(c_across(starts_with("currency_")))))) %>%
    ungroup()



  Y_hat_transformed$currency_0_value <- 1-(Y_hat_transformed$currency_1_value+Y_hat_transformed$currency_2_value+Y_hat_transformed$currency_3_value+Y_hat_transformed$currency_4_value)


  Y_hat_long <- Y_hat_transformed %>%
    pivot_longer(
      cols = starts_with("currency_"),
      names_to = "currency",
      values_to = "value"
    ) %>%
    mutate(currency = sub("currency_", "", currency),   # Remove "currency_" prefix
           currency = sub("_value", "", currency))



  Y_hat_long_raw_varma <-
    Y_hat_long |>
    group_by(t,currency) %>%
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      qlo = quantile(value, probs = 0.05, na.rm = TRUE, names = FALSE),
      qhi = quantile(value, probs = 0.95, na.rm = TRUE, names = FALSE),
    ) |>
    dplyr::select(all_of(c(ds_var, compositional_var)), mean)


  Y_hat_long_raw_varma$currency <- as.numeric(Y_hat_long_raw_varma$currency)+1

  Y_hat_long_raw_varma$currency <- paste0(Y_hat_long_raw_varma$currency)



  final_result <- Y_hat_long_raw_varma


  final_result_train <- Y_hat_long_raw_varm_train

  final_result_train$forecast <- "btvarma"
  final_result_train$t <- as.numeric(final_result_train$t)



  final_result$forecast <- "btvarma"
  final_result$t <- as.numeric(final_result$t)


  validation_r1_btvarma <- full_join(final_result,y_test_long)


  validation_r1_btvarma_train <- full_join(final_result_train,y_train_long)

  validation_r1_bdarch_bdarma_btvarma <- full_join(validation_r1_btvarma,validation_r1_bdarch_bdarma)



  validation_r1_bdarch_bdarma_btvarma_var <-
    validation_r1_bdarch_bdarma_btvarma %>%
    mutate(log_mean=log(mean)) %>%
    group_by(t,forecast) %>%
    summarise(total_var=mean(var(log_mean)))



  ggplot(validation_r1_bdarch_bdarma_btvarma_var, aes(x = t, y = total_var, color = forecast)) +
    geom_line() +
    labs(title = "DGP: DARMA(1,1)-DARCH(1,1)",
         subtitle="forecast vs actuals",
         x = "t",
         y = "proportion",
         color = "forecast") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  bdarch_valid <-  ggplot(validation_r1_bdarch_bdarma_btvarma, aes(x = t, y = mean, color = forecast)) +
    geom_line() +
    geom_ribbon(aes(ymin = qlo, ymax = qhi), alpha = 0.1) +
    facet_wrap(~currency, scales = 'free_y') +
    labs(title = "DGP: DARMA(1,1)-DARCH(1,1)",
         subtitle="forecast vs actuals",
         x = "t",
         y = "proportion",
         color = "forecast") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  bdarch_valid

  #### need to backtransform VARMA first


  forecasts_btvarma <- validation_r1_btvarma %>%
    filter(forecast == "btvarma") %>%
    dplyr::select(t, currency, mean) %>%
    rename(forecast_mean = mean)

  actuals_tvarma <- validation_r1_btvarma %>%
    filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    rename(actual_mean = mean)





  forecasts_btvarma_train <- validation_r1_btvarma_train %>%
    filter(forecast == "btvarma") %>%
    dplyr::select(t, currency, mean) %>%
    rename(forecast_mean = mean)

  actuals_tvarma_train <- validation_r1_btvarma_train %>%
    filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    rename(actual_mean = mean)

  head(forecasts_btvarma_train)



  forecasts_btvarma <- validation_r1_btvarma %>%
    filter(forecast == "btvarma") %>%
    dplyr::select(t, currency, mean) %>%
    rename(forecast_mean = mean)

  actuals_tvarma <- validation_r1_btvarma %>%
    filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    rename(actual_mean = mean)





  alr_transformed_forecast <- forecasts_btvarma_train %>%
    group_by(t) %>%
    mutate(alr_value_forecast = if_else(currency != "1",
                                        log(forecast_mean / forecast_mean[currency == "1"]),
                                        NA_real_)) %>%
    ungroup() %>% na.omit()


  alr_transformed_actuals <- actuals_tvarma_train %>%
    group_by(t) %>%
    mutate(alr_value_actual = if_else(currency != "1",
                                      log(actual_mean / actual_mean[currency == "1"]),
                                      NA_real_)) %>%
    ungroup() %>% na.omit()
  #


  alr_transformed_forecast$currency <- as.numeric(alr_transformed_forecast$currency)
  alr_transformed_actuals$currency <- as.numeric(alr_transformed_actuals$currency)
  Y_hat_varma_sigma$currency <- as.numeric(Y_hat_varma_sigma$currency)+1

  merged_df_btvarma <- actuals_tvarma %>%
    inner_join(forecasts_btvarma, by = c("t", "currency"))




  merged_df_btvarma_train <- actuals_tvarma_train %>%
    inner_join(forecasts_btvarma_train, by = c("t", "currency"))


  merged_df_alr_btvarma <- alr_transformed_forecast %>%
    inner_join(alr_transformed_actuals, by = c("t", "currency"))%>%
    full_join(Y_hat_varma_sigma)



  sigma_samples <- fit_varma$draws("sigma", format = "draws_matrix")
  L_Omega_samples <- fit_varma$draws("L_Omega", format = "draws_array")
  dim(L_Omega_samples)

  iterations <- dim(L_Omega_samples)[1]
  chains <- dim(L_Omega_samples)[2]
  # Number of categories

  # Combine iterations and chains
  total_draws <- iterations * chains

  # Reshape the array to combine iterations and chains
  L_Omega_samples_combined <- array(
    data = L_Omega_samples,
    dim = c(total_draws, 4, 4)
  )

  L_Omega_mean <- apply(L_Omega_samples_combined, c(2, 3), mean)
  # Compute posterior means
  sigma_mean <- colMeans(sigma_samples)
  Omega_mean <- L_Omega_mean %*% t(L_Omega_mean)
  dim(L_Omega_mean)
  dim(Omega_mean)
  length(sigma_mean)
  Sigma_mean <- diag(sigma_mean) %*% Omega_mean %*% diag(sigma_mean)

  # Step 6: Compute the Cholesky factor of 'Sigma_mean'
  L_Sigma <- t(chol(Sigma_mean))


  merged_df_alr_btvarma <- merged_df_alr_btvarma %>%
    mutate(residual = alr_value_forecast - alr_value_actual)
  merged_df_alr_btvarma <- merged_df_alr_btvarma %>%
    arrange(t, currency)

  grouped_residuals <- merged_df_alr_btvarma %>%
    group_by(t) %>%
    summarise(residuals = list(residual))

  grouped_residuals <- grouped_residuals %>%
    mutate(
      standardized_residuals = map(residuals, ~ as.numeric(solve(L_Sigma, .x)))
    )
  grouped_residuals <- grouped_residuals %>%
    mutate(
      squared_error = map_dbl(standardized_residuals, ~ sum(.x^2))
    )

  squared_residual_btvarma <- grouped_residuals %>%
    dplyr::select(t, squared_error)

  # squared_residual_btvarma <-  merged_df_alr_btvarma %>%
  #   mutate(squared_error = (alr_value_forecast - alr_value_actual)/sqrt(sd)) %>%
  #   mutate(squared_error=squared_error^2) %>%
  #   group_by(t) %>%
  #   summarise(value=sum(squared_error))


  metrics_varma <-  merged_df_btvarma %>%
    mutate(squared_error = (actual_mean - forecast_mean)^2) %>%
    mutate(absolute_error=abs((actual_mean - forecast_mean))) %>%
    ungroup() %>%
    summarise(SSR=sum(squared_error),MAE=mean(absolute_error),mSSR=sqrt(mean(squared_error)),TAE=sum(absolute_error))


  squared_residual_btvarma_series <- merged_df_btvarma_train %>%
    arrange(t,currency) %>%
    mutate(squared_error = (actual_mean - forecast_mean)^2)%>%
    group_by(t) %>% summarise(squared_error=sum(squared_error))


  pacf_results_btvama_series <- squared_residual_btvarma %>%
    do({
      # Extract squared errors for the current group
      squared_errors <- .$squared_error

      # Calculate PACF
      pacf_result <- pacf(squared_errors, plot = FALSE)

      # Create a tibble with PACF values and lags
      tibble(
        lag = 1:length(pacf_result$acf),
        pacf_value = pacf_result$acf
      )
    })

  pacf_results_btvama_series


  g3 <-ggplot(pacf_results_btvama_series, aes(x = lag, y = pacf_value)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, fill = "skyblue") +
    labs(title = "PACF of BtVARMA", subtitle = "DGP: DARCH(1,1)", x = "Lag", y = "PACF Value") +
    geom_hline(yintercept = 0.2, linetype = "dotted", color = "blue") +
    geom_hline(yintercept = -0.2, linetype = "dotted", color = "blue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(-.75, .75)

  g2 /g1/g3






  # -----------------------------
  # TVP-VAR (homoscedastic) block
  # -----------------------------

  # Stan data for TVP-VAR on the ALR scale (C-1)
  stan_data_tvpvar <- list(
    T          = length(proportions_list),
    C          = C - 1,                         # ALR dimension
    Y          = proportions_list,              # list of length T with numeric(C-1)
    N          = length(X_component[[1]]),      # matches your VARMA setup (intercepts per eq.)
    K          = purrr::map_int(X_component[[1]], length),  # typically rep(1, C-1)
    X          = X_component,
    P          = 1,
    T_new      = T_new,
    X_new      = X_new_component,
    prior_only = 0,
    beta_sd    = 1,
    rw_scale_sd = 0.1                           # RW scale prior; adjust if needed
  )

  stan_file_tvpvar <- file.path(stan_dir, "TVP.stan")
  mod_tvpvar <- cmdstanr::cmdstan_model(stan_file_tvpvar, include_paths = stan_dir)

  fit_tvpvar <- mod_tvpvar$sample(
    data = stan_data_tvpvar,
    seed = 1234,
    init = 0.25,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 250,
    iter_sampling = 250
  )

  # ---- Extract forecasts (ALR -> simplex) for TRAIN ----
  Y_hat_tvpvar_train <- fit_tvpvar$draws("Y_hat_train", format = "draws_df")
  Y_hat_long_raw_tvpvar_train <-
    Y_hat_tvpvar_train |>
    dplyr::select(.draw, starts_with("Y_hat")) |>
    tidyr::pivot_longer(starts_with("Y_hat")) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",")

  # exponentiate ALR and back-transform to simplex (add base component)
  Y_hat_long_raw_tvpvar_train_exp <- Y_hat_long_raw_tvpvar_train
  Y_hat_long_raw_tvpvar_train_exp$value <- exp(Y_hat_long_raw_tvpvar_train_exp$value)
  Y_hat_long_raw_tvpvar_train_exp <- Y_hat_long_raw_tvpvar_train_exp[, c(1,3,4,5)]

  Y_hat_wide_tvpvar_train <- Y_hat_long_raw_tvpvar_train_exp %>%
    tidyr::pivot_wider(
      names_from = currency,
      values_from = value,
      names_prefix = "currency_",
      values_fill = list(value = 0)
    ) %>%
    dplyr::rename_with(~ paste0(.x, "_value"), starts_with("currency_"))

  Y_hat_transformed_tvpvar_train <- Y_hat_wide_tvpvar_train %>%
    dplyr::group_by(.draw, t) %>%
    dplyr::mutate(across(starts_with("currency_"),
                         ~ .x / (1 + sum(c_across(starts_with("currency_")))))) %>%
    dplyr::ungroup()

  Y_hat_transformed_tvpvar_train$currency_0_value <-
    1 - (Y_hat_transformed_tvpvar_train$currency_1_value +
           Y_hat_transformed_tvpvar_train$currency_2_value +
           Y_hat_transformed_tvpvar_train$currency_3_value +
           Y_hat_transformed_tvpvar_train$currency_4_value)

  Y_hat_long_tvpvar_train <- Y_hat_transformed_tvpvar_train %>%
    tidyr::pivot_longer(cols = starts_with("currency_"), names_to = "currency", values_to = "value") %>%
    dplyr::mutate(currency = sub("currency_", "", currency),
                  currency = sub("_value", "", currency)) %>%
    dplyr::group_by(t, currency) %>%
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      qlo  = quantile(value, 0.05, na.rm = TRUE, names = FALSE),
      qhi  = quantile(value, 0.95, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )

  # reindex currencies to match your 1..5 labels
  Y_hat_long_tvpvar_train$currency <- as.numeric(Y_hat_long_tvpvar_train$currency) + 1
  Y_hat_long_tvpvar_train$currency <- paste0(Y_hat_long_tvpvar_train$currency)
  Y_hat_long_tvpvar_train$forecast <- "tvpvar"
  Y_hat_long_tvpvar_train$t <- as.numeric(Y_hat_long_tvpvar_train$t)

  # ---- Extract forecasts (ALR -> simplex) for TEST ----
  Y_hat_tvpvar <- fit_tvpvar$draws("Y_hat", format = "draws_df")
  Y_hat_long_raw_tvpvar <-
    Y_hat_tvpvar |>
    dplyr::select(.draw, starts_with("Y_hat")) |>
    tidyr::pivot_longer(starts_with("Y_hat")) |>
    dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
    tidyr::separate(idx, c("t", compositional_var), sep = ",")

  Y_hat_long_raw_tvpvar_exp <- Y_hat_long_raw_tvpvar
  Y_hat_long_raw_tvpvar_exp$value <- exp(Y_hat_long_raw_tvpvar_exp$value)
  Y_hat_long_raw_tvpvar_exp <- Y_hat_long_raw_tvpvar_exp[, c(1,3,4,5)]

  Y_hat_wide_tvpvar <- Y_hat_long_raw_tvpvar_exp %>%
    tidyr::pivot_wider(
      names_from = currency,
      values_from = value,
      names_prefix = "currency_",
      values_fill = list(value = 0)
    ) %>%
    dplyr::rename_with(~ paste0(.x, "_value"), starts_with("currency_"))

  Y_hat_transformed_tvpvar <- Y_hat_wide_tvpvar %>%
    dplyr::group_by(.draw, t) %>%
    dplyr::mutate(across(starts_with("currency_"),
                         ~ .x / (1 + sum(c_across(starts_with("currency_")))))) %>%
    dplyr::ungroup()

  Y_hat_transformed_tvpvar$currency_0_value <-
    1 - (Y_hat_transformed_tvpvar$currency_1_value +
           Y_hat_transformed_tvpvar$currency_2_value +
           Y_hat_transformed_tvpvar$currency_3_value +
           Y_hat_transformed_tvpvar$currency_4_value)

  Y_hat_long_tvpvar <- Y_hat_transformed_tvpvar %>%
    tidyr::pivot_longer(cols = starts_with("currency_"), names_to = "currency", values_to = "value") %>%
    dplyr::mutate(currency = sub("currency_", "", currency),
                  currency = sub("_value", "", currency)) %>%
    dplyr::group_by(t, currency) %>%
    dplyr::summarize(
      mean = mean(value, na.rm = TRUE),
      qlo  = quantile(value, 0.05, na.rm = TRUE, names = FALSE),
      qhi  = quantile(value, 0.95, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )

  Y_hat_long_tvpvar$currency <- as.numeric(Y_hat_long_tvpvar$currency) + 1
  Y_hat_long_tvpvar$currency <- paste0(Y_hat_long_tvpvar$currency)
  Y_hat_long_tvpvar$forecast <- "tvpvar"
  Y_hat_long_tvpvar$t <- as.numeric(Y_hat_long_tvpvar$t)

  # ---- Join with actuals
  validation_r1_tvpvar_train <- dplyr::full_join(Y_hat_long_tvpvar_train, y_train_long)
  validation_r1_tvpvar       <- dplyr::full_join(Y_hat_long_tvpvar, y_test_long)

  # ---- Metrics on the simplex
  forecasts_tvpvar <- validation_r1_tvpvar %>%
    dplyr::filter(forecast == "tvpvar") %>%
    dplyr::select(t, currency, mean) %>%
    dplyr::rename(forecast_mean = mean)

  actuals_tvpvar <- validation_r1_tvpvar %>%
    dplyr::filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    dplyr::rename(actual_mean = mean)

  merged_df_tvpvar <- forecasts_tvpvar %>%
    dplyr::inner_join(actuals_tvpvar, by = c("t", "currency"))

  metrics_tvpvar <- merged_df_tvpvar %>%
    dplyr::mutate(squared_error  = (actual_mean - forecast_mean)^2,
                  absolute_error = abs(actual_mean - forecast_mean)) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(
      SSR  = sum(squared_error),
      MAE  = mean(absolute_error),
      mSSR = sqrt(mean(squared_error)),
      TAE  = sum(absolute_error)
    )

  # ---- Squared standardized residuals on ALR (training) for PACF
  # Build ALR of forecasts/actuals on TRAIN window
  forecasts_tvpvar_train <- validation_r1_tvpvar_train %>%
    dplyr::filter(forecast == "tvpvar") %>%
    dplyr::select(t, currency, mean) %>%
    dplyr::rename(forecast_mean = mean)

  actuals_tvpvar_train <- validation_r1_tvpvar_train %>%
    dplyr::filter(forecast == "actuals") %>%
    dplyr::select(t, currency, mean) %>%
    dplyr::rename(actual_mean = mean)

  alr_forecast_tvpvar_train <- forecasts_tvpvar_train %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(alr_value_forecast = dplyr::if_else(
      currency != "1",
      log(forecast_mean / forecast_mean[currency == "1"]),
      NA_real_
    )) %>% dplyr::ungroup() %>% tidyr::drop_na()

  alr_actual_tvpvar_train <- actuals_tvpvar_train %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(alr_value_actual = dplyr::if_else(
      currency != "1",
      log(actual_mean / actual_mean[currency == "1"]),
      NA_real_
    )) %>% dplyr::ungroup() %>% tidyr::drop_na()

  merged_alr_tvpvar_train <- alr_forecast_tvpvar_train %>%
    dplyr::inner_join(alr_actual_tvpvar_train, by = c("t", "currency"))

  # Posterior mean Sigma on ALR from TVP-VAR
  sigma_samples_tvpvar   <- fit_tvpvar$draws("sigma", format = "draws_matrix")
  L_Omega_samples_tvpvar <- fit_tvpvar$draws("L_Omega", format = "draws_array")
  iters <- dim(L_Omega_samples_tvpvar)[1]; chs <- dim(L_Omega_samples_tvpvar)[2]
  L_Omega_mean_tvpvar <- array(L_Omega_samples_tvpvar, dim = c(iters * chs, C - 1, C - 1)) |>
    apply(c(2, 3), mean)
  sigma_mean_tvpvar <- colMeans(sigma_samples_tvpvar)
  Omega_mean_tvpvar <- L_Omega_mean_tvpvar %*% t(L_Omega_mean_tvpvar)
  Sigma_mean_tvpvar <- diag(sigma_mean_tvpvar) %*% Omega_mean_tvpvar %*% diag(sigma_mean_tvpvar)
  L_Sigma_tvpvar    <- t(chol(Sigma_mean_tvpvar))

  # Standardize ALR residuals and aggregate per t
  merged_alr_tvpvar_train <- merged_alr_tvpvar_train %>%
    dplyr::mutate(residual = alr_value_forecast - alr_value_actual) %>%
    dplyr::arrange(t, currency)

  grouped_residuals_tvpvar <- merged_alr_tvpvar_train %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(residuals = list(residual), .groups = "drop") %>%
    dplyr::mutate(
      standardized_residuals = purrr::map(residuals, ~ as.numeric(solve(L_Sigma_tvpvar, .x))),
      squared_error          = purrr::map_dbl(standardized_residuals, ~ sum(.x^2))
    ) %>%
    dplyr::select(t, squared_error)

  # PACF of squared standardized residuals
  pacf_results_tvpvar <- grouped_residuals_tvpvar %>%
    dplyr::do({
      se <- .$squared_error
      pc <- pacf(se, plot = FALSE)
      tibble::tibble(lag = 1:length(pc$acf), pacf_value = pc$acf)
    })

  # ---- Store results to match the rest
  pacf_results_tvpvar_list[[i]] <- tibble::tibble(lag = pacf_results_tvpvar$lag,
                                                  pacf_value = pacf_results_tvpvar$pacf_value)
  results_tvpvar[[i]] <- tibble::tibble(simulation = i, value = metrics_tvpvar)


  pacf_results_bdarma_list[[i]] <- tibble(lag = pacf_results_bdarma$lag, pacf_value = pacf_results_bdarma$pacf_value)

  pacf_results_bdarch_list[[i]] <- tibble(lag = pacf_results_bdarch_series$lag, pacf_value = pacf_results_bdarch_series$pacf_value)
  pacf_results_btvarma_list[[i]] <- tibble(lag = pacf_results_btvama_series$lag, pacf_value = pacf_results_btvama_series$pacf_value)


  results_btvarma[[i]] <- tibble(simulation = i, value = metrics_varma)
  results_bdarma[[i]] <- tibble(simulation = i, value = metrics_bdarma)
  results_bdarch[[i]] <- tibble(simulation = i, value = metrics_bdarch)

  print(metrics_varma)

  print(metrics_tvpvar)


  print(metrics_bdarma)
  print(metrics_bdarch)

  print(i)
}

results_tvpvar

combined_metrics <- bind_rows(
  bind_rows(results_btvarma, .id = "simulation") %>% mutate(model = "B-tVARMA"),
  bind_rows(results_tvpvar, .id = "simulation") %>% mutate(model = "B-tvpvar"),
  bind_rows(results_bdarch, .id = "simulation") %>% mutate(model = "B-DARMA-DARCH"),
  bind_rows(results_bdarma, .id = "simulation") %>% mutate(model = "B-DARMA")
)
saveRDS(combined_metrics,file="summ_varma_shock.RDS")


combined_metrics

summary_df <- combined_metrics %>%
  unnest(value) %>%
  group_by(model) %>%
  summarize(
    avg_SSR =100* mean(`SSR`, na.rm = TRUE),
    avg_MAE =100* mean(`MAE`, na.rm = TRUE),
    avg_mSSR =100* mean(`mSSR`, na.rm = TRUE),
    avg_TAE = 100*mean(`TAE`, na.rm = TRUE)
  )

# Print the summarized dataframe
print(summary_df,digits=10)
i



summary_df_tvarma_shock <- summary_df
library(purrr)




count_pacf_threshold <- function(tibble, threshold) {
  tibble %>%
    summarise(count = sum(abs(pacf_value) > threshold)) %>%
    pull(count) # Extract the numeric count from the summarise result
}


compute_counts <- function(pacf_list, thresholds) {
  map_dfr(thresholds, function(threshold) {
    counts <- map_dbl(pacf_list, count_pacf_threshold, threshold = threshold)
    tibble(threshold = threshold, sum_count = sum(counts))
  })
}



# Apply the function to each element in the list
thresholds <- c(.2,0.25, 0.3)

# Compute counts for each model and store in a data frame
pacf_counts_bdarma_df <- compute_counts(pacf_results_bdarma_list, thresholds)
pacf_counts_bdarch_df <- compute_counts(pacf_results_bdarch_list, thresholds)

pacf_counts_bdarch_df
pacf_counts_btvarma_df <- compute_counts(pacf_results_btvarma_list, thresholds)

pacf_counts_tvpvar_df <- compute_counts(pacf_results_tvpvar_list, thresholds)


# Combine the results into a single data frame
final_counts_df <- bind_rows(
  pacf_counts_bdarma_df %>% mutate(model = "BDARMA"),
  pacf_counts_bdarch_df %>% mutate(model = "BDARCH"),
  pacf_counts_btvarma_df %>% mutate(model = "BtVARMA"),
  pacf_counts_tvpvar_df  %>% mutate(model = "Btvpvar")
)

# View the results
final_counts_df$sum_count <- 100*final_counts_df$sum_count/(16*50)
final_counts_df


final_counts_df_tvarma_shock <- final_counts_df

combined_pacf_bdarma <- bind_rows(pacf_results_bdarma_list, .id = "list_id")

# Calculate the average pacf_value for each lag
average_pacf_bdarma <- combined_pacf_bdarma %>%
  group_by(lag) %>%
  summarize(mean_pacf_value = mean(pacf_value))



combined_pacf_bdarch <- bind_rows(pacf_results_bdarch_list, .id = "list_id")

combined_pacf_tvpvar <- bind_rows(pacf_results_tvpvar_list, .id = "list_id")

# Calculate the average pacf_value for each lag
average_pacf_bdarch <- combined_pacf_bdarch %>%
  group_by(lag) %>%
  summarize(mean_pacf_value = mean(pacf_value))

average_pacf_tvpvar <- combined_pacf_tvpvar %>%
  group_by(lag) %>%
  summarize(mean_pacf_value = mean(pacf_value))

combined_pacf_btvarma <- bind_rows(pacf_results_btvarma_list, .id = "list_id")

# Calculate the average pacf_value for each lag
average_pacf_btvarma <- combined_pacf_btvarma %>%
  group_by(lag) %>%
  summarize(mean_pacf_value = mean(pacf_value))

average_pacf_btvarma <- combined_pacf_btvarma %>%
  group_by(lag) %>%
  summarize(mean_pacf_value = mean(pacf_value))



average_pacf_btvarma <- as.data.frame(average_pacf_btvarma)
average_pacf_bdarma <- as.data.frame(average_pacf_bdarma)
average_pacf_bdarch <- as.data.frame(average_pacf_bdarch)

average_pacf_tvpvar <- as.data.frame(average_pacf_tvpvar)


average_pacf_bdarch$model <- "B-DARCH"
average_pacf_bdarma$model <- "B-DARMA"
average_pacf_btvarma$model <- "B-tVARMA"

average_pacf_tvpvar$model <- "B-tvpvar"


combined_pacf <- full_join(average_pacf_bdarch,average_pacf_bdarma)
combined_pacf <- full_join(combined_pacf,average_pacf_btvarma)
combined_pacf <- full_join(combined_pacf,average_pacf_tvpvar)


combined_pacf_tvarma_shock <- combined_pacf

saveRDS(combined_pacf_tvarma_shock,file="pacf_sim_tvarma_shock.RDS")

ggplot(combined_pacf, aes(x = lag, y = mean_pacf_value, color = model)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "Comparison of PACF Values for Different Models",
       x = "Lag",
       y = "Mean PACF Value") +
  theme_minimal()

ggplot(combined_pacf, aes(x = factor(lag), y = mean_pacf_value, fill = model)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Side-by-side bars
  labs(title = "PACF Comparison of Models",
       x = "Lag",
       y = "Mean PACF Value",
       fill = "Model") +  # Label for the fill legend
  theme_minimal() + ylim(-.2,.2)+
  theme(legend.position = "bottom")




simulation_setup <- bind_rows(
  bind_rows(A_sim, .id = "simulation")%>% mutate(parameter = "A"),
  bind_rows(beta_sim, .id = "simulation") %>% mutate(parameter = "beta"),
  bind_rows(y_sim, .id = "simulation") %>% mutate(parameter = "y_init"),
  bind_rows(phi_beta_sim, .id = "simulation") %>% mutate(parameter = "phi_beta")
)




# Function to calculate the proportion of simulations where each lag exceeds the threshold
compute_lag_proportions <- function(pacf_list, thresholds) {
  # Get the number of simulations
  num_sims <- length(pacf_list)

  # Assuming each simulation returns a tibble with pacf_value and lag columns
  map_dfr(thresholds, function(threshold) {
    # Create a data frame to store proportions for each lag
    lag_proportions <- map_dfr(1:nrow(pacf_list[[1]]), function(lag_index) {
      # Count the number of simulations where the PACF value for this lag exceeds the threshold
      exceed_count <- sum(map_lgl(pacf_list, function(tibble) {
        abs(tibble$pacf_value[lag_index]) > threshold
      }))

      # Compute proportion of simulations where this lag exceeds the threshold
      proportion <- exceed_count / num_sims

      # Return a tibble with the lag, threshold, and proportion
      tibble(lag = lag_index, threshold = threshold, proportion = proportion)
    })

    lag_proportions
  })
}

# Define your thresholds
thresholds <- c( 0.2, 0.25, 0.3, 0.35)

# Compute proportions for each model
pacf_proportions_bdarma <- compute_lag_proportions(pacf_results_bdarma_list, thresholds)
pacf_proportions_bdarch <- compute_lag_proportions(pacf_results_bdarch_list, thresholds)
pacf_proportions_btvarma <- compute_lag_proportions(pacf_results_btvarma_list, thresholds)

# Add a column to indicate the model and combine results into a single data frame
final_proportions_df <- bind_rows(
  pacf_proportions_bdarma %>% mutate(model = "BDARMA"),
  pacf_proportions_bdarch %>% mutate(model = "BDARCH"),
  pacf_proportions_btvarma %>% mutate(model = "BtVARMA")
)

library(ggplot2)

# Create the ggplot for PACF proportions, faceting by lag
ggplot(final_proportions_df, aes(x = threshold, y = proportion, color = model)) +
  geom_line() +                # Draw lines connecting the points for each model
  geom_point() +               # Add points for each threshold
  facet_wrap(~ lag) +          # Separate plots for each lag
  labs(
    title = "PACF Proportions by Lag",
    x = "Threshold",
    y = "Proportion of PACF values > Threshold"
  ) +
  theme_minimal() +            # Use a minimal theme for clean visualization
  theme(
    strip.text = element_text(size = 10),    # Customize facet labels
    legend.position = "right"                # Show the legend on the right
  )

