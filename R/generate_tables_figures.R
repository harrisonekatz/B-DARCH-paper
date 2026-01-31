#!/usr/bin/env Rscript
# =============================================================================
# generate_tables_figures.R
# Generate Tables 1-2 and Figures S1-S2 from simulation study outputs
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Create output directory if it doesn't exist
if (!dir.exists("output")) dir.create("output")

# Helper function to standardize model names
standardize_model_names <- function(df) {
  df %>%
    mutate(
      model = case_when(
        grepl("DARCH", model, ignore.case = TRUE) & !grepl("DARMA", model, ignore.case = TRUE) ~ "B-DARCH",
        grepl("DARMA-DARCH", model, ignore.case = TRUE) ~ "B-DARCH",
        grepl("DARMA", model, ignore.case = TRUE) & !grepl("DARCH", model, ignore.case = TRUE) ~ "B-DARMA",
        grepl("tvpvar", model, ignore.case = TRUE) ~ "B-TVP-tVARMA",
        grepl("TVP", model, ignore.case = TRUE) ~ "B-TVP-tVARMA",
        grepl("tVARMA", model, ignore.case = TRUE) ~ "B-tVARMA",
        TRUE ~ model
      )
    )
}

# =============================================================================
# Table 1: Simulation Studies 1-2 (Shock Perturbations)
# =============================================================================

generate_table1 <- function() {
  cat("Generating Table 1...\n")
  
  # Load simulation results
  # Study 1: DARMA DGP with shocks
  if (file.exists("summ_darma_shock.RDS")) {
    sim1_metrics <- readRDS("summ_darma_shock.RDS")
  } else {
    warning("summ_darma_shock.RDS not found. Run simulation_study_shock_DARMA.R first.")
    return(NULL)
  }
  
  # Study 2: tVARMA DGP with shocks
  if (file.exists("summ_varma_shock.RDS")) {
    sim2_metrics <- readRDS("summ_varma_shock.RDS")
  } else {
    warning("summ_varma_shock.RDS not found. Run simulation_tvarma_shock.R first.")
    return(NULL)
  }
  
  # Process Study 1
  # Structure: tibble with columns simulation, model, value (nested with SSR, MAE, etc.)
  summary_sim1 <- sim1_metrics %>%
    unnest(value) %>%
    group_by(model) %>%
    summarize(
      FRMSE_sim1 = 100 * sqrt(mean(SSR, na.rm = TRUE)),
      FMAE_sim1 = 100 * mean(MAE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    standardize_model_names()
  
  # Process Study 2
  summary_sim2 <- sim2_metrics %>%
    unnest(value) %>%
    group_by(model) %>%
    summarize(
      FRMSE_sim2 = 100 * sqrt(mean(SSR, na.rm = TRUE)),
      FMAE_sim2 = 100 * mean(MAE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    standardize_model_names()
  
  # Combine into Table 1 format
  table1 <- summary_sim1 %>%
    left_join(summary_sim2, by = "model") %>%
    select(
      Model = model,
      `Sim 1 FRMSE` = FRMSE_sim1,
      `Sim 1 FMAE` = FMAE_sim1,
      `Sim 2 FRMSE` = FRMSE_sim2,
      `Sim 2 FMAE` = FMAE_sim2
    ) %>%
    arrange(match(Model, c("B-DARMA", "B-DARCH", "B-tVARMA", "B-TVP-tVARMA")))
  
  # Save as CSV
  write.csv(table1, "output/table1.csv", row.names = FALSE)
  
  # Print formatted table
  cat("\nTable 1: Summary of model performance metrics (x100) - Shock Simulations\n")
  cat("=========================================================================\n")
  print(table1, digits = 2)
  
  return(table1)
}

# =============================================================================
# Table 2: Simulation Studies 3-4 (Regime Shift Perturbations)
# =============================================================================

generate_table2 <- function() {
  cat("\nGenerating Table 2...\n")
  
  # Load simulation results
  # Study 3: DARMA DGP with regime shifts
  if (file.exists("summ_darma_regime_shift.RDS")) {
    sim3_metrics <- readRDS("summ_darma_regime_shift.RDS")
  } else {
    warning("summ_darma_regime_shift.RDS not found. Run simulation_study_bdarma_regime_shift.R first.")
    return(NULL)
  }
  
  # Study 4: tVARMA DGP with regime shifts
  if (file.exists("summ_varma_regime_shift.RDS")) {
    sim4_metrics <- readRDS("summ_varma_regime_shift.RDS")
  } else {
    warning("summ_varma_regime_shift.RDS not found. Run simulation_tvarma_regime_shift.R first.")
    return(NULL)
  }
  
  # Process Study 3
  summary_sim3 <- sim3_metrics %>%
    unnest(value) %>%
    group_by(model) %>%
    summarize(
      FRMSE_sim3 = 100 * sqrt(mean(SSR, na.rm = TRUE)),
      FMAE_sim3 = 100 * mean(MAE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    standardize_model_names()
  
  # Process Study 4
  summary_sim4 <- sim4_metrics %>%
    unnest(value) %>%
    group_by(model) %>%
    summarize(
      FRMSE_sim4 = 100 * sqrt(mean(SSR, na.rm = TRUE)),
      FMAE_sim4 = 100 * mean(MAE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    standardize_model_names()
  
  # Combine into Table 2 format
  table2 <- summary_sim3 %>%
    left_join(summary_sim4, by = "model") %>%
    select(
      Model = model,
      `Sim 3 FRMSE` = FRMSE_sim3,
      `Sim 3 FMAE` = FMAE_sim3,
      `Sim 4 FRMSE` = FRMSE_sim4,
      `Sim 4 FMAE` = FMAE_sim4
    ) %>%
    arrange(match(Model, c("B-DARMA", "B-DARCH", "B-tVARMA", "B-TVP-tVARMA")))
  
  # Save as CSV
  write.csv(table2, "output/table2.csv", row.names = FALSE)
  
  # Print formatted table
  cat("\nTable 2: Summary of model performance metrics (x100) - Regime Shift Simulations\n")
  cat("================================================================================\n")
  print(table2, digits = 2)
  
  return(table2)
}

# =============================================================================
# Figure S1: PACF for Simulation Studies 1-2 (Shock Perturbations)
# =============================================================================

generate_figure_s1 <- function() {
  cat("\nGenerating Figure S1...\n")
  
  # Load PACF results
  if (!file.exists("pacf_sim_darma_shock.RDS") || !file.exists("pacf_sim_tvarma_shock.RDS")) {
    warning("PACF files not found. Run simulation scripts first.")
    return(NULL)
  }
  
  # Structure: tibble with columns lag, mean_pacf_value, model
  pacf_darma <- readRDS("pacf_sim_darma_shock.RDS") %>%
    mutate(dgp = "DARMA (Shock)") %>%
    standardize_model_names()
  
  pacf_tvarma <- readRDS("pacf_sim_tvarma_shock.RDS") %>%
    mutate(dgp = "tVARMA (Shock)") %>%
    standardize_model_names()
  
  # Combine
  combined_pacf <- bind_rows(pacf_darma, pacf_tvarma)
  
  # Define colors matching paper
  model_colors <- c(
    "B-DARCH" = "#E69F00",
    "B-DARMA" = "#56B4E9",
    "B-tVARMA" = "#009E73",
    "B-TVP-tVARMA" = "#CC79A7"
  )
  
  # Create plot for DARMA DGP
  p1 <- ggplot(combined_pacf %>% filter(dgp == "DARMA (Shock)"), 
               aes(x = factor(lag), y = mean_pacf_value, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = model_colors, name = "Model") +
    labs(
      title = "PACF - DGP: DARMA (Shock)",
      x = "Lag",
      y = "Mean PACF Value"
    ) +
    coord_cartesian(ylim = c(-0.3, 0.3)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12)
    )
  
  # Create plot for tVARMA DGP
  p2 <- ggplot(combined_pacf %>% filter(dgp == "tVARMA (Shock)"), 
               aes(x = factor(lag), y = mean_pacf_value, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = model_colors, name = "Model") +
    labs(
      title = "PACF - DGP: tVARMA (Shock)",
      x = "Lag",
      y = "Mean PACF Value"
    ) +
    coord_cartesian(ylim = c(-0.3, 0.3)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12)
    )
  
  # Combine plots
  fig_s1 <- p1 + p2 + 
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  # Save
  ggsave("output/figure_s1.pdf", fig_s1, width = 12, height = 5)
  ggsave("output/figure_s1.png", fig_s1, width = 12, height = 5, dpi = 300)
  cat("  Saved: output/figure_s1.pdf, output/figure_s1.png\n")
  
  return(fig_s1)
}

# =============================================================================
# Figure S2: PACF for Simulation Studies 3-4 (Regime Shift Perturbations)
# =============================================================================

generate_figure_s2 <- function() {
  cat("\nGenerating Figure S2...\n")
  
  # Load PACF results
  if (!file.exists("pacf_sim_darma_regime_shift.RDS") || !file.exists("pacf_sim_tvarma_regime_shift.RDS")) {
    warning("PACF files not found. Run simulation scripts first.")
    return(NULL)
  }
  
  pacf_darma <- readRDS("pacf_sim_darma_regime_shift.RDS") %>%
    mutate(dgp = "DARMA (Regime Shift)") %>%
    standardize_model_names()
  
  pacf_tvarma <- readRDS("pacf_sim_tvarma_regime_shift.RDS") %>%
    mutate(dgp = "tVARMA (Regime Shift)") %>%
    standardize_model_names()
  
  # Combine
  combined_pacf <- bind_rows(pacf_darma, pacf_tvarma)
  
  # Define colors matching paper
  model_colors <- c(
    "B-DARCH" = "#E69F00",
    "B-DARMA" = "#56B4E9",
    "B-tVARMA" = "#009E73",
    "B-TVP-tVARMA" = "#CC79A7"
  )
  
  # Create plot for DARMA DGP
  p1 <- ggplot(combined_pacf %>% filter(dgp == "DARMA (Regime Shift)"), 
               aes(x = factor(lag), y = mean_pacf_value, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = model_colors, name = "Model") +
    labs(
      title = "PACF - DGP: DARMA (Regime Shift)",
      x = "Lag",
      y = "Mean PACF Value"
    ) +
    coord_cartesian(ylim = c(-0.3, 0.3)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12)
    )
  
  # Create plot for tVARMA DGP
  p2 <- ggplot(combined_pacf %>% filter(dgp == "tVARMA (Regime Shift)"), 
               aes(x = factor(lag), y = mean_pacf_value, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = model_colors, name = "Model") +
    labs(
      title = "PACF - DGP: tVARMA (Regime Shift)",
      x = "Lag",
      y = "Mean PACF Value"
    ) +
    coord_cartesian(ylim = c(-0.3, 0.3)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12)
    )
  
  # Combine plots
  fig_s2 <- p1 + p2 + 
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  # Save
  ggsave("output/figure_s2.pdf", fig_s2, width = 12, height = 5)
  ggsave("output/figure_s2.png", fig_s2, width = 12, height = 5, dpi = 300)
  cat("  Saved: output/figure_s2.pdf, output/figure_s2.png\n")
  
  return(fig_s2)
}

# =============================================================================
# Generate LaTeX Tables
# =============================================================================

generate_latex_tables <- function(table1, table2) {
  cat("\nGenerating LaTeX tables...\n")
  
  if (is.null(table1) || is.null(table2)) {
    warning("Cannot generate LaTeX tables: missing data.")
    return(NULL)
  }
  
  # Table 1 LaTeX
  latex1 <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{Summary of model performance metrics ($\\times 100$) on the test set ($T=40$) across 50 simulations for Simulation Studies 1--2. ",
    "FRMSE ($\\times 100$) equals the across-simulation average of RMSE $= \\sqrt{\\frac{1}{TJ}\\sum_{t=1}^{T}\\sum_{j=1}^{J}(y_{tj} - \\hat{y}_{tj})^2}$; ",
    "FMAE ($\\times 100$) equals the across-simulation average of MAE $= \\frac{1}{TJ}\\sum_{t=1}^{T}\\sum_{j=1}^{J}|y_{tj} - \\hat{y}_{tj}|$.}\n",
    "\\label{tab:sim12}\n",
    "\\begin{tabular}{lcccc}\n",
    "\\hline\n",
    "Model & \\multicolumn{2}{c}{Sim 1 (DARMA Shock)} & \\multicolumn{2}{c}{Sim 2 (tVARMA Shock)} \\\\\n",
    " & Average FRMSE & Average FMAE & Average FRMSE & Average FMAE \\\\\n",
    "\\hline\n"
  )
  
  for (i in 1:nrow(table1)) {
    latex1 <- paste0(latex1, 
                     table1$Model[i], " & ",
                     sprintf("%.2f", table1$`Sim 1 FRMSE`[i]), " & ",
                     sprintf("%.2f", table1$`Sim 1 FMAE`[i]), " & ",
                     sprintf("%.2f", table1$`Sim 2 FRMSE`[i]), " & ",
                     sprintf("%.2f", table1$`Sim 2 FMAE`[i]), " \\\\\n")
  }
  
  latex1 <- paste0(latex1, "\\hline\n\\end{tabular}\n\\end{table}\n")
  
  writeLines(latex1, "output/table1.tex")
  
  # Table 2 LaTeX
  latex2 <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{Summary of model performance metrics ($\\times 100$) on the test set ($T=40$) across 50 simulations for Simulation Studies 3--4. ",
    "FRMSE ($\\times 100$) equals the across-simulation average of RMSE $= \\sqrt{\\frac{1}{TJ}\\sum_{t=1}^{T}\\sum_{j=1}^{J}(y_{tj} - \\hat{y}_{tj})^2}$; ",
    "FMAE ($\\times 100$) equals the across-simulation average of MAE $= \\frac{1}{TJ}\\sum_{t=1}^{T}\\sum_{j=1}^{J}|y_{tj} - \\hat{y}_{tj}|$.}\n",
    "\\label{tab:sim34}\n",
    "\\begin{tabular}{lcccc}\n",
    "\\hline\n",
    "Model & \\multicolumn{2}{c}{Sim 3 (DARMA Regime Shift)} & \\multicolumn{2}{c}{Sim 4 (tVARMA Regime Shift)} \\\\\n",
    " & Average FRMSE & Average FMAE & Average FRMSE & Average FMAE \\\\\n",
    "\\hline\n"
  )
  
  for (i in 1:nrow(table2)) {
    latex2 <- paste0(latex2, 
                     table2$Model[i], " & ",
                     sprintf("%.2f", table2$`Sim 3 FRMSE`[i]), " & ",
                     sprintf("%.2f", table2$`Sim 3 FMAE`[i]), " & ",
                     sprintf("%.2f", table2$`Sim 4 FRMSE`[i]), " & ",
                     sprintf("%.2f", table2$`Sim 4 FMAE`[i]), " \\\\\n")
  }
  
  latex2 <- paste0(latex2, "\\hline\n\\end{tabular}\n\\end{table}\n")
  
  writeLines(latex2, "output/table2.tex")
  
  cat("  Saved: output/table1.tex, output/table2.tex\n")
}

# =============================================================================
# Main Execution
# =============================================================================

if (!interactive()) {
  # Running as script
  table1 <- generate_table1()
  table2 <- generate_table2()
  fig_s1 <- generate_figure_s1()
  fig_s2 <- generate_figure_s2()
  
  if (!is.null(table1) && !is.null(table2)) {
    generate_latex_tables(table1, table2)
  }
  
  cat("\n=== Table and Figure Generation Complete ===\n")
} else {
  # Interactive mode
  cat("Functions available:\n")
  cat("  generate_table1()    - Generate Table 1 (Shock simulations)\n")
  cat("  generate_table2()    - Generate Table 2 (Regime shift simulations)\n")
  cat("  generate_figure_s1() - Generate Figure S1 (PACF for shocks)\n")
  cat("  generate_figure_s2() - Generate Figure S2 (PACF for regime shifts)\n")
  cat("  generate_latex_tables(table1, table2) - Generate LaTeX versions\n")
}
