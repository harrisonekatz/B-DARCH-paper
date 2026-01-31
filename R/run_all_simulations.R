#!/usr/bin/env Rscript
# =============================================================================
# run_all_simulations.R
# Master script to run all simulation studies and produce Tables 1-2, Figures S1-S2
# =============================================================================

# Setup
library(here)

# Create output directory
if (!dir.exists("output")) dir.create("output")

cat("=== B-DARCH Paper Reproducibility Package ===\n\n")

# -----------------------------------------------------------------------------
# Run Simulation Studies
# -----------------------------------------------------------------------------

cat("Running Simulation Study 1: DARMA DGP with Shocks...\n")
source("R/simulation_study_shock_DARMA.R")
cat("  -> Output: summ_darma_shock.RDS, pacf_sim_darma_shock.RDS\n\n")

cat("Running Simulation Study 2: tVARMA DGP with Shocks...\n")
source("R/simulation_tvarma_shock.R")
cat("  -> Output: summ_varma_shock.RDS, pacf_sim_tvarma_shock.RDS\n\n")
cat("Running Simulation Study 3: DARMA DGP with Regime Shifts...\n")
source("R/simulation_study_bdarma_regime_shift.R")
cat("  -> Output: summ_darma_regime_shift.RDS, pacf_sim_darma_regime_shift.RDS\n\n")

cat("Running Simulation Study 4: tVARMA DGP with Regime Shifts...\n")
source("R/simulation_tvarma_regime_shift.R")
cat("  -> Output: summ_varma_regime_shift.RDS, pacf_sim_tvarma_regime_shift.RDS\n\n")

# -----------------------------------------------------------------------------
# Generate Tables and Figures
# -----------------------------------------------------------------------------

cat("Generating Tables 1-2 and Figures S1-S2...\n")
source("R/generate_tables_figures.R")
table1 <- generate_table1()
table2 <- generate_table2()
generate_figure_s1()
generate_figure_s2()
generate_latex_tables(table1, table2)
cat("\n=== All outputs generated successfully ===\n")
cat("Check the output/ directory for:\n")
cat("  - table1.csv (Table 1: Shock simulations)\n")
cat("  - table2.csv (Table 2: Regime shift simulations)\n")
cat("  - figure_s1.pdf (Figure S1: PACF for shock simulations)\n")
cat("  - figure_s2.pdf (Figure S2: PACF for regime shift simulations)\n")




