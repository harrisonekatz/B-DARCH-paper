# B-DARCH: Bayesian Dirichlet Auto-Regressive Conditional Heteroskedasticity Model

Reproducibility package for "A Bayesian Dirichlet Auto-Regressive Conditional Heteroskedasticity Model for Forecasting Currency Shares" (International Journal of Forecasting, forthcoming).

**Repository:** https://github.com/harrisonekatz/B-DARCH-paper  
**Contact:** Harrison Katz (harrison.katz@airbnb.com)  
**Package assembled:** January 2026

---

## Data Availability

**Simulation study (Section 3):** Fully reproducible. All DGPs are specified in code and synthetic data is generated programmatically.

**Application (Section 4):** The Airbnb currency share data are proprietary and cannot be shared due to legal and compliance restrictions. The simulation study serves as the public, reproducible evidence base, isolating the same phenomena (mean drift, volatility clustering, regime shifts) under controlled conditions.

---

## Repository Structure

```
B-DARCH-paper/
├── R/
│   ├── compositional.R                         # Helper functions (ALR/ILR transforms, softmax, Fourier terms)
│   ├── simulation_study_shock_DARMA.R          # Sim 1: Dirichlet DGP with misreported observations
│   ├── simulation_study_bdarma_regime_shift.R  # Sim 3: Dirichlet DGP with temporary regime shifts
│   ├── simulation_tvarma_shock.R               # Sim 2: Transformed-Gaussian DGP with shocks
│   ├── simulation_tvarma_regime_shift.R        # Sim 4: Transformed-Gaussian DGP with regime shifts
│   ├── generate_tables_figures.R               # Generates Tables 1-2, Figures S1-S2 from RDS outputs
│   └── run_all_simulations.R                   # Master script to run everything
├── stan/
│   ├── bdarma.stan              # B-DARMA model (Dirichlet, constant precision)
│   ├── bdarch_pq.stan           # B-DARCH(P,Q) model (Dirichlet, time-varying precision)
│   ├── bdarch_pq_sim_trend.stan # B-DARCH variant for simulation with trend
│   ├── varma.stan               # B-tVARMA (Gaussian on ALR scale)
│   ├── TVP.stan                 # B-TVP-tVARMA (time-varying parameters)
│   ├── alr.stan                 # ALR transform functions (included by other models)
│   └── component_product.stan   # Helper functions (included by other models)
├── output/                      # Generated tables and figures (after running scripts)
└── README.md
```

---

## Computing Environment

| Component | Version |
|-----------|---------|
| R | 4.3.x or later |
| Stan / CmdStan | 2.33+ |
| cmdstanr | 0.6+ |
| Operating System | macOS / Linux / Windows |

### Required R Packages

```r
install.packages(c(
  "cmdstanr",
  "brms",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "patchwork",
  "MASS",
  "compositions",
  "checkmate",
  "tibble",
  "lubridate",
  "tidyselect",
  "rlang"
))

# Install CmdStan (if not already installed)
cmdstanr::install_cmdstan()
```

---

## Reproducing Results

### Quick Start (Full Reproduction)

```r
# From the repository root, run the master script
source("R/run_all_simulations.R")

# This will:
# 1. Run all four simulation studies (Sims 1-4)
# 2. Generate Tables 1-2 as CSV and LaTeX
# 3. Generate Figures S1-S2 as PDF
# Output files will be in the output/ directory
```

**Warning:** Full reproduction takes approximately 20-30 hours. See "Reduced Runtime" below for faster testing.

### Simulation Study Scripts

Each simulation script runs 50 replications across four models (B-DARCH, B-DARMA, B-tVARMA, B-TVP-tVARMA).

| Script | Paper Section | Output Files |
|--------|---------------|--------------|
| `simulation_study_shock_DARMA.R` | Sim 1 (Table 1, left) | `summ_darma_shock.RDS`, `pacf_sim_darma_shock.RDS` |
| `simulation_tvarma_shock.R` | Sim 2 (Table 1, right) | `summ_varma_shock.RDS`, `pacf_sim_tvarma_shock.RDS` |
| `simulation_study_bdarma_regime_shift.R` | Sim 3 (Table 2, left) | `summ_darma_regime_shift.RDS`, `pacf_sim_darma_regime_shift.RDS` |
| `simulation_tvarma_regime_shift.R` | Sim 4 (Table 2, right) | `summ_varma_regime_shift.RDS`, `pacf_sim_tvarma_regime_shift.RDS` |

### Generating Tables and Figures

After running simulations (or using saved RDS files):

```r
source("R/generate_tables_figures.R")

# Or call individual functions:
table1 <- generate_table1()    # Table 1: Shock simulations
table2 <- generate_table2()    # Table 2: Regime shift simulations
fig_s1 <- generate_figure_s1() # Figure S1: PACF for shocks
fig_s2 <- generate_figure_s2() # Figure S2: PACF for regime shifts
```

### Output Mapping

| Paper Output | Generated File | Source Script |
|--------------|----------------|---------------|
| Table 1 (Shock simulations) | `output/table1.csv`, `output/table1.tex` | `generate_tables_figures.R` |
| Table 2 (Regime shift simulations) | `output/table2.csv`, `output/table2.tex` | `generate_tables_figures.R` |
| Figure S1 (PACF - shocks) | `output/figure_s1.pdf` | `generate_tables_figures.R` |
| Figure S2 (PACF - regime shifts) | `output/figure_s2.pdf` | `generate_tables_figures.R` |

### Reduced Runtime for Testing

To verify the code works without running all 50 replications, modify `n_sim` at the top of each simulation script:

```r
n_sim <- 5  # instead of 50 (reduces runtime by ~90%)
```

---

## Runtime Expectations

Tested on Apple M2 Pro (12-core), 32GB RAM.

Each simulation script runs 50 replications. Per-replication timing:

| Model | Iterations | Approx. Time per Rep |
|-------|------------|---------------------|
| B-DARMA | 4 chains × 500 (250 warmup) | ~2-3 min |
| B-DARCH | 4 chains × 500 (250 warmup) | ~5-8 min |
| B-tVARMA | 4 chains × 500 (250 warmup) | ~2-3 min |
| B-TVP-tVARMA | 4 chains × 500 (250 warmup) | ~8-12 min |

**Total runtime per simulation script:** approximately 8-15 hours (50 reps × 4 models).

To reduce runtime for testing, modify `n_sim` at the top of each script:
```r
n_sim <- 10  # instead of 50
```

---

## MCMC Settings

All models use:
- 4 chains in parallel (`parallel_chains = 4`)
- 250 warmup + 250 sampling iterations per chain
- Seed: 1234 for reproducibility

The paper's final results used longer runs; these settings are sufficient to verify the methodology.

---

## Stan Model Details

### Include Structure

The Stan models use `#include` directives:
- `alr.stan`: Additive log-ratio transform and inverse
- `component_product.stan`: Helper functions for componentwise operations

When compiling, CmdStan resolves includes from the `stan/` directory via the `include_paths` argument.

### Model Summary

| File | Model | Key Features |
|------|-------|--------------|
| `bdarma.stan` | B-DARMA | Dirichlet likelihood, VAR(P) on ALR scale, constant precision |
| `bdarch_pq.stan` | B-DARCH(P,Q) | Dirichlet likelihood, VARMA on ALR, ARCH dynamics on precision |
| `varma.stan` | B-tVARMA | Gaussian likelihood on ALR scale, VARMA structure |
| `TVP.stan` | B-TVP-tVARMA | Gaussian on ALR with time-varying mean coefficients |

---

## License

Code: MIT License  
Paper content: Copyright retained by authors / Elsevier

---

## Citation

```bibtex
@article{katz2026bdarch,
  title={A {B}ayesian {D}irichlet Auto-Regressive Conditional Heteroskedasticity 
         Model for Forecasting Currency Shares},
  author={Katz, Harrison E. and Weiss, Robert E.},
  journal={International Journal of Forecasting},
  year={2026},
  note={Forthcoming}
}
```
