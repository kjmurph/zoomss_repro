# =============================================================================
# 04_run_calibration.R
# Run optimisation to calibrate revised ZooMSS against original baseline
# =============================================================================
#
# Purpose: Use numerical optimisation to find parameter values that minimise
#          the difference between revised and original ZooMSS outputs across
#          a chlorophyll gradient.
#
# Prerequisites:
#   - Baseline from 01_generate_baseline.R
#   - Sensitivity results from 02_sensitivity_analysis.R (to inform parameter selection)
#   - Objective function from 03_objective_function.R
#   - Revised zoomss package installed (devtools::load_all("."))
#
# Output: calibration/calibration_results.rds
# =============================================================================

library(zoomss)
library(future.apply)

# ── Set up parallel execution across chlorophyll levels ──
plan(multisession, workers = parallelly::availableCores() - 1)
cat("Using", parallelly::availableCores() - 1, "parallel workers\n")

# ── Load baseline from original model ──
baseline <- readRDS("calibration/baseline_original_zoomss.rds")
chl_levels <- baseline$chl_levels

# ── Load sensitivity analysis results to inform parameter selection ──
if (file.exists("calibration/sensitivity_summary.rds")) {
  sensitivity <- readRDS("calibration/sensitivity_summary.rds")
  cat("Sensitivity ranking (by zoo proportions):\n")
  print(sensitivity[order(-sensitivity$delta_zoo_prop), ])
  cat("\nUse this to decide which parameters to include in optimisation.\n")
  cat("Parameters with negligible delta_zoo_prop can be fixed at defaults.\n\n")
} else {
  cat("WARNING: No sensitivity results found. Run 02_sensitivity_analysis.R first.\n")
  cat("Proceeding with all candidate parameters.\n\n")
}

# ── Load default Groups as starting point ──
Groups <- getGroups()

# ── Source the objective function ──
source("calibration/03_objective_function.R")

# ── Define parameter bounds ──
# These bounds enforce the energy budget constraint f_M + K_growth <= 1
par_lower <- c(
  f_M                = 0.30,   # Minimum metabolic fraction
  K_growth_zoo_base  = 0.15,   # Minimum zoo growth fraction
  K_growth_fish      = 0.15,   # Minimum fish growth fraction
  repro_eff          = 1e-6    # Minimum reproductive efficiency
)

par_upper <- c(
  f_M                = 0.70,   # Maximum metabolic fraction
  K_growth_zoo_base  = 0.50,   # Maximum zoo growth fraction
  K_growth_fish      = 0.45,   # Maximum fish growth fraction
  repro_eff          = 0.01    # Maximum reproductive efficiency
)

# Starting values (from current defaults)
par_start <- c(
  f_M                = 0.50,
  K_growth_zoo_base  = 0.35,
  K_growth_fish      = 0.30,
  repro_eff          = 0.002
)

cat("Parameter bounds:\n")
print(data.frame(
  parameter = names(par_start),
  lower     = par_lower,
  start     = par_start,
  upper     = par_upper
))

# ── Test objective function at starting values ──
cat("\nTesting objective function at starting values...\n")
obj_start <- calibration_objective(par_start, baseline, chl_levels, Groups, verbose = TRUE)
cat("Starting objective:", obj_start, "\n\n")

# =============================================================================
# Option A: Nelder-Mead (derivative-free, good for noisy objectives)
# =============================================================================
cat("=== Running Nelder-Mead optimisation ===\n")
t_start <- Sys.time()

result_nm <- optim(
  par     = par_start,
  fn      = calibration_objective,
  method  = "Nelder-Mead",
  baseline   = baseline,
  chl_levels = chl_levels,
  Groups     = Groups,
  verbose    = TRUE,
  control = list(
    maxit  = 200,
    trace  = 1,
    reltol = 1e-4
  )
)

t_nm <- difftime(Sys.time(), t_start, units = "mins")
cat(sprintf("\nNelder-Mead complete in %.1f minutes\n", as.numeric(t_nm)))
cat("Best objective:", result_nm$value, "\n")
cat("Best parameters:\n")
print(result_nm$par)

# =============================================================================
# Option B: L-BFGS-B (box-constrained, more efficient if objective is smooth)
# =============================================================================
cat("\n=== Running L-BFGS-B optimisation ===\n")
t_start <- Sys.time()

result_bfgs <- optim(
  par     = par_start,
  fn      = calibration_objective,
  method  = "L-BFGS-B",
  lower   = par_lower,
  upper   = par_upper,
  baseline   = baseline,
  chl_levels = chl_levels,
  Groups     = Groups,
  verbose    = TRUE,
  control = list(
    maxit = 100,
    trace = 1
  )
)

t_bfgs <- difftime(Sys.time(), t_start, units = "mins")
cat(sprintf("\nL-BFGS-B complete in %.1f minutes\n", as.numeric(t_bfgs)))
cat("Best objective:", result_bfgs$value, "\n")
cat("Best parameters:\n")
print(result_bfgs$par)

# =============================================================================
# Option C: DEoptim (global optimiser, avoids local minima)
# Recommended for initial exploration because the objective landscape may be
# multimodal (multiple parameter combinations produce similar patterns)
# =============================================================================
result_de <- NULL
if (requireNamespace("DEoptim", quietly = TRUE)) {
  cat("\n=== Running DEoptim (global optimisation) ===\n")
  cat("This may take 2-4 hours with parallelised chl evaluations.\n")
  t_start <- Sys.time()

  result_de <- DEoptim::DEoptim(
    fn    = calibration_objective,
    lower = par_lower,
    upper = par_upper,
    baseline   = baseline,
    chl_levels = chl_levels,
    Groups     = Groups,
    verbose    = TRUE,
    control = DEoptim::DEoptim.control(
      NP       = 40,   # Population size (10x number of parameters)
      itermax  = 100,  # Maximum iterations
      trace    = 10,   # Print every 10 iterations
      parallelType = 0 # No internal parallelism (we parallelise across chl)
    )
  )

  t_de <- difftime(Sys.time(), t_start, units = "mins")
  cat(sprintf("\nDEoptim complete in %.1f minutes\n", as.numeric(t_de)))
  cat("Best objective:", result_de$optim$bestval, "\n")
  cat("Best parameters:\n")
  print(result_de$optim$bestmem)
} else {
  cat("\nDEoptim package not installed. Skipping global optimisation.\n")
  cat("Install with: install.packages('DEoptim')\n")
}

# =============================================================================
# Compare results across methods
# =============================================================================
cat("\n=== OPTIMISATION SUMMARY ===\n\n")

comparison <- data.frame(
  method     = c("Start", "Nelder-Mead", "L-BFGS-B"),
  objective  = c(obj_start, result_nm$value, result_bfgs$value),
  f_M        = c(par_start["f_M"], result_nm$par["f_M"], result_bfgs$par["f_M"]),
  K_zoo      = c(par_start["K_growth_zoo_base"], result_nm$par["K_growth_zoo_base"],
                 result_bfgs$par["K_growth_zoo_base"]),
  K_fish     = c(par_start["K_growth_fish"], result_nm$par["K_growth_fish"],
                 result_bfgs$par["K_growth_fish"]),
  repro_eff  = c(par_start["repro_eff"], result_nm$par["repro_eff"],
                 result_bfgs$par["repro_eff"])
)

if (!is.null(result_de)) {
  comparison <- rbind(comparison, data.frame(
    method    = "DEoptim",
    objective = result_de$optim$bestval,
    f_M       = result_de$optim$bestmem["f_M"],
    K_zoo     = result_de$optim$bestmem["K_growth_zoo_base"],
    K_fish    = result_de$optim$bestmem["K_growth_fish"],
    repro_eff = result_de$optim$bestmem["repro_eff"]
  ))
}

print(comparison)

# ── Save all results ──
saveRDS(list(
  nelder_mead = result_nm,
  lbfgsb      = result_bfgs,
  deoptim     = result_de,
  baseline    = baseline,
  par_bounds  = list(lower = par_lower, upper = par_upper),
  par_start   = par_start,
  comparison  = comparison
), "calibration/calibration_results.rds")

cat("\nResults saved to calibration/calibration_results.rds\n")
cat("Proceed to 05_evaluate_calibration.R for diagnostic plots.\n")
