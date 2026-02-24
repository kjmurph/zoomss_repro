# Runner script for OAT sensitivity analysis
# Output logged to calibration/oat_run_log.txt
setwd("C:/Users/kjmurphy/OneDrive - The University of Queensland/Documents/GitHub/zoomss_repro")
source("calibration/06_fish_param_sensitivity.R")
results <- run_full_sensitivity(run_mc = FALSE)
cat("\nOAT run complete. Results saved to calibration/fish_sensitivity_full_results.rds\n")
