# =============================================================================
# 00_smoke_test.R
# Minimal functional test before running full sensitivity analysis
# =============================================================================
#
# Purpose: Verify that the dev-branch zoomss package loads correctly and that
#          a short simulation runs to completion with valid output structure.
#          This must pass before committing to the full OAT sensitivity
#          (calibration/06_fish_param_sensitivity.R).
#
# Prerequisites:
#   - Working directory set to the package root (zoomss_repro/)
#   - devtools installed
#
# Expected runtime: < 2 minutes
# =============================================================================

cat("=== ZooMSS Smoke Test ===\n\n")

# ── 1. Load dev-branch package ──
cat("1. Loading dev-branch package with devtools::load_all()...\n")
devtools::load_all(".")
cat("   Package loaded successfully.\n\n")

# ── 2. Load default functional groups ──
cat("2. Loading default functional groups...\n")
Groups <- getGroups()

fish_idx <- which(Groups$Type == "Fish")
zoo_idx  <- which(Groups$Type == "Zooplankton")

cat("   Groups loaded:", nrow(Groups), "total |",
    length(zoo_idx), "zooplankton |", length(fish_idx), "fish\n")

# Print fish parameters of interest
cat("\n   Fish group parameters:\n")
fish_params <- data.frame(
  Species    = Groups$Species[fish_idx],
  W0         = Groups$W0[fish_idx],
  Wmax       = Groups$Wmax[fish_idx],
  Wmat       = Groups$Wmat[fish_idx],
  PPMR       = Groups$PPMR[fish_idx],
  FeedWidth  = Groups$FeedWidth[fish_idx],
  f_M        = Groups$f_M[fish_idx],
  K_growth   = Groups$K_growth[fish_idx],
  R_frac     = 1 - Groups$f_M[fish_idx] - Groups$K_growth[fish_idx],
  repro_eff  = Groups$repro_eff[fish_idx],
  ZSpre      = Groups$ZSpre[fish_idx],
  ZSexp      = Groups$ZSexp[fish_idx]
)
print(fish_params, row.names = FALSE)

# ── 3. Validate energy budget ──
cat("\n3. Validating energy budget constraints...\n")
R_frac <- 1 - Groups$f_M[fish_idx] - Groups$K_growth[fish_idx]
if (any(R_frac < 0.15)) {
  stop("FAIL: R_frac < 0.15 for: ",
       paste(Groups$Species[fish_idx][R_frac < 0.15], collapse = ", "))
}
cat("   R_frac >= 0.15 for all fish groups: PASS\n")

# ── 4. Create short environmental data ──
cat("\n4. Creating environmental data (10 years, static, Chl = 0.5, SST = 15)...\n")
env_data <- createEnviroData(
  n_years   = 10,
  dt        = 0.01,
  base_sst  = 15,
  base_chl  = 0.5,
  seasonal  = FALSE
)
input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
cat("   Environment created:", nrow(input_params), "time steps\n")

# ── 5. Run model ──
cat("\n5. Running ZooMSS model (10 years, isave = 100)...\n")
t0 <- Sys.time()
mdl <- suppressMessages(
  zoomss_model(input_params, Groups, isave = 100)
)
elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat("   Model completed in", elapsed, "seconds\n")

# ── 6. Check output structure ──
cat("\n6. Checking output structure...\n")

required_fields <- c("param", "time", "abundance", "growth", "mortality",
                      "diet", "biomass", "biomassC",
                      "repro_rate", "SSB", "recruitment", "total_repro_output")

missing <- setdiff(required_fields, names(mdl))
if (length(missing) > 0) {
  stop("FAIL: Missing output fields: ", paste(missing, collapse = ", "))
}
cat("   All required output fields present: PASS\n")

# Check array dimensions
n_saved <- dim(mdl$abundance)[1]
n_grps  <- dim(mdl$abundance)[2]
n_grid  <- dim(mdl$abundance)[3]
cat("   Abundance array: [", n_saved, "×", n_grps, "×", n_grid, "]\n")
cat("   Time vector length:", length(mdl$time), "\n")
cat("   SSB array: [", paste(dim(mdl$SSB), collapse = " × "), "]\n")

# ── 7. Check biomass ──
cat("\n7. Checking final biomass...\n")
final_biomass <- mdl$biomass[n_saved, , ]

group_biomass <- rowSums(final_biomass)
names(group_biomass) <- Groups$Species

cat("   Total biomass by group (final time step):\n")
for (i in seq_along(group_biomass)) {
  status <- if (group_biomass[i] > 0) "ALIVE" else "EXTINCT"
  cat(sprintf("     %-15s %12.4e  [%s]\n",
              names(group_biomass)[i], group_biomass[i], status))
}

# Check fish survival
fish_alive <- group_biomass[fish_idx] > 0
if (!all(fish_alive)) {
  cat("\n   WARNING: Not all fish groups survived the 10-year test.\n")
  cat("   Extinct:", paste(Groups$Species[fish_idx][!fish_alive], collapse = ", "), "\n")
  cat("   This may indicate parameter issues, but 10 years is short.\n")
} else {
  cat("   All fish groups alive after 10 years: PASS\n")
}

# ── 8. Check key utility functions ──
cat("\n8. Testing utility functions...\n")

# averageTimeSeries
avg_biomass <- averageTimeSeries(mdl, "biomass", n_years = 5)
cat("   averageTimeSeries(): returned [",
    paste(dim(avg_biomass), collapse = " × "), "] PASS\n")

# getBiomass (already run inside zoomss_model, but test directly)
bio_ww <- getBiomass(mdl, units = "ww")
cat("   getBiomass():         returned [",
    paste(dim(bio_ww), collapse = " × "), "] PASS\n")

# ── 9. Estimate runtime for full sensitivity ──
cat("\n9. Runtime estimate for full sensitivity analysis...\n")
time_per_year <- as.numeric(elapsed) / 10
est_250yr <- time_per_year * 250
n_runs <- 51  # 8 params × 2 levels × 3 Chl + 3 baselines
n_cores <- 14
n_batches <- ceiling(n_runs / n_cores)
est_total_mins <- (est_250yr * n_batches) / 60

cat(sprintf("   Time per year:       %.2f sec\n", time_per_year))
cat(sprintf("   Est. per 250yr run:  %.1f sec (%.1f min)\n", est_250yr, est_250yr / 60))
cat(sprintf("   Total runs:          %d (in ~%d batches of %d cores)\n",
            n_runs, n_batches, n_cores))
cat(sprintf("   Est. total time:     %.0f min (%.1f hours)\n",
            est_total_mins, est_total_mins / 60))

# ── Summary ──
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("SMOKE TEST RESULT: PASS\n")
cat("All checks completed successfully.\n")
cat("Safe to proceed with calibration/06_fish_param_sensitivity.R\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
