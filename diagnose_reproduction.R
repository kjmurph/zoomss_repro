# Comprehensive ZooMSS Reproduction Diagnostic Suite
# Investigate why SSB is showing as 0 and debug reproduction implementation

library(devtools)
load_all()
source("R/zoomss_reproduction.R")
source("R/update_groups_for_reproduction.R")

cat("=== COMPREHENSIVE REPRODUCTION DIAGNOSTICS ===\n")
cat("Date:", Sys.time(), "\n\n")

# Test 1: Parameter Setup Verification
cat("=== TEST 1: PARAMETER SETUP ===\n")
Groups <- getGroups()
Groups <- add_reproduction_params(Groups)

fish_grps <- which(Groups$Type == "Fish")
cat("Fish groups found:", length(fish_grps), "\n")
cat("Fish group indices:", paste(fish_grps, collapse = ", "), "\n")
cat("Fish species:", paste(Groups$Species[fish_grps], collapse = ", "), "\n")

for(i in seq_along(fish_grps)) {
  grp <- fish_grps[i]
  cat("Fish", i, "(", Groups$Species[grp], "):\n")
  cat("  W0 (egg size):", Groups$W0[grp], "\n")
  cat("  Wmat (maturity):", Groups$Wmat[grp], "\n")
  cat("  Wmax (max size):", Groups$Wmax[grp], "\n")
  cat("  ReproInvest:", Groups$ReproInvest[grp], "\n")
  cat("  ReproExp:", Groups$ReproExp[grp], "\n")
  cat("  Repro flag:", Groups$Repro[grp], "\n")
}
cat("\n")

# Test 2: Size Grid and Maturity Analysis
cat("=== TEST 2: SIZE GRID AND MATURITY ANALYSIS ===\n")

# Create test environment for longer simulation
test_years <- 5
test_time <- seq(0, test_years, length.out = test_years * 12)
test_env <- data.frame(
  time = test_time,
  sst = 15 + 3 * sin(2 * pi * test_time),
  chl = exp(-3 + 0.3 * sin(2 * pi * test_time))
)

input_params <- createInputParams(test_env$time, test_env$sst, test_env$chl)
params <- zoomss_params(Groups, input_params, isave = 6)  # Save every 6 time steps

cat("Size grid information:\n")
cat("Number of size classes:", params$ngrid, "\n")
cat("Weight range:", round(range(params$w), 6), "g\n")
cat("Log10 weight range:", round(range(log10(params$w)), 2), "\n")

# Check size alignment for fish groups
for(i in seq_along(fish_grps)) {
  grp <- fish_grps[i]
  cat("\nFish", i, "(", Groups$Species[grp], ") size analysis:\n")
  
  # Find size class indices
  egg_idx <- which.min(abs(log10(params$w) - Groups$W0[grp]))
  mat_idx <- which.min(abs(log10(params$w) - Groups$Wmat[grp]))
  max_idx <- which.min(abs(log10(params$w) - Groups$Wmax[grp]))
  
  cat("  Egg size index:", egg_idx, "(weight =", round(params$w[egg_idx], 6), "g)\n")
  cat("  Maturity index:", mat_idx, "(weight =", round(params$w[mat_idx], 3), "g)\n")
  cat("  Max size index:", max_idx, "(weight =", round(params$w[max_idx], 0), "g)\n")
  
  # Test maturity ogive
  mat_ogive <- calculate_maturity_ogive(params$w, Groups$Wmat[grp])
  mature_sizes <- which(mat_ogive > 0.5)
  cat("  Mature size classes:", length(mature_sizes), "out of", params$ngrid, "\n")
  cat("  Maturity range:", round(range(mat_ogive[mature_sizes]), 3), "\n")
  
  # Test reproductive investment
  repro_inv <- calculate_reproductive_investment(params$w, Groups$ReproInvest[grp], 
                                               Groups$ReproExp[grp], Groups$Wmat[grp])
  investing_sizes <- which(repro_inv > 0)
  cat("  Size classes with repro investment:", length(investing_sizes), "\n")
  if(length(investing_sizes) > 0) {
    cat("  Repro investment range:", round(range(repro_inv[investing_sizes]), 4), "\n")
  }
}
cat("\n")

# Test 3: Initial Abundance Analysis  
cat("=== TEST 3: INITIAL ABUNDANCE ANALYSIS ===\n")
model <- zoomss_setup(params)

cat("Initial abundance dimensions:", paste(dim(model$N), collapse = " x "), "\n")
initial_N <- model$N[1,,]  # First time step

for(i in seq_along(fish_grps)) {
  grp <- fish_grps[i]
  cat("\nFish", i, "(", Groups$Species[grp], ") initial abundances:\n")
  
  # Get size range for this group
  min_idx <- which.min(abs(log10(params$w) - Groups$W0[grp]))
  max_idx <- which.min(abs(log10(params$w) - Groups$Wmax[grp]))
  size_range <- min_idx:max_idx
  
  fish_abundance <- initial_N[grp, size_range]
  cat("  Size range indices:", min(size_range), "to", max(size_range), "\n")
  cat("  Total initial abundance:", sum(fish_abundance), "\n")
  cat("  Abundance range:", paste(round(range(fish_abundance), 6), collapse = " to "), "\n")
  
  # Check for mature individuals
  mat_idx <- which.min(abs(log10(params$w) - Groups$Wmat[grp]))
  if(mat_idx <= max_idx) {
    mature_abundance <- sum(initial_N[grp, mat_idx:max_idx])
    cat("  Mature individual abundance:", mature_abundance, "\n")
  } else {
    cat("  WARNING: No mature individuals in initial population!\n")
  }
}
cat("\n")

# Test 4: Step-by-step Reproduction Calculation
cat("=== TEST 4: STEP-BY-STEP REPRODUCTION CALCULATION ===\n")

# Simulate one time step of reproduction calculation manually
N_test <- initial_N
net_production <- matrix(runif(params$ngrps * params$ngrid, 0.01, 0.1), 
                        nrow = params$ngrps, ncol = params$ngrid)

repro_energy_test <- matrix(0, nrow = params$ngrps, ncol = params$ngrid)

for(i in seq_along(fish_grps)) {
  grp <- fish_grps[i]
  cat("\nFish", i, "(", Groups$Species[grp], ") reproduction calculation:\n")
  
  # Calculate maturity ogive
  maturity_ogive <- calculate_maturity_ogive(params$w, Groups$Wmat[grp])
  mature_count <- sum(maturity_ogive > 0.5)
  cat("  Mature size classes:", mature_count, "\n")
  
  # Calculate reproductive investment
  repro_invest <- calculate_reproductive_investment(params$w, Groups$ReproInvest[grp], 
                                                  Groups$ReproExp[grp], Groups$Wmat[grp])
  investing_count <- sum(repro_invest > 0)
  cat("  Size classes investing in reproduction:", investing_count, "\n")
  
  # Calculate reproductive energy
  repro_energy_test[grp,] <- net_production[grp,] * repro_invest * maturity_ogive
  total_repro <- sum(repro_energy_test[grp,])
  cat("  Total reproductive energy:", round(total_repro, 6), "\n")
  
  # Calculate SSB contribution
  SSB_contrib <- sum(repro_energy_test[grp,] * N_test[grp,] * params$w)
  cat("  SSB contribution:", round(SSB_contrib, 6), "\n")
}

# Test total SSB calculation
total_SSB <- calculate_spawning_stock_biomass(N_test, repro_energy_test, params$w, fish_grps)
cat("\nTotal SSB by fish group:", paste(round(total_SSB, 6), collapse = ", "), "\n")
cat("\n")

# Test 5: Full Model Run with Detailed Tracking
cat("=== TEST 5: FULL MODEL RUN WITH DETAILED TRACKING ===\n")

cat("Running ZooMSS with reproduction for", test_years, "years...\n")
tryCatch({
  results <- zoomss_run(model)
  
  cat("Model run completed successfully!\n")
  cat("Reproduction enabled:", ifelse(is.null(results$reproduction_enabled), "Unknown", results$reproduction_enabled), "\n")
  
  if(!is.null(results$SSB)) {
    cat("SSB array dimensions:", paste(dim(results$SSB), collapse = " x "), "\n")
    cat("Time points saved:", ncol(results$SSB), "\n")
    
    # Analyze SSB time series
    for(i in seq_along(fish_grps)) {
      grp_name <- Groups$Species[fish_grps[i]]
      SSB_timeseries <- results$SSB[i, ]
      cat("\n", grp_name, "SSB analysis:\n")
      cat("  SSB range:", paste(round(range(SSB_timeseries, na.rm = TRUE), 6), collapse = " to "), "\n")
      cat("  Mean SSB:", round(mean(SSB_timeseries, na.rm = TRUE), 6), "\n")
      cat("  Final SSB:", round(SSB_timeseries[length(SSB_timeseries)], 6), "\n")
      cat("  Non-zero SSB values:", sum(SSB_timeseries > 1e-10, na.rm = TRUE), "\n")
    }
  } else {
    cat("WARNING: SSB array not found in results!\n")
  }
  
  # Check final abundances
  final_N <- results$N[dim(results$N)[1],,]
  cat("\nFinal abundance analysis:\n")
  for(i in seq_along(fish_grps)) {
    grp <- fish_grps[i]
    grp_name <- Groups$Species[grp]
    
    min_idx <- which.min(abs(log10(params$w) - Groups$W0[grp]))
    max_idx <- which.min(abs(log10(params$w) - Groups$Wmax[grp]))
    mat_idx <- which.min(abs(log10(params$w) - Groups$Wmat[grp]))
    
    total_abundance <- sum(final_N[grp, min_idx:max_idx])
    if(mat_idx <= max_idx) {
      mature_abundance <- sum(final_N[grp, mat_idx:max_idx])
    } else {
      mature_abundance <- 0
    }
    
    cat("  ", grp_name, "- Total:", round(total_abundance, 6), 
        "Mature:", round(mature_abundance, 6), "\n")
  }
  
}, error = function(e) {
  cat("ERROR in model run:", e$message, "\n")
  cat("Stack trace:\n")
  traceback()
})

cat("\n=== DIAGNOSTIC COMPLETE ===\n")
cat("Check results above to identify why SSB = 0\n")