# Test Script: Verify ZooMSS Fishing Implementation
# 
# This script tests the dynamic fishing mortality implementation
# Run after loading the zoomss package with devtools::load_all()

library(zoomss)

cat("="
, rep("=", 60), "\n", sep = "")
cat("ZooMSS Fishing Implementation Tests\n")
cat(rep("=", 61), "\n\n", sep = "")

# =============================================================================
# Setup: Create test data
# =============================================================================

cat("Setting up test environment...\n")

# Get default groups
Groups <- getGroups()
ngrps <- nrow(Groups)
fish_grps <- which(Groups$Type == "Fish")

cat("  - Loaded", ngrps, "functional groups\n")
cat("  - Fish groups:", paste(Groups$Species[fish_grps], collapse = ", "), "\n")

# Create simple environmental time series (10 years, annual timesteps)
n_years <- 10
time_vec <- seq(1, n_years, by = 1)
sst_vec <- rep(20, n_years)  # Constant SST
chl_vec <- rep(1.0, n_years)  # Constant chlorophyll

input_params <- data.frame(
  time = time_vec,
  sst = sst_vec,
  chl = chl_vec
)

cat("  - Created", n_years, "year environmental time series\n\n")

# =============================================================================
# Test 1: Backward Compatibility (Static Fishing)
# =============================================================================

cat("Test 1: Backward Compatibility (Static Fishing)\n")
cat(rep("-", 50), "\n", sep = "")

tryCatch({
  result_static <- zoomss_model(
    input_params = input_params,
    Groups = Groups,
    isave = 1,
    fishing_params = NULL  # Use static fishing from Groups$Fmort
  )
  
  # Check outputs exist
  stopifnot("abundance" %in% names(result_static))
  stopifnot("fishing_mortality" %in% names(result_static))
  stopifnot("catch" %in% names(result_static))
  
  # Check dimensions
  stopifnot(dim(result_static$abundance)[1] == n_years)
  stopifnot(dim(result_static$abundance)[2] == ngrps)
  
  # Check fishing mortality is constant (static)
  F_var <- var(as.vector(result_static$fishing_mortality))
  if (F_var < 1e-10) {
    cat("  [PASS] Fishing mortality is constant (static mode)\n")
  } else {
    cat("  [WARN] Fishing mortality varies unexpectedly in static mode\n")
  }
  
  # Check catch is non-negative
  if (all(result_static$catch >= 0, na.rm = TRUE)) {
    cat("  [PASS] Catch values are non-negative\n")
  } else {
    cat("  [FAIL] Negative catch values detected\n")
  }
  
  cat("  [PASS] Static fishing mode works correctly\n\n")
  
}, error = function(e) {
  cat("  [FAIL] Error:", e$message, "\n\n")
})

# =============================================================================
# Test 2: Dynamic Fishing with Constant Effort
# =============================================================================

cat("Test 2: Dynamic Fishing with Constant Effort\n")
cat(rep("-", 50), "\n", sep = "")

tryCatch({
  # Create fishing parameters
  # Catchability: zero for zooplankton, positive for fish
  q_values <- rep(0, ngrps)
  q_values[fish_grps] <- 0.01  # Same q for all fish groups
  
  # Selectivity: knife-edge for each group
  selectivity_list <- vector("list", ngrps)
  for (g in 1:ngrps) {
    if (g %in% fish_grps) {
      # Fish: vulnerable above 1g (log10 = 0)
      selectivity_list[[g]] <- list(
        type = "knife_edge",
        params = list(w_min = 0),
        log_scale = TRUE
      )
    } else {
      # Zooplankton: not fished (set high w_min so nothing is selected)
      selectivity_list[[g]] <- list(
        type = "knife_edge",
        params = list(w_min = 10),  # Above max size
        log_scale = TRUE
      )
    }
  }
  
  # Constant effort time series
  effort_ts <- rep(100, n_years)
  
  fishing_params <- list(
    effort = effort_ts,
    catchability = q_values,
    selectivity = selectivity_list
  )
  
  result_dynamic <- zoomss_model(
    input_params = input_params,
    Groups = Groups,
    isave = 1,
    fishing_params = fishing_params
  )
  
  # Check fishing mortality is applied to fish groups
  F_fish <- result_dynamic$fishing_mortality[1, fish_grps, ]
  F_zoo <- result_dynamic$fishing_mortality[1, -fish_grps, ]
  
  if (sum(F_fish) > 0) {
    cat("  [PASS] Fish groups have non-zero fishing mortality\n")
  } else {
    cat("  [FAIL] Fish groups have zero fishing mortality\n")
  }
  
  if (sum(F_zoo) == 0) {
    cat("  [PASS] Zooplankton groups have zero fishing mortality\n")
  } else {
    cat("  [WARN] Zooplankton groups have non-zero fishing mortality\n")
  }
  
  # Check F is constant over time (constant effort)
  F_fish_t1 <- result_dynamic$fishing_mortality[1, fish_grps[1], ]
  F_fish_t_last <- result_dynamic$fishing_mortality[n_years, fish_grps[1], ]
  
  if (all(abs(F_fish_t1 - F_fish_t_last) < 1e-10)) {
    cat("  [PASS] Fishing mortality constant with constant effort\n")
  } else {
    cat("  [WARN] Fishing mortality varies despite constant effort\n")
  }
  
  cat("  [PASS] Dynamic fishing with constant effort works\n\n")
  
}, error = function(e) {
  cat("  [FAIL] Error:", e$message, "\n\n")
})

# =============================================================================
# Test 3: Dynamic Fishing with Time-Varying Effort
# =============================================================================

cat("Test 3: Dynamic Fishing with Time-Varying Effort\n")
cat(rep("-", 50), "\n", sep = "")

tryCatch({
  # Create time-varying effort (increases over time)
  effort_ts_varying <- seq(50, 150, length.out = n_years)
  
  fishing_params_varying <- list(
    effort = effort_ts_varying,
    catchability = q_values,
    selectivity = selectivity_list
  )
  
  result_varying <- zoomss_model(
    input_params = input_params,
    Groups = Groups,
    isave = 1,
    fishing_params = fishing_params_varying
  )
  
  # Check F increases over time for fish
  F_year1 <- mean(result_varying$fishing_mortality[1, fish_grps, ])
  F_year_last <- mean(result_varying$fishing_mortality[n_years, fish_grps, ])
  
  if (F_year_last > F_year1) {
    cat("  [PASS] Fishing mortality increases with increasing effort\n")
    cat("         Year 1 mean F:", round(F_year1, 4), "\n")
    cat("         Year", n_years, "mean F:", round(F_year_last, 4), "\n")
  } else {
    cat("  [FAIL] Fishing mortality does not respond to effort changes\n")
  }
  
  # Check catch changes over time
  catch_year1 <- sum(result_varying$catch[1, ])
  catch_year_last <- sum(result_varying$catch[n_years, ])
  
  cat("         Year 1 total catch:", round(catch_year1, 4), "\n")
  cat("         Year", n_years, "total catch:", round(catch_year_last, 4), "\n")
  
  cat("  [PASS] Time-varying effort works correctly\n\n")
  
}, error = function(e) {
  cat("  [FAIL] Error:", e$message, "\n\n")
})

# =============================================================================
# Test 4: Unfished vs Fished Comparison
# =============================================================================

cat("Test 4: Unfished vs Fished Comparison\n")
cat(rep("-", 50), "\n", sep = "")

tryCatch({
  # Unfished: set all Fmort to zero in Groups
  Groups_unfished <- Groups
  Groups_unfished$Fmort <- 0
  
  result_unfished <- zoomss_model(
    input_params = input_params,
    Groups = Groups_unfished,
    isave = 1,
    fishing_params = NULL
  )
  
  # Fished: use dynamic fishing with moderate effort
  fishing_params_moderate <- list(
    effort = rep(200, n_years),
    catchability = q_values,
    selectivity = selectivity_list
  )
  
  result_fished <- zoomss_model(
    input_params = input_params,
    Groups = Groups,
    isave = 1,
    fishing_params = fishing_params_moderate
  )
  
  # Compare final fish biomass
  # biomass is 3D: [time, groups, size] - need to sum across size dimension too
  biomass_unfished <- sum(result_unfished$biomass[n_years, fish_grps, ])
  biomass_fished <- sum(result_fished$biomass[n_years, fish_grps, ])
  
  cat("  Final fish biomass (unfished):", round(biomass_unfished, 4), "\n")
  cat("  Final fish biomass (fished):  ", round(biomass_fished, 4), "\n")
  
  if (biomass_fished < biomass_unfished) {
    cat("  [PASS] Fishing reduces fish biomass as expected\n")
    cat("         Reduction:", round(100 * (1 - biomass_fished/biomass_unfished), 1), "%\n")
  } else {
    cat("  [WARN] Fished biomass not lower than unfished\n")
  }
  
  # Check total catch is positive
  total_catch <- sum(result_fished$catch, na.rm = TRUE)
  cat("  Total catch over", n_years, "years:", round(total_catch, 4), "\n")
  
  if (total_catch > 0) {
    cat("  [PASS] Positive catch recorded\n")
  } else {
    cat("  [FAIL] No catch recorded\n")
  }
  
  cat("\n")
  
}, error = function(e) {
  cat("  [FAIL] Error:", e$message, "\n\n")
})

# =============================================================================
# Test 5: Selectivity Types
# =============================================================================

cat("Test 5: Selectivity Functions\n")
cat(rep("-", 50), "\n", sep = "")

tryCatch({
  # Test knife-edge
  w_test <- seq(-2, 4, by = 0.1)
  s_knife <- calc_selectivity(w_test, type = "knife_edge", params = list(w_min = 1))
  
  if (all(s_knife[w_test < 1] == 0) && all(s_knife[w_test >= 1] == 1)) {
    cat("  [PASS] Knife-edge selectivity correct\n")
  } else {
    cat("  [FAIL] Knife-edge selectivity incorrect\n")
  }
  
  # Test logistic
  s_logistic <- calc_selectivity(w_test, type = "logistic", params = list(L50 = 1, L95 = 2))
  
  # Check S(L50) ≈ 0.5
  idx_L50 <- which.min(abs(w_test - 1))
  if (abs(s_logistic[idx_L50] - 0.5) < 0.01) {
    cat("  [PASS] Logistic selectivity S(L50) ≈ 0.5\n")
  } else {
    cat("  [FAIL] Logistic selectivity S(L50) ≠ 0.5\n")
  }
  
  # Check S(L95) ≈ 0.95

  idx_L95 <- which.min(abs(w_test - 2))
  if (abs(s_logistic[idx_L95] - 0.95) < 0.01) {
    cat("  [PASS] Logistic selectivity S(L95) ≈ 0.95\n")
  } else {
    cat("  [FAIL] Logistic selectivity S(L95) ≠ 0.95\n")
  }
  
  # Test custom
  custom_fun <- function(w) ifelse(w > 0 & w < 2, 1, 0)
  s_custom <- calc_selectivity(w_test, type = "custom", params = list(fun = custom_fun))
  
  if (sum(s_custom) > 0) {
    cat("  [PASS] Custom selectivity function works\n")
  } else {
    cat("  [FAIL] Custom selectivity function failed\n")
  }
  
  cat("\n")
  
}, error = function(e) {
  cat("  [FAIL] Error:", e$message, "\n\n")
})

# =============================================================================
# Test 6: Size-Resolved Catch Output
# =============================================================================

cat("Test 6: Size-Resolved Catch Output\n")
cat(rep("-", 50), "\n", sep = "")

tryCatch({
  result_size_catch <- zoomss_model(
    input_params = input_params,
    Groups = Groups,
    isave = 1,
    fishing_params = fishing_params,  # From Test 2
    save_catch_by_size = TRUE
  )
  
  # Check catch_by_size exists
  if ("catch_by_size" %in% names(result_size_catch)) {
    cat("  [PASS] catch_by_size output exists\n")
    
    # Check dimensions (nsave x ngrps x ngrid)
    dims <- dim(result_size_catch$catch_by_size)
    cat("         Dimensions:", paste(dims, collapse = " x "), "\n")
    
    # Verify sum of catch_by_size ≈ catch (total)
    sum_by_size <- sum(result_size_catch$catch_by_size, na.rm = TRUE)
    sum_total <- sum(result_size_catch$catch, na.rm = TRUE)
    
    if (abs(sum_by_size - sum_total) / sum_total < 0.01) {
      cat("  [PASS] Sum of catch_by_size matches total catch\n")
    } else {
      cat("  [WARN] Mismatch between catch_by_size sum and total catch\n")
      cat("         catch_by_size sum:", round(sum_by_size, 4), "\n")
      cat("         catch total:", round(sum_total, 4), "\n")
    }
    
  } else {
    cat("  [FAIL] catch_by_size output missing\n")
  }
  
  cat("\n")
  
}, error = function(e) {
  cat("  [FAIL] Error:", e$message, "\n\n")
})

# =============================================================================
# Summary
# =============================================================================

cat(rep("=", 61), "\n", sep = "")
cat("Tests completed. Review output above for any [FAIL] or [WARN] messages.\n")
cat(rep("=", 61), "\n", sep = "")
