# =============================================================================
# 03_objective_function.R
# Calibration objective function for revised ZooMSS
# =============================================================================
#
# Purpose: Define the objective function that compares revised model outputs
#          against the original ZooMSS baseline across a chlorophyll gradient.
#          Designed to be sourced by 04_run_calibration.R.
#
# Usage:
#   source("calibration/03_objective_function.R")
#   obj_value <- calibration_objective(par, baseline, chl_levels, Groups)
#
# Prerequisites:
#   - Baseline from 01_generate_baseline.R
#   - future.apply::plan() already set up externally
# =============================================================================

#' Calibration Objective Function for Revised ZooMSS
#'
#' @param par Named numeric vector of parameters to optimise:
#'   - f_M: Metabolic fraction (single value applied to all groups)
#'   - K_growth_zoo_base: Base zooplankton growth fraction (scaled per group)
#'   - K_growth_fish: Fish growth fraction (single value for all fish)
#'   - repro_eff: Fish reproductive efficiency (egg-to-recruit survival)
#' @param baseline List from 01_generate_baseline.R
#' @param chl_levels Chlorophyll gradient to evaluate
#' @param Groups Default Groups data frame (modified by par)
#' @param verbose Logical, print progress
#' @return Scalar objective value (lower = better fit)
calibration_objective <- function(par, baseline, chl_levels, Groups, verbose = FALSE) {

  # ── Unpack and apply parameters ──
  Groups$f_M <- par["f_M"]

  zoo_idx  <- which(Groups$Type == "Zooplankton")
  fish_idx <- which(Groups$Type == "Fish")
  n_zoo  <- length(zoo_idx)
  n_fish <- length(fish_idx)

  # K_growth for zooplankton: group-specific scaling from a base value
  # Preserves relative rank order between groups
  K_base_zoo <- par["K_growth_zoo_base"]
  default_K_zoo <- Groups$K_growth[zoo_idx]
  K_relative <- default_K_zoo / mean(default_K_zoo)  # preserve relative pattern
  Groups$K_growth[zoo_idx] <- K_base_zoo * K_relative
  # Clamp to valid range
  Groups$K_growth[zoo_idx] <- pmin(pmax(Groups$K_growth[zoo_idx], 0.05), 0.49)

  # K_growth for fish
  Groups$K_growth[fish_idx] <- par["K_growth_fish"]

  # Reproductive efficiency for fish
  Groups$repro_eff[fish_idx] <- par["repro_eff"]

  # ── Validate energy budget closure ──
  R_frac <- 1 - Groups$f_M - Groups$K_growth
  if (any(R_frac < 0) || any(R_frac > 1)) {
    return(1e6)  # Penalty for invalid parameter combinations
  }

  # ── Run model across chlorophyll gradient (parallelised) ──
  n_chl <- length(chl_levels)

  run_single_chl <- function(chl, Groups) {
    tryCatch({
      env <- createInputParams(
        time = seq(0, 400, by = 0.1),
        sst  = 15,
        chl  = chl
      )

      mdl <- zoomss_model(
        input_params = env,
        Groups = Groups,
        isave = 10  # coarser saving for speed during calibration
      )

      # Extract steady-state biomass using getBiomass (final 100 years)
      Biomass <- getBiomass(mdl, units = "ww")
      time_vec <- mdl$time
      time_idx <- which(time_vec >= max(time_vec) - 100)
      avg_biomass <- apply(Biomass[time_idx, , , drop = FALSE], c(2, 3), mean)
      group_biomass <- rowSums(avg_biomass)

      list(group_biomass = group_biomass, success = TRUE)
    }, error = function(e) {
      list(group_biomass = NULL, success = FALSE, error = e$message)
    })
  }

  # Run in parallel using future.apply (plan set up externally)
  if (requireNamespace("future.apply", quietly = TRUE)) {
    chl_results <- future.apply::future_lapply(
      chl_levels, run_single_chl, Groups = Groups,
      future.seed = TRUE
    )
  } else {
    chl_results <- lapply(chl_levels, run_single_chl, Groups = Groups)
  }

  # ── Unpack results ──
  zoo_proportions <- matrix(NA, nrow = n_chl, ncol = n_zoo)
  fish_biomass    <- matrix(NA, nrow = n_chl, ncol = n_fish)
  total_biomass   <- numeric(n_chl)

  for (i in seq_along(chl_levels)) {
    res <- chl_results[[i]]
    if (res$success) {
      gb <- res$group_biomass

      zoo_bm <- gb[zoo_idx]
      zoo_proportions[i, ] <- zoo_bm / sum(zoo_bm)
      fish_biomass[i, ]    <- gb[fish_idx]
      total_biomass[i]     <- sum(gb)
    } else {
      if (verbose) cat("Error at chl =", chl_levels[i], ":", res$error, "\n")
    }
  }

  # ── Calculate objective ──
  # If any runs failed, return large penalty
  if (any(is.na(zoo_proportions))) {
    return(1e6 - sum(!is.na(zoo_proportions)))  # Partial credit for partial success
  }

  # Component 1: Zooplankton proportions across gradient (primary target)
  ss_zoo_prop <- sum((zoo_proportions - baseline$zoo_proportions)^2)

  # Component 2: Zoo:Fish biomass ratio across gradient (secondary target)
  zoo_total_bm <- rowSums(zoo_proportions * total_biomass)  # reconstruct zoo biomass
  fish_total_bm <- rowSums(fish_biomass)
  zoo_fish_ratio <- zoo_total_bm / pmax(fish_total_bm, 1e-20)
  # Use log-ratio to handle scale differences
  ss_ratio <- sum((log10(pmax(zoo_fish_ratio, 1e-20)) -
                    log10(pmax(baseline$zoo_fish_ratio, 1e-20)))^2, na.rm = TRUE)

  # Component 3: Total biomass magnitude (tertiary target)
  ss_total <- sum((log10(pmax(total_biomass, 1e-20)) -
                    log10(pmax(baseline$total_biomass, 1e-20)))^2, na.rm = TRUE)

  # Weighted objective
  w_prop  <- 1.0    # Primary: zooplankton community composition
  w_ratio <- 0.3    # Secondary: zoo:fish ratio
  w_total <- 0.1    # Tertiary: absolute biomass magnitude

  objective <- w_prop * ss_zoo_prop + w_ratio * ss_ratio + w_total * ss_total

  if (verbose) {
    cat(sprintf("f_M=%.3f K_zoo=%.3f K_fish=%.3f repro_eff=%.4f | obj=%.4f (prop=%.4f ratio=%.4f total=%.4f)\n",
                par["f_M"], par["K_growth_zoo_base"], par["K_growth_fish"],
                par["repro_eff"], objective, ss_zoo_prop, ss_ratio, ss_total))
  }

  return(objective)
}
