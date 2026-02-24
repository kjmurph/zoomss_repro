#!/usr/bin/env Rscript
# =============================================================================
# 06_fish_param_sensitivity.R
# Sensitivity Analysis: Fish Group Parameter Differentiation
# =============================================================================
#
# Tests the effect of differentiating fish functional group parameters
# in ZooMSS vs the current uniform parameterisation.
#
# Design:
#   Phase 0: Smoke test (short run to verify everything functions)
#   Phase 1: One-at-a-time (OAT) sensitivity — vary each parameter independently
#   Phase 2: Monte Carlo uncertainty analysis with Latin Hypercube Sampling
#
# Outputs measured:
#   - Biomass by fish group (equilibrium average)
#   - Fish biomass ratios (Fish_Small:Med:Large)
#   - Small:Large ratio (ecosystem indicator)
#   - Zoo:Fish ratio
#   - Size spectrum slope
#   - Steady-state diagnostics (CV, trend slope, group persistence)
#
# Prerequisites:
#   - Development branch of zoomss with ZSpre/ZSexp in GroupInputs
#   - Packages: future.apply, lhs, parallelly, ggplot2, patchwork, dplyr
#
# Output files:
#   - calibration/fish_sensitivity_smoke_test.rds
#   - calibration/fish_sensitivity_oat_results.rds
#   - calibration/fish_sensitivity_mc_results.rds
#   - calibration/fish_sensitivity_plots.pdf
#   - calibration/fish_sensitivity_full_results.rds
#
# Usage:
#   source("calibration/06_fish_param_sensitivity.R")
#   run_smoke_test()                              # Verify setup
#   results <- run_full_sensitivity(run_mc = FALSE)  # OAT only
#   results <- run_full_sensitivity(run_mc = TRUE)   # Full analysis
# =============================================================================


# =============================================================================
# 0. SETUP & CONFIGURATION
# =============================================================================

# --- Package loading ---
# Development branch: use devtools::load_all() for main session,
# then install so parallel workers (separate R processes) can also access it.
cat("Setting up environment...\n")
pkg_path <- normalizePath(".")

devtools::load_all(".")  # Load dev branch functions for interactive use

# Install dev version so parallel workers can load via library(zoomss)
devtools::install(pkg = pkg_path, quiet = TRUE, upgrade = "never",
                  dependencies = FALSE)
library(zoomss)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future.apply)
library(lhs)

# --- Parallel setup ---
n_workers <- min(parallelly::availableCores() - 2, 14)
plan(multisession, workers = n_workers)
cat("Using", n_workers, "of", parallelly::availableCores(), "available cores\n")

# Increase global size limit for future workers (model objects can be large)
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB

# --- Simulation settings ---
SIM_YEARS   <- 250          # Simulation length (years)
DT          <- 0.01         # Time step (years) — matches standard ZooMSS dt
ISAVE       <- 100          # Save every 100 steps (= 1 year with dt=0.01)
N_AVG_YEARS <- 50           # Final years for equilibrium metrics (years 200-250)
CHL_LEVELS  <- c(0.01, 0.1, 1.0)  # Chlorophyll gradient (mg m-3)
SST         <- 15           # Constant SST (degrees C)
N_MC        <- 100          # Monte Carlo replicates (LHS)
SEED        <- 42           # Random seed for reproducibility

# Smoke test settings (short run to verify)
SMOKE_YEARS <- 10
SMOKE_ISAVE <- 100          # Annual output (matching ISAVE / DT ratio)

# Minimum reproduction fraction (R_frac = 1 - f_M - K_growth)
MIN_R_FRAC <- 0.15

# Trend threshold for steady-state assessment (log10 biomass units per year)
# 0.005 corresponds to ~1.2% change per year
TREND_THRESHOLD <- 0.005

# Output directory
out_dir <- "calibration"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# =============================================================================
# 1. PARAMETER DEFINITIONS
# =============================================================================

# --- Baseline: current uniform fish parameters ---
# These match the existing GroupInputs defaults for fish groups
baseline_params <- list(
  PPMR       = c(Fish_Small = 100,   Fish_Med = 100,   Fish_Large = 100),
  FeedWidth  = c(Fish_Small = 1.3,   Fish_Med = 1.3,   Fish_Large = 1.3),
  f_M        = c(Fish_Small = 0.50,  Fish_Med = 0.50,  Fish_Large = 0.50),
  K_growth   = c(Fish_Small = 0.30,  Fish_Med = 0.30,  Fish_Large = 0.30),
  repro_eff  = c(Fish_Small = 0.001, Fish_Med = 0.001, Fish_Large = 0.001),
  ZSpre      = c(Fish_Small = 0.10,  Fish_Med = 0.10,  Fish_Large = 0.10),
  ZSexp      = c(Fish_Small = 0.30,  Fish_Med = 0.30,  Fish_Large = 0.30),
  Wmat       = c(Fish_Small = 0.0,   Fish_Med = 2.0,   Fish_Large = 4.0)
)

# --- Differentiated: literature-informed hypothesis ---
# Moderate differentiation for initial testing (can be widened later).
#
# Energy budget check (R_frac = 1 - f_M - K_growth >= 0.15):
#   Fish_Small: 1 - 0.50 - 0.30 = 0.20  >=  0.15  OK
#   Fish_Med:   1 - 0.50 - 0.25 = 0.25  >=  0.15  OK
#   Fish_Large: 1 - 0.45 - 0.22 = 0.33  >=  0.15  OK
#
# Wmat bounds (W0 < Wmat < Wmax):
#   Fish_Small: -3 < -0.5 < 2  OK
#   Fish_Med:   -3 <  1.5 < 4  OK
#   Fish_Large: -3 <  3.5 < 6  OK
#
differentiated_params <- list(
  PPMR       = c(Fish_Small = 100,    Fish_Med = 250,    Fish_Large = 500),
  FeedWidth  = c(Fish_Small = 1.1,    Fish_Med = 1.4,    Fish_Large = 1.2),
  f_M        = c(Fish_Small = 0.50,   Fish_Med = 0.50,   Fish_Large = 0.45),
  K_growth   = c(Fish_Small = 0.30,   Fish_Med = 0.25,   Fish_Large = 0.22),
  repro_eff  = c(Fish_Small = 0.002,  Fish_Med = 0.001,  Fish_Large = 0.0005),
  ZSpre      = c(Fish_Small = 0.12,   Fish_Med = 0.10,   Fish_Large = 0.15),
  ZSexp      = c(Fish_Small = 0.30,   Fish_Med = 0.30,   Fish_Large = 0.40),
  Wmat       = c(Fish_Small = -0.5,   Fish_Med = 1.5,    Fish_Large = 3.5)
)

# --- Monte Carlo uncertainty distributions ---
# Tighter ranges for initial testing; can be widened if needed.
# R_frac >= 0.15 constraint enforced post-hoc in draw_lhs_params().
param_ranges <- list(
  PPMR = list(
    Fish_Small = list(dist = "lognormal", meanlog = log(100),  sdlog = 0.20),
    Fish_Med   = list(dist = "lognormal", meanlog = log(250),  sdlog = 0.25),
    Fish_Large = list(dist = "lognormal", meanlog = log(500),  sdlog = 0.25)
  ),
  FeedWidth = list(
    Fish_Small = list(dist = "normal", mean = 1.1, sd = 0.10),
    Fish_Med   = list(dist = "normal", mean = 1.4, sd = 0.10),
    Fish_Large = list(dist = "normal", mean = 1.2, sd = 0.10)
  ),
  f_M = list(
    Fish_Small = list(dist = "normal", mean = 0.50, sd = 0.02),
    Fish_Med   = list(dist = "normal", mean = 0.50, sd = 0.02),
    Fish_Large = list(dist = "normal", mean = 0.45, sd = 0.02)
  ),
  K_growth = list(
    Fish_Small = list(dist = "normal", mean = 0.30, sd = 0.02),
    Fish_Med   = list(dist = "normal", mean = 0.25, sd = 0.02),
    Fish_Large = list(dist = "normal", mean = 0.22, sd = 0.02)
  ),
  repro_eff = list(
    Fish_Small = list(dist = "lognormal", meanlog = log(0.002),  sdlog = 0.30),
    Fish_Med   = list(dist = "lognormal", meanlog = log(0.001),  sdlog = 0.30),
    Fish_Large = list(dist = "lognormal", meanlog = log(0.0005), sdlog = 0.30)
  ),
  ZSpre = list(
    Fish_Small = list(dist = "lognormal", meanlog = log(0.12), sdlog = 0.20),
    Fish_Med   = list(dist = "lognormal", meanlog = log(0.10), sdlog = 0.20),
    Fish_Large = list(dist = "lognormal", meanlog = log(0.15), sdlog = 0.20)
  ),
  ZSexp = list(
    Fish_Small = list(dist = "normal", mean = 0.30, sd = 0.03),
    Fish_Med   = list(dist = "normal", mean = 0.30, sd = 0.03),
    Fish_Large = list(dist = "normal", mean = 0.40, sd = 0.05)
  ),
  Wmat = list(
    Fish_Small = list(dist = "normal", mean = -0.5, sd = 0.30),
    Fish_Med   = list(dist = "normal", mean =  1.5, sd = 0.30),
    Fish_Large = list(dist = "normal", mean =  3.5, sd = 0.30)
  )
)

# Parameter bounds (for clipping extreme draws)
param_bounds <- list(
  PPMR      = c(lower = 30,   upper = 2000),
  FeedWidth = c(lower = 0.5,  upper = 2.5),
  f_M       = c(lower = 0.30, upper = 0.60),
  K_growth  = c(lower = 0.10, upper = 0.40),
  repro_eff = c(lower = 1e-5, upper = 0.01),
  ZSpre     = c(lower = 0.01, upper = 0.50),
  ZSexp     = c(lower = 0.10, upper = 0.80)
)

# Group-specific Wmat bounds (W0 + buffer < Wmat < Wmax - buffer)
wmat_bounds <- list(
  Fish_Small = c(lower = -2.5, upper = 1.5),   # W0 = -3, Wmax = 2
  Fish_Med   = c(lower = -2.0, upper = 3.5),   # W0 = -3, Wmax = 4
  Fish_Large = c(lower = -1.5, upper = 5.5)    # W0 = -3, Wmax = 6
)


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

#' Apply fish parameter set to a Groups data frame
#'
#' Modifies only fish group rows, preserving all zooplankton parameters.
#' Validates energy budget closure (R_frac >= MIN_R_FRAC) and Wmat bounds.
#'
#' @param Groups Default Groups data frame from getGroups()
#' @param params Named list of parameter vectors, each with Fish_Small/Med/Large elements
#' @return Modified Groups data frame
apply_fish_params <- function(Groups, params) {

  fish_idx <- which(Groups$Type == "Fish")
  fish_names <- Groups$Species[fish_idx]

  for (param_name in names(params)) {
    vals <- params[[param_name]]
    for (i in seq_along(fish_idx)) {
      fg <- fish_idx[i]
      fn <- fish_names[i]
      if (fn %in% names(vals)) {
        Groups[[param_name]][fg] <- vals[[fn]]
      }
    }
  }

  # Validate energy budget: R_frac = 1 - f_M - K_growth >= MIN_R_FRAC
  for (fg in fish_idx) {
    R_frac <- 1 - Groups$f_M[fg] - Groups$K_growth[fg]
    if (R_frac < MIN_R_FRAC) {
      warning("Energy budget constraint for ", Groups$Species[fg],
              ": R_frac = ", round(R_frac, 3), " < ", MIN_R_FRAC,
              ". Reducing K_growth to enforce constraint.")
      Groups$K_growth[fg] <- 1 - Groups$f_M[fg] - MIN_R_FRAC
    }
  }

  # Validate Wmat bounds: W0 < Wmat < Wmax
  if ("Wmat" %in% names(params)) {
    for (fg in fish_idx) {
      if (Groups$Wmat[fg] <= Groups$W0[fg]) {
        Groups$Wmat[fg] <- Groups$W0[fg] + 0.5
        warning("Wmat clipped to W0 + 0.5 for ", Groups$Species[fg])
      }
      if (Groups$Wmat[fg] >= Groups$Wmax[fg]) {
        Groups$Wmat[fg] <- Groups$Wmax[fg] - 0.5
        warning("Wmat clipped to Wmax - 0.5 for ", Groups$Species[fg])
      }
    }
  }

  return(Groups)
}


#' Run a single ZooMSS simulation and extract metrics
#'
#' @param Groups Groups data frame (already modified with test parameters)
#' @param chl Chlorophyll level (mg m-3)
#' @param sim_years Simulation duration (years)
#' @param isave Save frequency (time steps)
#' @param n_avg_years Years at end of simulation for averaging
#' @return Named list with $metrics, $diagnostics, $success, and optionally $error
run_single_sim <- function(Groups, chl, sim_years = SIM_YEARS, isave = ISAVE,
                           n_avg_years = N_AVG_YEARS) {

  tryCatch({
    env <- createInputParams(
      time = seq(0, sim_years, by = DT),
      sst  = SST,
      chl  = chl
    )

    mdl <- suppressMessages(
      zoomss_model(input_params = env, Groups = Groups, isave = isave)
    )

    metrics     <- extract_metrics(mdl, n_avg_years = n_avg_years)
    diagnostics <- assess_steady_state(mdl, n_years = n_avg_years)

    list(metrics = metrics, diagnostics = diagnostics, success = TRUE)

  }, error = function(e) {
    list(metrics = NULL, diagnostics = NULL, success = FALSE, error = e$message)
  })
}


#' Extract equilibrium metrics from model output
#'
#' Computes biomass, ratios, steady-state indicators, and size spectrum
#' statistics from the final n_avg_years of a simulation.
#'
#' @param mdl ZooMSS model output from zoomss_model()
#' @param n_avg_years Years to average over (from end of simulation)
#' @return Named numeric vector of metrics
extract_metrics <- function(mdl, n_avg_years = 50) {

  fish_grps  <- mdl$param$fish_grps
  fish_names <- mdl$param$Groups$Species[fish_grps]
  zoo_grps   <- mdl$param$zoo_grps

  time_vec  <- mdl$time
  avg_start <- max(time_vec) - n_avg_years
  time_idx  <- which(time_vec >= avg_start)

  # --- Group biomass (wet weight, time-averaged over final years) ---
  biomass_3d    <- mdl$biomass   # nsave x ngrps x ngrid
  group_bm_ts   <- apply(biomass_3d[time_idx, , , drop = FALSE], c(1, 2), sum)
  avg_group_bm  <- colMeans(group_bm_ts)

  fish_bm    <- avg_group_bm[fish_grps]
  names(fish_bm) <- fish_names
  total_fish <- sum(fish_bm)
  total_zoo  <- sum(avg_group_bm[zoo_grps])

  # Fish biomass ratios
  fish_ratios <- fish_bm / max(total_fish, 1e-30)
  names(fish_ratios) <- paste0("ratio_", fish_names)

  # Small:Large ratio (key ecosystem indicator)
  SL_ratio <- fish_bm["Fish_Small"] / max(fish_bm["Fish_Large"], 1e-30)

  # Zoo:Fish ratio
  ZF_ratio <- total_zoo / max(total_fish, 1e-30)

  # --- Steady-state diagnostics per fish group ---
  # Coefficient of variation over final years
  cv_fish <- sapply(fish_grps, function(j) {
    bm_ts <- group_bm_ts[, j]
    if (mean(bm_ts) > 1e-30) sd(bm_ts) / mean(bm_ts) else NA
  })
  names(cv_fish) <- paste0("CV_", fish_names)

  # Trend: slope of log10(biomass) vs time over final years
  final_times <- time_vec[time_idx]
  trend_fish <- sapply(fish_grps, function(j) {
    bm_ts <- group_bm_ts[, j]
    if (all(bm_ts > 0)) {
      coef(lm(log10(bm_ts) ~ final_times))[2]
    } else {
      NA
    }
  })
  names(trend_fish) <- paste0("trend_", fish_names)

  # Group persistence (biomass above a minimal threshold)
  persists <- as.numeric(fish_bm > 1e-10)
  names(persists) <- paste0("persists_", fish_names)

  # --- Size spectrum slope (fish component) ---
  avg_abundance <- apply(mdl$abundance[time_idx, , , drop = FALSE], c(2, 3), mean)
  fish_abundance <- colSums(avg_abundance[fish_grps, , drop = FALSE])
  w_log10 <- mdl$param$w_log10

  fish_present <- fish_abundance > 0
  if (sum(fish_present) > 5) {
    lm_fit <- lm(log10(fish_abundance[fish_present]) ~ w_log10[fish_present])
    spectrum_slope     <- coef(lm_fit)[2]
    spectrum_intercept <- coef(lm_fit)[1]
  } else {
    spectrum_slope     <- NA
    spectrum_intercept <- NA
  }

  # --- Combine all metrics ---
  metrics <- c(
    fish_bm,
    fish_ratios,
    SL_ratio           = unname(SL_ratio),
    ZF_ratio           = unname(ZF_ratio),
    total_fish         = total_fish,
    total_zoo          = total_zoo,
    cv_fish,
    trend_fish,
    persists,
    spectrum_slope     = unname(spectrum_slope),
    spectrum_intercept = unname(spectrum_intercept)
  )

  return(metrics)
}


#' Assess steady-state for all groups
#'
#' Tests whether each functional group has achieved an oscillating steady
#' state over the final n_years: low trend in log10(biomass) and persistence.
#'
#' @param mdl ZooMSS model output
#' @param n_years Years to assess (from end of simulation)
#' @return Data frame with per-group diagnostics
assess_steady_state <- function(mdl, n_years = 50) {

  time_vec    <- mdl$time
  time_idx    <- which(time_vec >= max(time_vec) - n_years)
  final_times <- time_vec[time_idx]

  biomass_3d   <- mdl$biomass
  group_bm_ts  <- apply(biomass_3d[time_idx, , , drop = FALSE], c(1, 2), sum)
  n_grps       <- ncol(group_bm_ts)

  diagnostics <- data.frame(
    group        = mdl$param$Groups$Species,
    type         = mdl$param$Groups$Type,
    mean_biomass = colMeans(group_bm_ts),
    sd_biomass   = apply(group_bm_ts, 2, sd),
    stringsAsFactors = FALSE
  )

  diagnostics$cv <- diagnostics$sd_biomass / pmax(diagnostics$mean_biomass, 1e-30)

  # Linear trend in log10(biomass) over final years
  diagnostics$trend_slope <- sapply(1:n_grps, function(j) {
    bm <- group_bm_ts[, j]
    if (all(bm > 0)) {
      coef(lm(log10(bm) ~ final_times))[2]
    } else {
      NA
    }
  })

  diagnostics$persists <- diagnostics$mean_biomass > 1e-10

  # Steady state: persists AND |trend| < threshold
  diagnostics$is_steady <- diagnostics$persists &
    !is.na(diagnostics$trend_slope) &
    abs(diagnostics$trend_slope) < TREND_THRESHOLD

  return(diagnostics)
}


#' Generate parameter sets using Latin Hypercube Sampling
#'
#' Uses lhs::randomLHS for space-filling design, then maps [0,1] samples
#' to the specified parameter distributions via inverse CDF. Enforces
#' R_frac >= MIN_R_FRAC and Wmat bounds post-hoc.
#'
#' @param n Number of parameter sets to generate
#' @param ranges Distribution specifications (param_ranges format)
#' @param seed Random seed for reproducibility
#' @return List of n parameter sets (each structured like differentiated_params)
draw_lhs_params <- function(n, ranges, seed = SEED) {

  set.seed(seed)

  fish_groups <- c("Fish_Small", "Fish_Med", "Fish_Large")
  param_names <- names(ranges)

  # Dimension mapping: each column of LHS = one (param, group) combination
  dim_map <- expand.grid(
    group = fish_groups,
    param = param_names,
    stringsAsFactors = FALSE
  )
  n_dim <- nrow(dim_map)

  # Generate LHS design (n x n_dim matrix of [0, 1] values)
  lhs_design <- lhs::randomLHS(n, n_dim)

  param_sets <- vector("list", n)

  for (i in 1:n) {
    drawn <- list()

    for (pn in param_names) {
      vals <- numeric(3)
      names(vals) <- fish_groups

      for (gi in seq_along(fish_groups)) {
        gn <- fish_groups[gi]

        # Find the LHS column for this (param, group)
        d_idx <- which(dim_map$param == pn & dim_map$group == gn)
        u <- lhs_design[i, d_idx]

        spec <- ranges[[pn]][[gn]]

        # Inverse CDF: map [0,1] to parameter distribution
        val <- switch(spec$dist,
          "normal"    = qnorm(u, mean = spec$mean, sd = spec$sd),
          "lognormal" = qlnorm(u, meanlog = spec$meanlog, sdlog = spec$sdlog),
          stop("Unknown distribution: ", spec$dist)
        )

        # Apply bounds
        if (pn == "Wmat") {
          bounds <- wmat_bounds[[gn]]
          val <- max(bounds["lower"], min(val, bounds["upper"]))
        } else {
          bounds <- param_bounds[[pn]]
          val <- max(bounds["lower"], min(val, bounds["upper"]))
        }

        vals[gn] <- val
      }
      drawn[[pn]] <- vals
    }

    # Enforce energy budget: f_M + K_growth <= 1 - MIN_R_FRAC
    max_budget <- 1.0 - MIN_R_FRAC
    for (gn in fish_groups) {
      budget <- drawn$f_M[gn] + drawn$K_growth[gn]
      if (budget > max_budget) {
        # Scale both down proportionally to preserve relative magnitudes
        scale_factor <- max_budget / budget
        drawn$f_M[gn]     <- drawn$f_M[gn] * scale_factor
        drawn$K_growth[gn] <- drawn$K_growth[gn] * scale_factor
      }
    }

    param_sets[[i]] <- drawn
  }

  return(param_sets)
}


# =============================================================================
# 3. PHASE 0: SMOKE TEST
# =============================================================================

#' Run minimal tests to verify the full pipeline functions correctly
#'
#' Tests: baseline run, differentiated run, metric extraction,
#' parallel execution, and LHS parameter generation.
#'
#' @return TRUE (invisible) if all tests pass; stops on failure
run_smoke_test <- function() {

  cat("\n",
      strrep("=", 56), "\n",
      "  Phase 0: Smoke Test\n",
      strrep("=", 56), "\n\n", sep = "")

  Groups_default <- getGroups()

  # Test 1: Baseline configuration, short run
  cat("Test 1: Baseline config (", SMOKE_YEARS, "yr, chl = 0.1)...\n")
  Groups_baseline <- apply_fish_params(Groups_default, baseline_params)
  result_baseline <- run_single_sim(
    Groups_baseline, chl = 0.1,
    sim_years = SMOKE_YEARS, isave = SMOKE_ISAVE, n_avg_years = 5
  )

  if (!result_baseline$success) {
    stop("SMOKE TEST FAILED (baseline): ", result_baseline$error)
  }
  cat("  OK — Baseline ran successfully\n")
  cat("  Fish biomass: ",
      paste(names(result_baseline$metrics[c("Fish_Small", "Fish_Med", "Fish_Large")]),
            "=",
            formatC(result_baseline$metrics[c("Fish_Small", "Fish_Med", "Fish_Large")],
                    format = "e", digits = 3),
            collapse = ", "), "\n")

  # Test 2: Differentiated configuration
  cat("Test 2: Differentiated config (", SMOKE_YEARS, "yr, chl = 0.1)...\n")
  Groups_diff <- apply_fish_params(Groups_default, differentiated_params)
  result_diff <- run_single_sim(
    Groups_diff, chl = 0.1,
    sim_years = SMOKE_YEARS, isave = SMOKE_ISAVE, n_avg_years = 5
  )

  if (!result_diff$success) {
    stop("SMOKE TEST FAILED (differentiated): ", result_diff$error)
  }
  cat("  OK — Differentiated ran successfully\n")
  cat("  Fish biomass: ",
      paste(names(result_diff$metrics[c("Fish_Small", "Fish_Med", "Fish_Large")]),
            "=",
            formatC(result_diff$metrics[c("Fish_Small", "Fish_Med", "Fish_Large")],
                    format = "e", digits = 3),
            collapse = ", "), "\n")

  # Test 3: Verify metric fields
  cat("Test 3: Metric field verification...\n")
  expected_fields <- c("Fish_Small", "Fish_Med", "Fish_Large",
                        "ratio_Fish_Small", "ratio_Fish_Med", "ratio_Fish_Large",
                        "SL_ratio", "ZF_ratio",
                        "total_fish", "total_zoo",
                        "CV_Fish_Small", "CV_Fish_Med", "CV_Fish_Large",
                        "trend_Fish_Small", "trend_Fish_Med", "trend_Fish_Large",
                        "persists_Fish_Small", "persists_Fish_Med", "persists_Fish_Large",
                        "spectrum_slope", "spectrum_intercept")
  missing_fields <- setdiff(expected_fields, names(result_baseline$metrics))
  if (length(missing_fields) > 0) {
    stop("SMOKE TEST FAILED: missing metric fields: ",
         paste(missing_fields, collapse = ", "))
  }
  cat("  OK — All", length(expected_fields), "expected metric fields present\n")

  # Test 4: Diagnostics structure
  cat("Test 4: Steady-state diagnostics structure...\n")
  diag <- result_baseline$diagnostics
  diag_cols <- c("group", "type", "mean_biomass", "sd_biomass",
                  "cv", "trend_slope", "persists", "is_steady")
  missing_diag <- setdiff(diag_cols, names(diag))
  if (length(missing_diag) > 0) {
    stop("SMOKE TEST FAILED: missing diagnostic columns: ",
         paste(missing_diag, collapse = ", "))
  }
  cat("  OK — Diagnostics structure valid (", nrow(diag), " groups)\n")

  # Test 5: Parallel execution
  cat("Test 5: Parallel execution (2 runs on", n_workers, "workers)...\n")
  par_results <- future_lapply(c(0.01, 1.0), function(chl_val) {
    run_single_sim(
      Groups_baseline, chl = chl_val,
      sim_years = SMOKE_YEARS, isave = SMOKE_ISAVE, n_avg_years = 5
    )
  }, future.seed = TRUE, future.packages = "zoomss")

  n_success <- sum(sapply(par_results, function(r) r$success))
  if (n_success != 2) {
    failed_msgs <- sapply(par_results[!sapply(par_results, function(r) r$success)],
                          function(r) r$error)
    stop("SMOKE TEST FAILED: parallel execution (", n_success, "/2 succeeded). ",
         "Errors: ", paste(failed_msgs, collapse = "; "))
  }
  cat("  OK — Parallel execution successful\n")

  # Test 6: LHS parameter generation
  cat("Test 6: LHS parameter generation (5 samples, 24 dimensions)...\n")
  test_params <- draw_lhs_params(5, param_ranges, seed = SEED)

  for (i in seq_along(test_params)) {
    for (gn in c("Fish_Small", "Fish_Med", "Fish_Large")) {
      R_frac <- 1 - test_params[[i]]$f_M[gn] - test_params[[i]]$K_growth[gn]
      if (R_frac < MIN_R_FRAC) {
        stop("SMOKE TEST FAILED: R_frac constraint violated for ", gn,
             " in LHS sample ", i, " (R_frac = ", round(R_frac, 3), ")")
      }

      wmat_val <- test_params[[i]]$Wmat[gn]
      wb <- wmat_bounds[[gn]]
      if (wmat_val < wb["lower"] || wmat_val > wb["upper"]) {
        stop("SMOKE TEST FAILED: Wmat out of bounds for ", gn,
             " in LHS sample ", i, " (Wmat = ", round(wmat_val, 2), ")")
      }
    }
  }
  cat("  OK — LHS generation valid (all R_frac >=", MIN_R_FRAC,
      ", all Wmat within bounds)\n")

  # Save smoke test results
  smoke_output <- list(
    baseline_metrics = result_baseline$metrics,
    diff_metrics     = result_diff$metrics,
    baseline_diag    = result_baseline$diagnostics,
    diff_diag        = result_diff$diagnostics,
    parallel_ok      = TRUE,
    lhs_ok           = TRUE,
    timestamp        = Sys.time()
  )
  saveRDS(smoke_output, file.path(out_dir, "fish_sensitivity_smoke_test.rds"))

  cat("\n", strrep("=", 56), "\n",
      "  All smoke tests passed\n",
      strrep("=", 56), "\n\n", sep = "")

  return(invisible(TRUE))
}


# =============================================================================
# 4. PHASE 1: OAT SENSITIVITY
# =============================================================================

#' Run One-at-a-Time sensitivity analysis
#'
#' For each parameter, switches it from uniform to differentiated values
#' while all other parameters remain at baseline. Runs across 3 CHL levels.
#' Also runs baseline (all uniform) and fully differentiated configurations.
#'
#' @return List with tasks, config_names, configs, results, elapsed_minutes
run_oat_sensitivity <- function() {

  cat("\n",
      strrep("=", 56), "\n",
      "  Phase 1: One-at-a-Time Sensitivity Analysis\n",
      strrep("=", 56), "\n\n", sep = "")

  Groups_default <- getGroups()

  # Build OAT configurations
  params_to_test <- c("PPMR", "FeedWidth", "f_M", "K_growth",
                       "repro_eff", "ZSpre", "ZSexp", "Wmat")

  configs      <- list()
  config_names <- character()

  # Config 1: Baseline (all uniform)
  configs[[1]]      <- baseline_params
  config_names[1]   <- "baseline"

  # Configs 2-9: Each parameter differentiated individually
  for (i in seq_along(params_to_test)) {
    pn <- params_to_test[i]
    test_params <- baseline_params
    test_params[[pn]] <- differentiated_params[[pn]]
    configs[[i + 1]]      <- test_params
    config_names[i + 1]   <- pn
  }

  # Config 10: Fully differentiated
  n_last <- length(configs) + 1
  configs[[n_last]]      <- differentiated_params
  config_names[n_last]   <- "all_differentiated"

  n_configs <- length(configs)

  # Build flat task list for parallelisation
  tasks <- expand.grid(
    config_id = 1:n_configs,
    chl       = CHL_LEVELS,
    stringsAsFactors = FALSE
  )
  tasks$config_name <- config_names[tasks$config_id]

  cat("Configurations: ", n_configs, "\n")
  cat("CHL levels:     ", paste(CHL_LEVELS, collapse = ", "), " mg/m3\n")
  cat("Total runs:     ", nrow(tasks), "\n")
  cat("Workers:        ", n_workers, "\n")
  cat("Sim years:      ", SIM_YEARS, "\n\n")

  # Prepare Groups data frames for each configuration
  groups_list <- lapply(configs, function(params) {
    apply_fish_params(Groups_default, params)
  })

  cat("Launching", nrow(tasks), "simulations...\n")
  t_start <- Sys.time()

  results <- future_lapply(1:nrow(tasks), function(task_idx) {
    cfg_id  <- tasks$config_id[task_idx]
    chl_val <- tasks$chl[task_idx]

    run_single_sim(
      Groups     = groups_list[[cfg_id]],
      chl        = chl_val,
      sim_years  = SIM_YEARS,
      isave      = ISAVE,
      n_avg_years = N_AVG_YEARS
    )
  }, future.seed = TRUE, future.packages = "zoomss")

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")

  n_success <- sum(sapply(results, function(r) r$success))
  n_failed  <- sum(!sapply(results, function(r) r$success))
  cat("\nOAT phase complete in", round(as.numeric(t_elapsed), 1), "minutes\n")
  cat("Successful:", n_success, "/ Failed:", n_failed, "\n")

  if (n_failed > 0) {
    failed_idx <- which(!sapply(results, function(r) r$success))
    for (fi in failed_idx) {
      cat("  FAILED: config=", tasks$config_name[fi],
          ", chl=", tasks$chl[fi],
          ", error=", results[[fi]]$error, "\n")
    }
  }

  # Organise output
  oat_output <- list(
    tasks           = tasks,
    config_names    = config_names,
    configs         = configs,
    results         = results,
    elapsed_minutes = as.numeric(t_elapsed)
  )

  # Checkpoint
  saveRDS(oat_output, file.path(out_dir, "fish_sensitivity_oat_results.rds"))
  cat("Saved OAT results to", file.path(out_dir, "fish_sensitivity_oat_results.rds"), "\n")

  return(oat_output)
}


#' Calculate sensitivity indices from OAT results
#'
#' For each (parameter, response, chl) combination, computes the absolute
#' and relative change from baseline.
#'
#' @param oat_output Output from run_oat_sensitivity()
#' @return Data frame with columns: parameter, response, chl, baseline_value,
#'         test_value, abs_change, rel_change
calculate_oat_sensitivity <- function(oat_output) {

  tasks        <- oat_output$tasks
  results      <- oat_output$results
  config_names <- oat_output$config_names

  response_vars <- c("Fish_Small", "Fish_Med", "Fish_Large",
                      "SL_ratio", "ZF_ratio", "total_fish", "total_zoo",
                      "spectrum_slope")

  sensitivity_rows <- list()

  for (chl_val in CHL_LEVELS) {
    # Baseline for this chl
    base_idx <- which(tasks$config_name == "baseline" & tasks$chl == chl_val)
    if (length(base_idx) == 0 || !results[[base_idx]]$success) next
    base_metrics <- results[[base_idx]]$metrics

    # Compare each OAT configuration to baseline
    param_configs <- setdiff(config_names, "baseline")

    for (cfg_name in param_configs) {
      test_idx <- which(tasks$config_name == cfg_name & tasks$chl == chl_val)
      if (length(test_idx) == 0 || !results[[test_idx]]$success) next
      test_metrics <- results[[test_idx]]$metrics

      for (rv in response_vars) {
        if (!(rv %in% names(base_metrics) && rv %in% names(test_metrics))) next

        base_val <- base_metrics[[rv]]
        test_val <- test_metrics[[rv]]

        abs_change <- test_val - base_val
        rel_change <- if (abs(base_val) > 1e-30) {
          100 * (test_val - base_val) / base_val
        } else {
          NA
        }

        sensitivity_rows <- c(sensitivity_rows, list(data.frame(
          parameter      = cfg_name,
          response       = rv,
          chl            = chl_val,
          baseline_value = base_val,
          test_value     = test_val,
          abs_change     = abs_change,
          rel_change     = rel_change,
          stringsAsFactors = FALSE
        )))
      }
    }
  }

  sensitivity_df <- dplyr::bind_rows(sensitivity_rows)
  return(sensitivity_df)
}


# =============================================================================
# 5. PHASE 2: MONTE CARLO SENSITIVITY (LHS)
# =============================================================================

#' Run Monte Carlo sensitivity analysis with Latin Hypercube Sampling
#'
#' Generates n_mc parameter sets using LHS for space-filling design, then
#' runs each across the CHL gradient. Includes checkpointing for robustness.
#'
#' @param n_mc Number of LHS replicates
#' @return List with mc_df, param_sets, n_mc, chl_levels, elapsed_minutes
run_mc_sensitivity <- function(n_mc = N_MC) {

  cat("\n",
      strrep("=", 56), "\n",
      "  Phase 2: Monte Carlo Sensitivity (LHS)\n",
      strrep("=", 56), "\n\n", sep = "")

  Groups_default <- getGroups()

  # Generate LHS parameter sets
  n_dim <- length(names(param_ranges)) * 3
  cat("LHS design: ", n_mc, " samples x ", n_dim, " dimensions\n")
  param_sets <- draw_lhs_params(n_mc, param_ranges, seed = SEED)

  # Build flat task list
  tasks <- expand.grid(
    replicate = 1:n_mc,
    chl       = CHL_LEVELS,
    stringsAsFactors = FALSE
  )
  n_tasks <- nrow(tasks)

  cat("Total MC runs:  ", n_tasks, "\n")
  cat("Workers:        ", n_workers, "\n")
  cat("Est. batches:   ", ceiling(n_tasks / n_workers), "\n\n")

  # Prepare Groups for each replicate
  groups_list <- lapply(param_sets, function(params) {
    apply_fish_params(Groups_default, params)
  })

  # --- Run with checkpointing ---
  batch_size   <- n_workers * 3   # ~3 rounds per checkpoint
  n_batches    <- ceiling(n_tasks / batch_size)
  all_results  <- vector("list", n_tasks)

  # Check for existing checkpoint (resume capability)
  mc_checkpoint <- file.path(out_dir, "fish_sensitivity_mc_checkpoint.rds")
  start_batch <- 1

  if (file.exists(mc_checkpoint)) {
    checkpoint  <- readRDS(mc_checkpoint)
    all_results <- checkpoint$results
    start_batch <- checkpoint$next_batch
    cat("Resuming from batch", start_batch, "of", n_batches, "\n")
  }

  cat("Running MC simulations...\n")
  t_start <- Sys.time()

  for (b in start_batch:n_batches) {
    batch_start <- (b - 1) * batch_size + 1
    batch_end   <- min(b * batch_size, n_tasks)
    batch_idx   <- batch_start:batch_end

    cat(sprintf("  Batch %d/%d (tasks %d-%d)...",
                b, n_batches, batch_start, batch_end))

    batch_results <- future_lapply(batch_idx, function(task_idx) {
      rep_id  <- tasks$replicate[task_idx]
      chl_val <- tasks$chl[task_idx]

      result <- run_single_sim(
        Groups      = groups_list[[rep_id]],
        chl         = chl_val,
        sim_years   = SIM_YEARS,
        isave       = ISAVE,
        n_avg_years = N_AVG_YEARS
      )

      # Attach parameter values for PRCC analysis
      if (result$success) {
        result$param_values <- unlist(param_sets[[rep_id]])
      }
      result$replicate <- rep_id
      result$chl       <- chl_val

      return(result)
    }, future.seed = TRUE, future.packages = "zoomss")

    all_results[batch_idx] <- batch_results

    # Save checkpoint
    saveRDS(list(results = all_results, next_batch = b + 1), mc_checkpoint)

    n_success <- sum(sapply(batch_results, function(r) r$success))
    cat(sprintf(" done (%d/%d successful)\n", n_success, length(batch_results)))
  }

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  cat("\nMC phase complete in", round(as.numeric(t_elapsed), 1), "minutes\n")

  # --- Build results data frame ---
  mc_rows <- lapply(1:n_tasks, function(i) {
    r <- all_results[[i]]
    if (is.null(r) || !r$success) {
      return(data.frame(replicate = r$replicate, chl = r$chl,
                        success = FALSE, stringsAsFactors = FALSE))
    }

    row <- data.frame(
      replicate = r$replicate,
      chl       = r$chl,
      success   = TRUE,
      as.data.frame(t(r$param_values)),
      as.data.frame(t(r$metrics)),
      stringsAsFactors = FALSE
    )
    return(row)
  })

  mc_df <- dplyr::bind_rows(mc_rows)

  n_total_success <- sum(mc_df$success, na.rm = TRUE)
  n_total_failed  <- sum(!mc_df$success, na.rm = TRUE)
  cat("Total successful: ", n_total_success, " / Failed: ", n_total_failed, "\n")

  mc_output <- list(
    mc_df           = mc_df,
    param_sets      = param_sets,
    n_mc            = n_mc,
    chl_levels      = CHL_LEVELS,
    elapsed_minutes = as.numeric(t_elapsed)
  )

  saveRDS(mc_output, file.path(out_dir, "fish_sensitivity_mc_results.rds"))
  cat("Saved MC results to", file.path(out_dir, "fish_sensitivity_mc_results.rds"), "\n")

  # Clean up checkpoint
  if (file.exists(mc_checkpoint)) file.remove(mc_checkpoint)

  return(mc_output)
}


#' Calculate Partial Rank Correlation Coefficients
#'
#' Computes PRCC between each input parameter and a response variable,
#' with proper bootstrap confidence intervals (resampling raw data,
#' re-ranking, and re-computing partial correlations each replicate).
#'
#' @param mc_df Results data frame from MC analysis
#' @param response_var Name of the response variable column
#' @param chl_val CHL level to filter by (NULL for pooled analysis)
#' @param n_boot Number of bootstrap replicates
#' @return Data frame with PRCC, CI_lower, CI_upper, p_value per parameter
calculate_PRCC <- function(mc_df, response_var, chl_val = NULL, n_boot = 999) {

  # Filter by chl if specified
  if (!is.null(chl_val)) {
    mc_df <- mc_df[mc_df$chl == chl_val, ]
  }

  # Filter to successful runs with valid response
  valid <- mc_df[mc_df$success & !is.na(mc_df[[response_var]]), ]

  # Input columns (parameter draws: param.group format)
  input_cols <- grep("^(PPMR|FeedWidth|f_M|K_growth|repro_eff|ZSpre|ZSexp|Wmat)\\.",
                     names(valid), value = TRUE)

  if (nrow(valid) < 30 || length(input_cols) < 2) {
    warning("Insufficient data for PRCC (n = ", nrow(valid),
            ", p = ", length(input_cols), ")")
    return(NULL)
  }

  # Ensure all input columns are numeric
  for (col in input_cols) {
    valid[[col]] <- as.numeric(valid[[col]])
  }
  valid[[response_var]] <- as.numeric(valid[[response_var]])

  # Function to compute all PRCCs from a dataset
  compute_all_prccs <- function(data) {
    ranked <- as.data.frame(lapply(data[c(input_cols, response_var)], rank))

    prccs <- sapply(input_cols, function(ic) {
      others <- setdiff(input_cols, ic)
      formula_in  <- reformulate(others, ic)
      formula_out <- reformulate(others, response_var)

      tryCatch({
        resid_in  <- residuals(lm(formula_in,  data = ranked))
        resid_out <- residuals(lm(formula_out, data = ranked))
        cor(resid_in, resid_out)
      }, error = function(e) NA)
    })

    return(prccs)
  }

  # Point estimates
  prcc_point <- compute_all_prccs(valid)

  # Bootstrap: resample full data, re-rank, re-compute
  set.seed(SEED + 1)
  boot_matrix <- replicate(n_boot, {
    idx <- sample(nrow(valid), replace = TRUE)
    compute_all_prccs(valid[idx, ])
  })
  # boot_matrix: n_params x n_boot

  prcc_df <- data.frame(
    parameter = input_cols,
    PRCC      = prcc_point,
    CI_lower  = apply(boot_matrix, 1, quantile, 0.025, na.rm = TRUE),
    CI_upper  = apply(boot_matrix, 1, quantile, 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # p-value: proportion of bootstrap replicates crossing zero
  prcc_df$p_value <- sapply(1:nrow(prcc_df), function(i) {
    bv <- boot_matrix[i, ]
    bv <- bv[!is.na(bv)]
    if (length(bv) < 10) return(NA)
    2 * min(mean(bv > 0), mean(bv < 0))
  })

  prcc_df <- prcc_df[order(-abs(prcc_df$PRCC)), ]
  rownames(prcc_df) <- NULL

  return(prcc_df)
}


# =============================================================================
# 6. ANALYSIS & PLOTTING
# =============================================================================

#' Tornado diagram for OAT sensitivity
#'
#' @param sensitivity_df Output from calculate_oat_sensitivity()
#' @param response_var Response variable to plot
#' @param chl_val CHL level to filter (NULL = average across CHL levels)
#' @return ggplot object
plot_tornado <- function(sensitivity_df, response_var, chl_val = NULL) {

  plot_data <- sensitivity_df[sensitivity_df$response == response_var, ]

  if (!is.null(chl_val)) {
    plot_data <- plot_data[plot_data$chl == chl_val, ]
  } else {
    # Average across CHL levels
    plot_data <- plot_data %>%
      dplyr::group_by(parameter, response) %>%
      dplyr::summarise(rel_change = mean(rel_change, na.rm = TRUE),
                       .groups = "drop")
  }

  plot_data <- plot_data %>%
    dplyr::filter(parameter != "all_differentiated") %>%
    dplyr::arrange(abs(rel_change)) %>%
    dplyr::mutate(parameter = factor(parameter, levels = parameter))

  ggplot(plot_data, aes(x = rel_change, y = parameter, fill = rel_change > 0)) +
    geom_col(show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C")) +
    theme_bw(base_size = 12) +
    labs(
      x = paste("% Change in", response_var),
      y = "Differentiated Parameter",
      title = paste("OAT Sensitivity:", response_var),
      subtitle = if (!is.null(chl_val)) {
        paste("Chl =", chl_val, "mg/m3")
      } else {
        "Averaged across CHL levels"
      }
    )
}


#' PRCC dot-and-whisker plot
#'
#' @param prcc_df Output from calculate_PRCC()
#' @param response_var Response variable name (for title)
#' @param chl_val CHL level (for subtitle, NULL if pooled)
#' @return ggplot object
plot_PRCC <- function(prcc_df, response_var, chl_val = NULL) {

  plot_data <- prcc_df %>%
    dplyr::mutate(
      parameter   = factor(parameter, levels = parameter[order(abs(PRCC))]),
      significant = p_value < 0.05
    )

  ggplot(plot_data, aes(x = PRCC, y = parameter)) +
    geom_point(aes(colour = significant), size = 3) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = expression(p >= 0.05))
    ) +
    theme_bw(base_size = 12) +
    labs(
      x = "Partial Rank Correlation Coefficient",
      y = "Parameter",
      title = paste("Parameter Sensitivity:", response_var),
      subtitle = if (!is.null(chl_val)) {
        paste("Chl =", chl_val, "mg/m3")
      } else {
        "Pooled across CHL levels"
      },
      colour = "Significance"
    )
}


#' Steady-state assessment plot
#'
#' Shows biomass trend slopes for each fish group across configurations
#' and CHL levels. Dashed lines indicate the steady-state threshold.
#'
#' @param oat_output Output from run_oat_sensitivity()
#' @return ggplot object
plot_steady_state <- function(oat_output) {

  tasks   <- oat_output$tasks
  results <- oat_output$results

  diag_rows <- list()
  for (i in 1:nrow(tasks)) {
    if (!results[[i]]$success) next
    d <- results[[i]]$diagnostics
    d$config <- tasks$config_name[i]
    d$chl    <- tasks$chl[i]
    diag_rows <- c(diag_rows, list(d))
  }

  diag_df <- dplyr::bind_rows(diag_rows)

  # Focus on fish groups
  fish_diag <- diag_df[diag_df$type == "Fish", ]

  ggplot(fish_diag, aes(x = config, y = trend_slope, colour = group)) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = c(-TREND_THRESHOLD, TREND_THRESHOLD),
               linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    facet_wrap(~ chl, labeller = label_both) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Configuration",
      y = "Trend slope (log10 biomass / year)",
      title = "Steady-State Assessment: Fish Group Biomass Trends",
      subtitle = paste0("Final ", N_AVG_YEARS,
                        " years; dashed lines = +/-", TREND_THRESHOLD,
                        " threshold"),
      colour = "Fish Group"
    )
}


#' CV summary plot for all configurations
#'
#' @param oat_output Output from run_oat_sensitivity()
#' @return ggplot object
plot_cv_summary <- function(oat_output) {

  tasks   <- oat_output$tasks
  results <- oat_output$results

  diag_rows <- list()
  for (i in 1:nrow(tasks)) {
    if (!results[[i]]$success) next
    d <- results[[i]]$diagnostics
    d$config <- tasks$config_name[i]
    d$chl    <- tasks$chl[i]
    diag_rows <- c(diag_rows, list(d))
  }

  diag_df <- dplyr::bind_rows(diag_rows)
  fish_diag <- diag_df[diag_df$type == "Fish", ]

  ggplot(fish_diag, aes(x = config, y = cv, fill = group)) +
    geom_col(position = "dodge") +
    facet_wrap(~ chl, labeller = label_both) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Configuration",
      y = "Coefficient of Variation",
      title = "Biomass Variability: Fish Groups",
      subtitle = paste0("CV over final ", N_AVG_YEARS, " years"),
      fill = "Fish Group"
    )
}


# =============================================================================
# 7. MASTER RUNNER
# =============================================================================

#' Run the complete fish parameter sensitivity analysis
#'
#' Orchestrates: Smoke test -> Phase 1 (OAT) -> Phase 2 (MC with LHS).
#' Each phase checkpoints results to disk for robustness.
#'
#' @param run_mc Logical; whether to run Monte Carlo phase (default TRUE)
#' @param n_mc Number of MC replicates (default N_MC = 100)
#' @return List with all results, sensitivity indices, PRCC, and plots
run_full_sensitivity <- function(run_mc = TRUE, n_mc = N_MC) {

  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  ZooMSS Fish Parameter Sensitivity Analysis\n")
  cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(strrep("-", 60), "\n")
  cat("  Sim years:   ", SIM_YEARS, "\n")
  cat("  Avg years:   ", N_AVG_YEARS, " (final years for equilibrium)\n")
  cat("  CHL levels:  ", paste(CHL_LEVELS, collapse = ", "), " mg/m3\n")
  cat("  Workers:     ", n_workers, "\n")
  cat("  Min R_frac:  ", MIN_R_FRAC, "\n")
  if (run_mc) cat("  MC reps:     ", n_mc, " (LHS)\n")
  cat(strrep("=", 60), "\n\n")

  # Log system info
  cat("System: ", Sys.info()["sysname"], Sys.info()["release"], "\n")
  cat("R:      ", R.version.string, "\n")
  cat("zoomss: ", as.character(packageVersion("zoomss")), "\n\n")

  # ── Phase 0: Smoke test ──
  run_smoke_test()

  # ── Phase 1: OAT ──
  oat_output     <- run_oat_sensitivity()
  sensitivity_df <- calculate_oat_sensitivity(oat_output)

  cat("\n--- OAT Sensitivity Summary (averaged across CHL levels) ---\n")
  summary_avg <- sensitivity_df %>%
    dplyr::filter(parameter != "all_differentiated") %>%
    dplyr::group_by(parameter, response) %>%
    dplyr::summarise(mean_rel_change = mean(rel_change, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(abs(mean_rel_change)))
  print(as.data.frame(head(summary_avg, 24)), digits = 3)

  # ── Phase 2: MC (optional) ──
  mc_output    <- NULL
  prcc_results <- NULL

  if (run_mc) {
    mc_output <- run_mc_sensitivity(n_mc = n_mc)

    # PRCC analysis for key responses, per CHL level
    prcc_results <- list()
    for (chl_val in CHL_LEVELS) {
      key <- paste0("chl_", chl_val)
      cat("\nCalculating PRCC for chl =", chl_val, "...\n")
      prcc_results[[key]] <- list(
        total_fish     = calculate_PRCC(mc_output$mc_df, "total_fish", chl_val),
        SL_ratio       = calculate_PRCC(mc_output$mc_df, "SL_ratio", chl_val),
        spectrum_slope = calculate_PRCC(mc_output$mc_df, "spectrum_slope", chl_val)
      )
    }
  }

  # ── Generate plots ──
  cat("\nGenerating plots...\n")

  plot_list <- list()

  # Tornado diagrams (averaged across CHL)
  for (rv in c("total_fish", "SL_ratio", "spectrum_slope",
               "Fish_Small", "Fish_Med", "Fish_Large")) {
    plot_list[[paste0("tornado_", rv)]] <- plot_tornado(sensitivity_df, rv)
  }

  # Per-CHL tornado for total_fish
  for (chl_val in CHL_LEVELS) {
    plot_list[[paste0("tornado_total_fish_chl", chl_val)]] <-
      plot_tornado(sensitivity_df, "total_fish", chl_val)
  }

  # Steady-state assessment
  plot_list[["steady_state"]] <- plot_steady_state(oat_output)
  plot_list[["cv_summary"]]   <- plot_cv_summary(oat_output)

  # PRCC plots (if MC was run)
  if (!is.null(prcc_results)) {
    for (chl_key in names(prcc_results)) {
      chl_val <- as.numeric(sub("chl_", "", chl_key))
      for (rv in names(prcc_results[[chl_key]])) {
        prcc_df <- prcc_results[[chl_key]][[rv]]
        if (!is.null(prcc_df)) {
          plot_list[[paste0("prcc_", chl_key, "_", rv)]] <-
            plot_PRCC(prcc_df, rv, chl_val)
        }
      }
    }
  }

  # Save all plots to PDF
  n_plots <- length(plot_list)
  if (n_plots > 0) {
    pdf_path <- file.path(out_dir, "fish_sensitivity_plots.pdf")
    pdf(pdf_path, width = 12, height = 8)
    for (p in plot_list) {
      tryCatch(print(p), error = function(e) {
        cat("  Warning: could not render plot:", e$message, "\n")
      })
    }
    dev.off()
    cat("Saved", n_plots, "plots to", pdf_path, "\n")
  }

  # ── Final output ──
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  Sensitivity Analysis Complete\n")
  cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(strrep("=", 60), "\n")

  final_output <- list(
    oat_output     = oat_output,
    sensitivity_df = sensitivity_df,
    mc_output      = mc_output,
    prcc_results   = prcc_results,
    plots          = plot_list,
    config         = list(
      SIM_YEARS   = SIM_YEARS,
      N_AVG_YEARS = N_AVG_YEARS,
      CHL_LEVELS  = CHL_LEVELS,
      SST         = SST,
      MIN_R_FRAC  = MIN_R_FRAC,
      SEED        = SEED,
      n_workers   = n_workers
    )
  )

  saveRDS(final_output, file.path(out_dir, "fish_sensitivity_full_results.rds"))
  cat("Saved complete results to",
      file.path(out_dir, "fish_sensitivity_full_results.rds"), "\n")

  # Reset parallel plan
  plan(sequential)

  return(final_output)
}


# =============================================================================
# USAGE
# =============================================================================
#
# Quick start:
#   source("calibration/06_fish_param_sensitivity.R")
#
# 1. Smoke test only:
#   run_smoke_test()
#
# 2. OAT sensitivity only (faster, no MC):
#   results <- run_full_sensitivity(run_mc = FALSE)
#
# 3. Full analysis with Monte Carlo:
#   results <- run_full_sensitivity(run_mc = TRUE, n_mc = 100)
#
# 4. Access results:
#   results$sensitivity_df                     # OAT sensitivity table
#   results$plots$tornado_total_fish           # Tornado diagram
#   results$plots$steady_state                 # Trend assessment
#   results$prcc_results$chl_0.1$total_fish    # PRCC rankings
#   results$plots$prcc_chl_0.1_total_fish      # PRCC plot
#
# 5. Resume interrupted MC (automatic):
#   # If MC was interrupted, re-running run_mc_sensitivity() or
#   # run_full_sensitivity() will detect the checkpoint and resume.
#
# 6. Widen parameter ranges (after initial results):
#   # Modify param_ranges and re-run MC phase only:
#   mc_output <- run_mc_sensitivity(n_mc = 200)
#
