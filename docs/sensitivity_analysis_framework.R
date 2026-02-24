#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity Analysis Framework: Fish Group Parameter Differentiation
# =============================================================================
# Tests the effect of differentiating fish functional group parameters
# in ZooMSS vs the current uniform parameterisation.
#
# Design:
#   Phase 1: One-at-a-time (OAT) sensitivity - vary each parameter independently
#   Phase 2: Full differentiation vs uniform baseline
#   Phase 3: Monte Carlo uncertainty analysis with correlated parameter draws
#
# Outputs measured:
#   - Biomass by group (total and size-resolved)
#   - Size spectrum slope and intercept
#   - Trophic level distribution
#   - Relative biomass proportions (Fish_Small:Med:Large ratios)
#   - SSB, recruitment, reproductive output
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
# library(zoomss)  # Your package

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

# Environmental scenario for sensitivity runs
# Use a stable environment to isolate parameter effects
ENVIRO <- list(
  n_years = 50,           # Long enough for equilibrium
  dt = 0.01,
  base_sst = 15,          # Moderate temperature
  base_chl = 0.5,         # Moderate productivity
  seasonal = FALSE,       # Static environment for clean sensitivity
  isave = 100             # Save every 100 steps (= 1 year)
)

# Years to average for equilibrium metrics (from end of simulation)
N_AVG_YEARS <- 10

# Number of Monte Carlo replicates for Phase 3
N_MC <- 200

# =============================================================================
# 1. PARAMETER DEFINITIONS
# =============================================================================

# Baseline (current uniform) parameters
baseline_params <- list(
  PPMR       = c(Fish_Small = 100,   Fish_Med = 100,   Fish_Large = 100),
  W0         = c(Fish_Small = -3.0,  Fish_Med = -3.0,  Fish_Large = -3.0),
  FeedWidth  = c(Fish_Small = 1.3,   Fish_Med = 1.3,   Fish_Large = 1.3),
  f_M        = c(Fish_Small = 0.50,  Fish_Med = 0.50,  Fish_Large = 0.50),
  K_growth   = c(Fish_Small = 0.30,  Fish_Med = 0.30,  Fish_Large = 0.30),
  repro_eff  = c(Fish_Small = 0.001, Fish_Med = 0.001, Fish_Large = 0.001),
  ZSpre      = c(Fish_Small = 0.10,  Fish_Med = 0.10,  Fish_Large = 0.10),
  ZSexp      = c(Fish_Small = 0.30,  Fish_Med = 0.30,  Fish_Large = 0.30)
)

# Differentiated (literature-informed) parameters
# NOTE: f_M + K_growth + R_frac = 1, so these must be considered as a triplet.
# Implied R_frac: Small = 0.20, Med = 0.25, Large = 0.35
# Recruitment differentiation via R_frac × repro_eff:
#   Small: 0.20 × 0.003  = 6e-4  (moderate allocation, high survival)
#   Med:   0.25 × 0.001  = 2.5e-4 (balanced)
#   Large: 0.35 × 0.0005 = 1.75e-4 (high allocation, low survival)
differentiated_params <- list(
  PPMR       = c(Fish_Small = 100,    Fish_Med = 300,    Fish_Large = 1000),
  W0         = c(Fish_Small = -3.3,   Fish_Med = -3.0,   Fish_Large = -3.0),
  FeedWidth  = c(Fish_Small = 1.1,    Fish_Med = 1.4,    Fish_Large = 1.2),
  f_M        = c(Fish_Small = 0.50,   Fish_Med = 0.50,   Fish_Large = 0.45),
  K_growth   = c(Fish_Small = 0.30,   Fish_Med = 0.25,   Fish_Large = 0.20),
  repro_eff  = c(Fish_Small = 0.003,  Fish_Med = 0.001,  Fish_Large = 0.0005),
  ZSpre      = c(Fish_Small = 0.15,   Fish_Med = 0.10,   Fish_Large = 0.20),
  ZSexp      = c(Fish_Small = 0.30,   Fish_Med = 0.30,   Fish_Large = 0.50)
)

# Uncertainty ranges for Monte Carlo (log-uniform or normal as appropriate)
# Format: list of (distribution, param1, param2) per group
# IMPORTANT: f_M and K_growth ranges must respect R_frac = 1 - f_M - K_growth >= 0
# The draw_random_params() function enforces this constraint post-hoc.
param_ranges <- list(
  PPMR = list(
    Fish_Small = list(dist = "lognormal", meanlog = log(100),  sdlog = 0.3),
    Fish_Med   = list(dist = "lognormal", meanlog = log(300),  sdlog = 0.4),
    Fish_Large = list(dist = "lognormal", meanlog = log(1000), sdlog = 0.5)
  ),
  FeedWidth = list(
    Fish_Small = list(dist = "normal", mean = 1.1, sd = 0.15),
    Fish_Med   = list(dist = "normal", mean = 1.4, sd = 0.15),
    Fish_Large = list(dist = "normal", mean = 1.2, sd = 0.15)
  ),
  f_M = list(
    Fish_Small = list(dist = "normal", mean = 0.50, sd = 0.03),
    Fish_Med   = list(dist = "normal", mean = 0.50, sd = 0.03),
    Fish_Large = list(dist = "normal", mean = 0.45, sd = 0.03)
  ),
  K_growth = list(
    Fish_Small = list(dist = "normal", mean = 0.30, sd = 0.03),
    Fish_Med   = list(dist = "normal", mean = 0.25, sd = 0.03),
    Fish_Large = list(dist = "normal", mean = 0.20, sd = 0.03)
  ),
  repro_eff = list(
    Fish_Small = list(dist = "lognormal", meanlog = log(0.003),  sdlog = 0.5),
    Fish_Med   = list(dist = "lognormal", meanlog = log(0.001),  sdlog = 0.5),
    Fish_Large = list(dist = "lognormal", meanlog = log(0.0005), sdlog = 0.5)
  ),
  ZSpre = list(
    Fish_Small = list(dist = "lognormal", meanlog = log(0.15), sdlog = 0.3),
    Fish_Med   = list(dist = "lognormal", meanlog = log(0.10), sdlog = 0.3),
    Fish_Large = list(dist = "lognormal", meanlog = log(0.20), sdlog = 0.3)
  ),
  ZSexp = list(
    Fish_Small = list(dist = "normal", mean = 0.30, sd = 0.05),
    Fish_Med   = list(dist = "normal", mean = 0.30, sd = 0.05),
    Fish_Large = list(dist = "normal", mean = 0.50, sd = 0.10)
  )
)

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

#' Apply parameter set to Groups data frame
#' @param Groups Default Groups data frame
#' @param params Named list of parameter vectors (names match GroupInputs columns)
#' @return Modified Groups data frame
apply_params_to_groups <- function(Groups, params) {

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

  # Ensure energy budget closure: R_frac = 1 - f_M - K_growth
  # (This is validated in zoomss_params, but worth checking here)
  for (fg in fish_idx) {
    R_frac <- 1 - Groups$f_M[fg] - Groups$K_growth[fg]
    if (R_frac < 0) {
      warning("Energy budget violation for ", Groups$Species[fg],
              ": f_M + K_growth = ", Groups$f_M[fg] + Groups$K_growth[fg],
              " > 1. Adjusting K_growth.")
      Groups$K_growth[fg] <- 1 - Groups$f_M[fg] - 0.01
    }
  }

  return(Groups)
}


#' Run a single ZooMSS simulation and extract metrics
#' @param Groups Modified Groups data frame
#' @param enviro_config Environment configuration list
#' @return Named list of equilibrium metrics
run_and_extract <- function(Groups, enviro_config = ENVIRO) {

  # Create environment
  env_data <- createEnviroData(
    n_years = enviro_config$n_years,
    dt = enviro_config$dt,
    base_sst = enviro_config$base_sst,
    base_chl = enviro_config$base_chl,
    seasonal = enviro_config$seasonal
  )

  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)

  # Run model
  mdl <- suppressMessages(
    zoomss_model(input_params, Groups, isave = enviro_config$isave)
  )

  # Extract equilibrium metrics
  metrics <- extract_metrics(mdl, n_years = N_AVG_YEARS)

  return(metrics)
}


#' Extract key metrics from model output
#' @param mdl ZooMSS model output
#' @param n_years Years to average over
#' @return Named list of metrics
extract_metrics <- function(mdl, n_years = 10) {

  fish_grps <- mdl$param$fish_grps
  fish_names <- mdl$param$Groups$Species[fish_grps]

  # Biomass by group (averaged over final n_years)
  avg_biomass <- averageTimeSeries(mdl, "biomass", n_years = n_years)
  fish_biomass <- rowSums(avg_biomass[fish_grps, , drop = FALSE])
  names(fish_biomass) <- fish_names

  total_fish_biomass <- sum(fish_biomass)

  # Biomass ratios
  biomass_ratios <- fish_biomass / total_fish_biomass
  names(biomass_ratios) <- paste0("ratio_", fish_names)

  # Small:Large ratio (key ecosystem indicator)
  SL_ratio <- fish_biomass["Fish_Small"] / max(fish_biomass["Fish_Large"], 1e-30)

  # Size spectrum slope (fish component only)
  avg_abundance <- averageTimeSeries(mdl, "abundance", n_years = n_years)
  fish_abundance <- colSums(avg_abundance[fish_grps, , drop = FALSE])
  w_log10 <- mdl$param$w_log10

  # Fit log-log slope where fish have non-zero abundance
  fish_present <- fish_abundance > 0
  if (sum(fish_present) > 5) {
    lm_fit <- lm(log10(fish_abundance[fish_present]) ~ w_log10[fish_present])
    spectrum_slope <- coef(lm_fit)[2]
    spectrum_intercept <- coef(lm_fit)[1]
  } else {
    spectrum_slope <- NA
    spectrum_intercept <- NA
  }

  # Trophic levels
  tl <- extractTrophicLevels(mdl)
  n_t <- nrow(tl)
  avg_start <- max(1, n_t - round(n_years / (mdl$param$isave * mdl$param$dt)))
  avg_tl <- colMeans(tl[avg_start:n_t, , drop = FALSE])
  fish_tl <- avg_tl[fish_grps]
  names(fish_tl) <- paste0("TL_", fish_names)

  # Reproduction metrics (final saved time step)
  final_idx <- dim(mdl$SSB)[1]
  ssb <- mdl$SSB[final_idx, ]
  names(ssb) <- paste0("SSB_", fish_names)
  recruitment <- mdl$recruitment[final_idx, ]
  names(recruitment) <- paste0("Recruit_", fish_names)

  # Combine all metrics
  metrics <- c(
    fish_biomass,
    biomass_ratios,
    SL_ratio = unname(SL_ratio),
    spectrum_slope = unname(spectrum_slope),
    spectrum_intercept = unname(spectrum_intercept),
    fish_tl,
    ssb,
    recruitment,
    total_fish_biomass = total_fish_biomass
  )

  return(metrics)
}


#' Draw random parameter set from uncertainty distributions
#' @param param_ranges List of distribution specifications
#' @return Named list of parameter vectors (same structure as differentiated_params)
draw_random_params <- function(param_ranges) {

  drawn <- list()

  for (param_name in names(param_ranges)) {
    vals <- numeric(3)
    names(vals) <- c("Fish_Small", "Fish_Med", "Fish_Large")

    for (group_name in names(param_ranges[[param_name]])) {
      spec <- param_ranges[[param_name]][[group_name]]

      val <- switch(spec$dist,
        "normal"    = rnorm(1, mean = spec$mean, sd = spec$sd),
        "lognormal" = rlnorm(1, meanlog = spec$meanlog, sdlog = spec$sdlog),
        stop("Unknown distribution: ", spec$dist)
      )

      # Apply bounds
      if (param_name == "FeedWidth") val <- max(0.5, min(val, 2.5))
      if (param_name == "f_M")      val <- max(0.2, min(val, 0.8))
      if (param_name == "K_growth") val <- max(0.10, min(val, 0.50))
      if (param_name == "repro_eff") val <- max(1e-5, min(val, 0.01))
      if (param_name == "PPMR")     val <- max(10, min(val, 5000))
      if (param_name == "ZSpre")    val <- max(0.01, min(val, 1.0))
      if (param_name == "ZSexp")    val <- max(0.1, min(val, 1.0))

      vals[group_name] <- val
    }
    drawn[[param_name]] <- vals
  }

  # Enforce energy budget constraint: f_M + K_growth < 1 AND R_frac reasonable
  # Target R_frac ranges: Small ~0.15-0.30, Med ~0.15-0.35, Large ~0.25-0.45
  min_R_frac <- c(Fish_Small = 0.10, Fish_Med = 0.10, Fish_Large = 0.15)

  for (gn in c("Fish_Small", "Fish_Med", "Fish_Large")) {
    budget <- drawn$f_M[gn] + drawn$K_growth[gn]
    max_allowed <- 1.0 - min_R_frac[gn]

    if (budget > max_allowed) {
      # Scale both down proportionally to preserve their relative magnitudes
      scale <- (max_allowed - 0.01) / budget
      drawn$f_M[gn] <- drawn$f_M[gn] * scale
      drawn$K_growth[gn] <- drawn$K_growth[gn] * scale
    }
  }

  return(drawn)
}


# =============================================================================
# 3. PHASE 1: ONE-AT-A-TIME SENSITIVITY
# =============================================================================

run_OAT_sensitivity <- function() {

  cat("=== Phase 1: One-at-a-Time Sensitivity Analysis ===\n")

  # Get default groups
  Groups_default <- getGroups()

  # Run baseline
  cat("Running baseline (uniform parameters)...\n")
  Groups_baseline <- apply_params_to_groups(Groups_default, baseline_params)
  baseline_metrics <- run_and_extract(Groups_baseline)

  # For each parameter, run with only that parameter differentiated
  oat_results <- list(baseline = baseline_metrics)

  params_to_test <- c("PPMR", "FeedWidth", "f_M", "K_growth",
                       "repro_eff", "ZSpre", "ZSexp")

  for (param_name in params_to_test) {
    cat("Testing differentiation of:", param_name, "...\n")

    # Start from baseline, change only this parameter
    test_params <- baseline_params
    test_params[[param_name]] <- differentiated_params[[param_name]]

    Groups_test <- apply_params_to_groups(Groups_default, test_params)
    test_metrics <- run_and_extract(Groups_test)

    oat_results[[param_name]] <- test_metrics
  }

  # Also run fully differentiated
  cat("Running fully differentiated...\n")
  Groups_diff <- apply_params_to_groups(Groups_default, differentiated_params)
  diff_metrics <- run_and_extract(Groups_diff)
  oat_results[["all_differentiated"]] <- diff_metrics

  return(oat_results)
}


# =============================================================================
# 4. PHASE 2: SENSITIVITY INDICES
# =============================================================================

#' Calculate sensitivity indices from OAT results
#' @param oat_results Output from run_OAT_sensitivity()
#' @return Data frame with sensitivity metrics
calculate_sensitivity_indices <- function(oat_results) {

  baseline <- oat_results$baseline
  param_names <- setdiff(names(oat_results), c("baseline", "all_differentiated"))

  # Key response variables to track
  response_vars <- c("Fish_Small", "Fish_Med", "Fish_Large",
                      "SL_ratio", "spectrum_slope", "total_fish_biomass")

  results <- data.frame()

  for (pn in param_names) {
    test <- oat_results[[pn]]

    for (rv in response_vars) {
      if (rv %in% names(baseline) && rv %in% names(test)) {
        base_val <- baseline[[rv]]
        test_val <- test[[rv]]

        # Absolute change
        abs_change <- test_val - base_val

        # Relative change (%)
        rel_change <- if (abs(base_val) > 1e-30) {
          100 * (test_val - base_val) / base_val
        } else { NA }

        results <- rbind(results, data.frame(
          parameter = pn,
          response = rv,
          baseline_value = base_val,
          test_value = test_val,
          abs_change = abs_change,
          rel_change = rel_change,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(results)
}


#' Plot OAT sensitivity tornado diagram
plot_tornado <- function(sensitivity_df, response_var = "total_fish_biomass") {

  plot_data <- sensitivity_df %>%
    filter(response == response_var) %>%
    arrange(abs(rel_change)) %>%
    mutate(parameter = factor(parameter, levels = parameter))

  ggplot(plot_data, aes(x = rel_change, y = parameter,
                        fill = rel_change > 0)) +
    geom_col(show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C")) +
    theme_bw() +
    labs(x = paste("% Change in", response_var),
         y = "Differentiated Parameter",
         title = paste("OAT Sensitivity:", response_var),
         subtitle = "Effect of differentiating each parameter independently")
}


# =============================================================================
# 5. PHASE 3: MONTE CARLO UNCERTAINTY ANALYSIS
# =============================================================================

run_MC_sensitivity <- function(n_replicates = N_MC, parallel = TRUE) {

  cat("=== Phase 3: Monte Carlo Uncertainty Analysis ===\n")
  cat("Running", n_replicates, "replicates...\n")

  Groups_default <- getGroups()

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- min(parallel::detectCores() - 1, 8)
    cat("Using", n_cores, "cores\n")

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))

    # Export required objects to workers
    parallel::clusterExport(cl, c("Groups_default", "param_ranges",
                                   "draw_random_params", "apply_params_to_groups",
                                   "run_and_extract", "extract_metrics",
                                   "ENVIRO", "N_AVG_YEARS"),
                           envir = environment())

    # Load package on workers
    parallel::clusterEvalQ(cl, {
      library(zoomss)
      library(dplyr)
    })

    results <- parallel::parLapply(cl, 1:n_replicates, function(i) {
      tryCatch({
        params <- draw_random_params(param_ranges)
        Groups_mc <- apply_params_to_groups(Groups_default, params)
        metrics <- run_and_extract(Groups_mc)
        c(replicate = i, unlist(params), metrics)
      }, error = function(e) {
        c(replicate = i, error = e$message)
      })
    })
  } else {
    # Sequential fallback
    results <- list()
    pb <- txtProgressBar(min = 0, max = n_replicates, style = 3)

    for (i in 1:n_replicates) {
      tryCatch({
        params <- draw_random_params(param_ranges)
        Groups_mc <- apply_params_to_groups(Groups_default, params)
        metrics <- run_and_extract(Groups_mc)
        results[[i]] <- c(replicate = i, unlist(params), metrics)
      }, error = function(e) {
        results[[i]] <- c(replicate = i, error = e$message)
        cat("\nError in replicate", i, ":", e$message, "\n")
      })
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  # Combine into data frame
  mc_df <- bind_rows(lapply(results, function(x) as.data.frame(t(x))))

  # Convert numeric columns
  numeric_cols <- setdiff(names(mc_df), "error")
  mc_df[numeric_cols] <- lapply(mc_df[numeric_cols], as.numeric)

  return(mc_df)
}


#' Analyse Monte Carlo results: partial rank correlation coefficients
#' @param mc_df Output from run_MC_sensitivity()
#' @param response_var Name of response variable
#' @return Data frame of PRCCs with confidence intervals
calculate_PRCC <- function(mc_df, response_var = "total_fish_biomass") {

  # Input parameters (columns with group-specific values)
  input_cols <- grep("^(PPMR|FeedWidth|f_M|K_growth|repro_eff|ZSpre|ZSexp)\\.",
                     names(mc_df), value = TRUE)

  # Filter out failed runs
  valid <- mc_df[!is.na(mc_df[[response_var]]), ]

  if (nrow(valid) < 30) {
    warning("Too few valid replicates (", nrow(valid), ") for PRCC analysis")
    return(NULL)
  }

  # Rank-transform all variables
  ranked <- as.data.frame(lapply(valid[c(input_cols, response_var)], rank))

  # Calculate partial correlations
  prcc_results <- data.frame()

  for (input_col in input_cols) {
    # Partial correlation: correlate input with response, controlling for all other inputs
    other_inputs <- setdiff(input_cols, input_col)

    # Residualize both input and response against other inputs
    resid_input <- residuals(lm(as.formula(paste(input_col, "~",
                                                  paste(other_inputs, collapse = "+"))),
                                data = ranked))
    resid_response <- residuals(lm(as.formula(paste(response_var, "~",
                                                     paste(other_inputs, collapse = "+"))),
                                   data = ranked))

    prcc <- cor(resid_input, resid_response)

    # Bootstrap CI
    boot_prccs <- replicate(1000, {
      idx <- sample(length(resid_input), replace = TRUE)
      cor(resid_input[idx], resid_response[idx])
    })

    prcc_results <- rbind(prcc_results, data.frame(
      parameter = input_col,
      PRCC = prcc,
      CI_lower = quantile(boot_prccs, 0.025),
      CI_upper = quantile(boot_prccs, 0.975),
      p_value = 2 * min(mean(boot_prccs > 0), mean(boot_prccs < 0)),
      stringsAsFactors = FALSE
    ))
  }

  prcc_results <- prcc_results %>%
    arrange(desc(abs(PRCC)))

  return(prcc_results)
}


#' Plot PRCC results
plot_PRCC <- function(prcc_df, response_var = "total_fish_biomass") {

  plot_data <- prcc_df %>%
    mutate(parameter = factor(parameter, levels = parameter[order(abs(PRCC))])) %>%
    mutate(significant = p_value < 0.05)

  ggplot(plot_data, aes(x = PRCC, y = parameter)) +
    geom_point(aes(colour = significant), size = 3) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
                       labels = c("p ≥ 0.05", "p < 0.05")) +
    theme_bw() +
    labs(x = "Partial Rank Correlation Coefficient",
         y = "Parameter",
         title = paste("Parameter Sensitivity:", response_var),
         subtitle = paste("PRCC from", nrow(prcc_df), "Monte Carlo replicates"),
         colour = "Significance")
}


# =============================================================================
# 6. MASTER RUNNER
# =============================================================================

run_full_sensitivity <- function(run_mc = TRUE, n_mc = N_MC) {

  cat("╔══════════════════════════════════════════════════╗\n")
  cat("║  ZooMSS Fish Parameter Sensitivity Analysis     ║\n")
  cat("╚══════════════════════════════════════════════════╝\n\n")

  # Phase 1: OAT
  oat_results <- run_OAT_sensitivity()
  sensitivity_df <- calculate_sensitivity_indices(oat_results)

  cat("\n--- OAT Sensitivity Summary ---\n")
  print(sensitivity_df %>%
          arrange(desc(abs(rel_change))) %>%
          head(20),
        digits = 3)

  # Save OAT plots
  response_vars <- c("total_fish_biomass", "Fish_Small", "Fish_Large", "SL_ratio")
  oat_plots <- lapply(response_vars, function(rv) plot_tornado(sensitivity_df, rv))
  names(oat_plots) <- response_vars

  # Phase 3: Monte Carlo (optional, computationally expensive)
  mc_results <- NULL
  prcc_results <- NULL
  mc_plots <- NULL

  if (run_mc) {
    mc_results <- run_MC_sensitivity(n_replicates = n_mc)

    # PRCC analysis for key responses
    prcc_biomass <- calculate_PRCC(mc_results, "total_fish_biomass")
    prcc_ratio   <- calculate_PRCC(mc_results, "SL_ratio")
    prcc_slope   <- calculate_PRCC(mc_results, "spectrum_slope")

    prcc_results <- list(
      total_fish_biomass = prcc_biomass,
      SL_ratio = prcc_ratio,
      spectrum_slope = prcc_slope
    )

    mc_plots <- list(
      prcc_biomass = plot_PRCC(prcc_biomass, "total_fish_biomass"),
      prcc_ratio   = plot_PRCC(prcc_ratio, "SL_ratio"),
      prcc_slope   = plot_PRCC(prcc_slope, "spectrum_slope")
    )
  }

  # Return everything
  return(list(
    oat_results = oat_results,
    sensitivity_df = sensitivity_df,
    oat_plots = oat_plots,
    mc_results = mc_results,
    prcc_results = prcc_results,
    mc_plots = mc_plots
  ))
}

# =============================================================================
# USAGE
# =============================================================================
# Quick OAT-only run:
#   results <- run_full_sensitivity(run_mc = FALSE)
#
# Full analysis with Monte Carlo:
#   results <- run_full_sensitivity(run_mc = TRUE, n_mc = 200)
#
# Access results:
#   results$sensitivity_df          # OAT sensitivity table
#   results$oat_plots$SL_ratio      # Tornado diagram
#   results$prcc_results$total_fish_biomass  # PRCC rankings
#   results$mc_plots$prcc_biomass   # PRCC plot
