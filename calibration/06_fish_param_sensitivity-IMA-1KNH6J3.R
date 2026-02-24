# =============================================================================
# 06_fish_param_sensitivity.R
# One-at-a-Time (OAT) Sensitivity Analysis for Fish Group Parameters
# =============================================================================
#
# Purpose: Test how perturbations of each fish-specific parameter affect model
#          stability and biomass patterns across a chlorophyll gradient.
#          This helps identify which parameters are most influential before
#          committing to differentiation or broader exploration.
#
# Design:
#   - 8 fish parameters tested: PPMR, FeedWidth, f_M, K_growth, repro_eff,
#     ZSpre, ZSexp, Wmat
#   - Each parameter perturbed ±10% (or ±20% for log-scale params) from
#     baseline, one at a time, applied to ALL fish groups simultaneously
#   - Tested at 3 chlorophyll levels: 0.01, 0.1, 1.0 mg/m³
#   - 250-year simulations; stability assessed over final 50 years
#   - Biomass CV, linear trend, and coexistence used as diagnostics
#
# Constraints:
#   - R_frac = 1 - f_M - K_growth >= 0.15 (enforced in perturbation design)
#   - W0 not perturbed (fixed)
#   - Wmat must remain within [W0, Wmax] for each group
#
# Prerequisites:
#   - Run calibration/00_smoke_test.R first
#   - devtools installed
#   - future.apply and parallelly installed
#
# Output:
#   - calibration/oat_sensitivity_results.rds  (full results)
#   - calibration/oat_sensitivity_summary.csv  (summary table)
#   - calibration/oat_sensitivity_plots.pdf    (diagnostic plots)
#
# Expected runtime: see smoke test estimate (typically 30-90 min)
# =============================================================================

cat("╔══════════════════════════════════════════════════════╗\n")
cat("║  Fish Parameter OAT Sensitivity Analysis            ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 0. SETUP
# =============================================================================

# Load dev-branch package
devtools::load_all(".")

library(future.apply)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Parallel backend ──
n_cores <- min(14, parallelly::availableCores() - 2)
plan(multisession, workers = n_cores)
cat("Parallel backend: future.apply with", n_cores, "workers\n\n")

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

# Environmental settings
CHL_LEVELS   <- c(0.01, 0.1, 1.0)  # Chlorophyll gradient (mg m⁻³)
SST_CONSTANT <- 15                   # Constant temperature (°C)
SIM_YEARS    <- 250                  # Total simulation length
DT           <- 0.01                 # Time step (years)
ISAVE        <- 100                  # Save every 100 steps (= 1 year)
ANALYSIS_YEARS <- 50                 # Final N years for stability analysis
MIN_R_FRAC   <- 0.15                 # Minimum allowed reproduction fraction

# Derived
ANALYSIS_START_YEAR <- SIM_YEARS - ANALYSIS_YEARS  # Year 200

cat("Configuration:\n")
cat("  Chl levels:     ", paste(CHL_LEVELS, collapse = ", "), "mg/m³\n")
cat("  SST:             ", SST_CONSTANT, "°C\n")
cat("  Sim length:      ", SIM_YEARS, "years\n")
cat("  Analysis window:  years", ANALYSIS_START_YEAR, "-", SIM_YEARS, "\n")
cat("  Min R_frac:      ", MIN_R_FRAC, "\n\n")

# =============================================================================
# 2. LOAD BASELINE GROUPS AND DEFINE PERTURBATIONS
# =============================================================================

Groups_default <- getGroups()
fish_idx <- which(Groups_default$Type == "Fish")
fish_names <- Groups_default$Species[fish_idx]

cat("Baseline fish parameters:\n")
cat(sprintf("  %-12s  PPMR=%4.0f  FW=%.2f  f_M=%.2f  Kg=%.2f  R_frac=%.2f  re=%.4f  ZSpre=%.2f  ZSexp=%.2f  Wmat=%+.1f\n",
            fish_names,
            Groups_default$PPMR[fish_idx],
            Groups_default$FeedWidth[fish_idx],
            Groups_default$f_M[fish_idx],
            Groups_default$K_growth[fish_idx],
            1 - Groups_default$f_M[fish_idx] - Groups_default$K_growth[fish_idx],
            Groups_default$repro_eff[fish_idx],
            Groups_default$ZSpre[fish_idx],
            Groups_default$ZSexp[fish_idx],
            Groups_default$Wmat[fish_idx]))

# ── Define perturbation table ──
# Each row: parameter name, column in Groups, low value(s), high value(s)
# Values applied UNIFORMLY to all fish groups unless noted.
# For parameters already differentiated (ZSpre, ZSexp, Wmat), perturbations
# are multiplicative or additive relative to each group's baseline.

# Helper: compute perturbed fish values, enforcing R_frac >= MIN_R_FRAC
build_perturbation_table <- function(Groups, fish_idx) {

  baseline <- list(
    PPMR      = Groups$PPMR[fish_idx],
    FeedWidth = Groups$FeedWidth[fish_idx],
    f_M       = Groups$f_M[fish_idx],
    K_growth  = Groups$K_growth[fish_idx],
    repro_eff = Groups$repro_eff[fish_idx],
    ZSpre     = Groups$ZSpre[fish_idx],
    ZSexp     = Groups$ZSexp[fish_idx],
    Wmat      = Groups$Wmat[fish_idx]
  )

  # ±10% multiplicative for most params; ±20% for repro_eff (very uncertain)
  # Additive ±0.3 for Wmat (log10 scale)
  perturbations <- list(
    PPMR = list(
      low  = round(baseline$PPMR * 0.9),
      high = round(baseline$PPMR * 1.1)
    ),
    FeedWidth = list(
      low  = round(baseline$FeedWidth * 0.9, 2),
      high = round(baseline$FeedWidth * 1.1, 2)
    ),
    f_M = list(
      # f_M - 0.05: R_frac increases (safe)
      # f_M + 0.04: R_frac = 1 - 0.54 - K_growth, check per group
      low  = baseline$f_M - 0.05,
      high = pmin(baseline$f_M + 0.04,
                  1 - baseline$K_growth - MIN_R_FRAC)
    ),
    K_growth = list(
      low  = baseline$K_growth - 0.03,
      high = pmin(baseline$K_growth + 0.03,
                  1 - baseline$f_M - MIN_R_FRAC)
    ),
    repro_eff = list(
      low  = baseline$repro_eff * 0.5,
      high = baseline$repro_eff * 2.0
    ),
    ZSpre = list(
      low  = round(baseline$ZSpre * 0.8, 3),
      high = round(baseline$ZSpre * 1.2, 3)
    ),
    ZSexp = list(
      low  = round(baseline$ZSexp * 0.9, 2),
      high = round(baseline$ZSexp * 1.1, 2)
    ),
    Wmat = list(
      low  = baseline$Wmat - 0.3,
      high = baseline$Wmat + 0.3
    )
  )

  # ── Validate R_frac for f_M and K_growth perturbations ──
  for (pname in c("f_M", "K_growth")) {
    for (dir in c("low", "high")) {
      test_f_M <- if (pname == "f_M") perturbations[[pname]][[dir]] else baseline$f_M
      test_Kg  <- if (pname == "K_growth") perturbations[[pname]][[dir]] else baseline$K_growth
      test_R   <- 1 - test_f_M - test_Kg
      if (any(test_R < MIN_R_FRAC)) {
        bad <- which(test_R < MIN_R_FRAC)
        stop("R_frac constraint violated for ", pname, " (", dir, ") in group(s): ",
             paste(fish_names[bad], collapse = ", "),
             " | R_frac = ", paste(round(test_R[bad], 3), collapse = ", "))
      }
    }
  }

  # ── Validate Wmat within [W0, Wmax] ──
  for (dir in c("low", "high")) {
    wmat_test <- perturbations$Wmat[[dir]]
    w0_vals   <- Groups$W0[fish_idx]
    wmax_vals <- Groups$Wmax[fish_idx]
    if (any(wmat_test < w0_vals) || any(wmat_test > wmax_vals)) {
      # Clamp rather than fail
      perturbations$Wmat[[dir]] <- pmax(w0_vals + 0.1,
                                         pmin(wmax_vals - 0.1, wmat_test))
      cat("  NOTE: Wmat (", dir, ") clamped to [W0+0.1, Wmax-0.1]\n")
    }
  }

  return(list(baseline = baseline, perturbations = perturbations))
}

pert <- build_perturbation_table(Groups_default, fish_idx)

# Print perturbation summary
cat("\nPerturbation summary:\n")
for (pname in names(pert$perturbations)) {
  bl  <- pert$baseline[[pname]]
  lo  <- pert$perturbations[[pname]]$low
  hi  <- pert$perturbations[[pname]]$high
  cat(sprintf("  %-10s  baseline: [%s]  low: [%s]  high: [%s]\n",
              pname,
              paste(round(bl, 4), collapse = ", "),
              paste(round(lo, 4), collapse = ", "),
              paste(round(hi, 4), collapse = ", ")))
}

# =============================================================================
# 3. BUILD DESIGN MATRIX
# =============================================================================

# Each run is defined by: (param_name, direction, chl_level)
design <- expand.grid(
  param     = names(pert$perturbations),
  direction = c("low", "high"),
  chl       = CHL_LEVELS,
  stringsAsFactors = FALSE
)

# Add baseline runs (one per Chl level)
baseline_rows <- data.frame(
  param     = "baseline",
  direction = "none",
  chl       = CHL_LEVELS,
  stringsAsFactors = FALSE
)

design <- rbind(baseline_rows, design)
design$run_id <- seq_len(nrow(design))

cat("\nDesign matrix:", nrow(design), "total runs\n")
cat("  Baselines:", sum(design$param == "baseline"), "\n")
cat("  Perturbed:", sum(design$param != "baseline"), "\n")
cat("  Batches (~", n_cores, "cores):", ceiling(nrow(design) / n_cores), "\n\n")

# =============================================================================
# 4. HELPER FUNCTIONS
# =============================================================================

#' Apply a single perturbation to the Groups data frame
#' @param Groups Default groups data frame
#' @param fish_idx Indices of fish groups
#' @param param_name Which parameter to perturb (or "baseline" for no change)
#' @param direction "low", "high", or "none"
#' @param perturbations List of perturbation values
#' @return Modified Groups data frame
apply_perturbation <- function(Groups, fish_idx, param_name, direction, perturbations) {

  G <- Groups  # copy

  if (param_name == "baseline" || direction == "none") {
    return(G)
  }

  vals <- perturbations[[param_name]][[direction]]
  G[[param_name]][fish_idx] <- vals

  return(G)
}


#' Run a single simulation and return stability metrics
#' @param Groups Modified groups data frame
#' @param chl Chlorophyll level
#' @return Named list of stability and biomass metrics
run_single_sim <- function(Groups, chl) {

  # Create static environment
  env_data <- createEnviroData(
    n_years  = SIM_YEARS,
    dt       = DT,
    base_sst = SST_CONSTANT,
    base_chl = chl,
    seasonal = FALSE
  )
  input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)

  # Run model
  mdl <- suppressMessages(
    zoomss_model(input_params, Groups, isave = ISAVE)
  )

  # Extract stability metrics from final ANALYSIS_YEARS

  metrics <- extract_stability_metrics(mdl)

  return(metrics)
}


#' Extract stability and biomass metrics from model output
#' @param mdl ZooMSS model output
#' @return Named list of metrics
extract_stability_metrics <- function(mdl) {

  # ── Dimensions ──
  time_vec   <- mdl$time
  n_saved    <- length(time_vec)
  fish_grps  <- mdl$param$fish_grps
  zoo_grps   <- mdl$param$zoo_grps
  fish_names <- mdl$param$Groups$Species[fish_grps]
  zoo_names  <- mdl$param$Groups$Species[zoo_grps]
  all_names  <- mdl$param$Groups$Species

  # ── Identify analysis window (final ANALYSIS_YEARS) ──
  dt_saved <- mdl$param$isave * mdl$param$dt  # years per saved step
  n_analysis_steps <- round(ANALYSIS_YEARS / dt_saved)
  idx_start <- max(1, n_saved - n_analysis_steps + 1)
  idx_end   <- n_saved
  analysis_idx <- idx_start:idx_end

  # ── Annual total biomass per group (summed across size classes) ──
  # biomass is [time × groups × size]
  biomass_ts <- apply(mdl$biomass[analysis_idx, , , drop = FALSE], c(1, 2), sum)
  # biomass_ts is [time × groups]
  colnames(biomass_ts) <- all_names

  # ── Fish and zoo biomass time series ──
  fish_bio_ts <- biomass_ts[, fish_grps, drop = FALSE]
  zoo_bio_ts  <- biomass_ts[, zoo_grps, drop = FALSE]

  # ── METRIC 1: Biomass coefficient of variation (CV) ──
  # CV = sd / mean; high CV with no trend = oscillating steady state
  fish_cv <- apply(fish_bio_ts, 2, function(x) {
    if (mean(x) > 0) sd(x) / mean(x) else NA
  })
  names(fish_cv) <- paste0("CV_", fish_names)

  zoo_cv <- apply(zoo_bio_ts, 2, function(x) {
    if (mean(x) > 0) sd(x) / mean(x) else NA
  })
  names(zoo_cv) <- paste0("CV_", zoo_names)

  # ── METRIC 2: Linear trend in log-biomass ──
  # Fit lm(log(biomass) ~ year) for each group
  # Slope near zero = no trend; negative = declining; positive = increasing
  year_vec <- seq_len(nrow(fish_bio_ts))

  fish_trend_slope  <- numeric(length(fish_grps))
  fish_trend_pvalue <- numeric(length(fish_grps))
  names(fish_trend_slope)  <- paste0("trend_slope_", fish_names)
  names(fish_trend_pvalue) <- paste0("trend_pval_", fish_names)

  for (j in seq_along(fish_grps)) {
    bio_j <- fish_bio_ts[, j]
    if (all(bio_j > 0)) {
      fit <- lm(log(bio_j) ~ year_vec)
      fish_trend_slope[j]  <- coef(fit)[2]
      fish_trend_pvalue[j] <- summary(fit)$coefficients[2, 4]
    } else {
      fish_trend_slope[j]  <- NA
      fish_trend_pvalue[j] <- NA
    }
  }

  zoo_trend_slope  <- numeric(length(zoo_grps))
  zoo_trend_pvalue <- numeric(length(zoo_grps))
  names(zoo_trend_slope)  <- paste0("trend_slope_", zoo_names)
  names(zoo_trend_pvalue) <- paste0("trend_pval_", zoo_names)

  for (j in seq_along(zoo_grps)) {
    bio_j <- zoo_bio_ts[, j]
    if (all(bio_j > 0)) {
      fit <- lm(log(bio_j) ~ year_vec)
      zoo_trend_slope[j]  <- coef(fit)[2]
      zoo_trend_pvalue[j] <- summary(fit)$coefficients[2, 4]
    } else {
      zoo_trend_slope[j]  <- NA
      zoo_trend_pvalue[j] <- NA
    }
  }

  # ── METRIC 3: Coexistence check ──
  # A group "persists" if its mean biomass in the analysis window > 1e-20
  # and it has non-zero biomass at the final time step
  fish_mean_bio <- colMeans(fish_bio_ts)
  names(fish_mean_bio) <- paste0("mean_bio_", fish_names)

  fish_final_bio <- fish_bio_ts[nrow(fish_bio_ts), ]
  names(fish_final_bio) <- paste0("final_bio_", fish_names)

  fish_persists <- as.numeric(fish_mean_bio > 1e-20 & fish_final_bio > 1e-20)
  names(fish_persists) <- paste0("persists_", fish_names)

  zoo_mean_bio <- colMeans(zoo_bio_ts)
  names(zoo_mean_bio) <- paste0("mean_bio_", zoo_names)

  zoo_persists <- as.numeric(zoo_mean_bio > 1e-20)
  names(zoo_persists) <- paste0("persists_", zoo_names)

  n_fish_alive <- sum(fish_persists)
  n_zoo_alive  <- sum(zoo_persists)

  # ── METRIC 4: Biomass ratios ──
  total_fish_bio <- sum(fish_mean_bio)
  total_zoo_bio  <- sum(zoo_mean_bio)

  fish_ratios <- fish_mean_bio / max(total_fish_bio, 1e-30)
  names(fish_ratios) <- paste0("ratio_", fish_names)

  SL_ratio <- fish_mean_bio[grep("Small", fish_names)] /
    max(fish_mean_bio[grep("Large", fish_names)], 1e-30)
  names(SL_ratio) <- "SL_ratio"

  zoo_fish_ratio <- total_zoo_bio / max(total_fish_bio, 1e-30)

  # ── METRIC 5: Steady-state classification ──
  # A group is in "oscillating steady state" if:
  #   (a) it persists, AND
  #   (b) no significant trend (p > 0.05 or |slope| < 0.001 per year)
  fish_steady <- as.numeric(
    fish_persists > 0 &
    (fish_trend_pvalue > 0.05 | abs(fish_trend_slope) < 0.001)
  )
  names(fish_steady) <- paste0("steady_", fish_names)

  # ── Combine all metrics ──
  metrics <- c(
    fish_cv,
    zoo_cv,
    fish_trend_slope,
    fish_trend_pvalue,
    zoo_trend_slope,
    zoo_trend_pvalue,
    fish_mean_bio,
    fish_final_bio,
    fish_persists,
    zoo_mean_bio,
    zoo_persists,
    fish_ratios,
    SL_ratio,
    total_fish_bio = unname(total_fish_bio),
    total_zoo_bio  = unname(total_zoo_bio),
    zoo_fish_ratio = unname(zoo_fish_ratio),
    n_fish_alive = n_fish_alive,
    n_zoo_alive  = n_zoo_alive,
    fish_steady
  )

  return(metrics)
}

# =============================================================================
# 5. RUN ALL SIMULATIONS IN PARALLEL
# =============================================================================

cat("Starting", nrow(design), "simulations...\n")
cat("Estimated runtime: see smoke test (calibration/00_smoke_test.R)\n\n")

t_start <- Sys.time()

results_list <- future_lapply(seq_len(nrow(design)), function(i) {

  row <- design[i, ]

  tryCatch({
    # Apply perturbation
    G <- apply_perturbation(
      Groups     = Groups_default,
      fish_idx   = fish_idx,
      param_name = row$param,
      direction  = row$direction,
      perturbations = pert$perturbations
    )

    # Run simulation
    metrics <- run_single_sim(G, chl = row$chl)

    # Return as named vector
    c(run_id   = row$run_id,
      param    = row$param,
      direction = row$direction,
      chl      = row$chl,
      status   = "success",
      metrics)

  }, error = function(e) {
    c(run_id   = row$run_id,
      param    = row$param,
      direction = row$direction,
      chl      = row$chl,
      status   = "error",
      error_msg = conditionMessage(e))
  })

}, future.seed = TRUE)

t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
cat("\nAll simulations completed in", round(t_elapsed, 1), "minutes.\n\n")

# =============================================================================
# 6. ASSEMBLE RESULTS
# =============================================================================

cat("Assembling results...\n")

# Separate successes and failures
success_idx <- sapply(results_list, function(x) x["status"] == "success")
n_success <- sum(success_idx)
n_fail    <- sum(!success_idx)

cat("  Successful:", n_success, "/", nrow(design), "\n")
if (n_fail > 0) {
  cat("  FAILED:", n_fail, "runs:\n")
  for (r in results_list[!success_idx]) {
    cat("    Run", r["run_id"], "(", r["param"], ",", r["direction"],
        ", chl=", r["chl"], "):", r["error_msg"], "\n")
  }
}

# Build data frame from successful runs
results_df <- bind_rows(lapply(results_list[success_idx], function(x) {
  as.data.frame(t(x), stringsAsFactors = FALSE)
}))

# Convert numeric columns (everything except param, direction, status)
char_cols <- c("param", "direction", "status")
num_cols  <- setdiff(names(results_df), char_cols)
results_df[num_cols] <- lapply(results_df[num_cols], as.numeric)

# =============================================================================
# 7. SENSITIVITY INDICES
# =============================================================================

cat("\nCalculating sensitivity indices...\n")

# For each (param, direction, chl), compute % change from baseline
baselines <- results_df %>% filter(param == "baseline")
perturbed <- results_df %>% filter(param != "baseline")

# Key response variables for sensitivity ranking
response_vars <- c("total_fish_bio", "total_zoo_bio", "SL_ratio",
                    paste0("mean_bio_", fish_names),
                    paste0("CV_", fish_names),
                    "n_fish_alive")

sensitivity_table <- perturbed %>%
  rowwise() %>%
  mutate(
    baseline_row = list(baselines[baselines$chl == chl, ]),
    .keep = "all"
  ) %>%
  ungroup()

# Calculate relative change for each response variable
sensitivity_long <- data.frame()

for (rv in response_vars) {
  if (!rv %in% names(perturbed)) next

  for (i in seq_len(nrow(perturbed))) {
    row <- perturbed[i, ]
    bl  <- baselines[baselines$chl == row$chl, ]

    if (nrow(bl) == 0) next

    base_val <- bl[[rv]]
    test_val <- row[[rv]]

    abs_change <- test_val - base_val
    rel_change <- if (!is.na(base_val) && abs(base_val) > 1e-30) {
      100 * (test_val - base_val) / base_val
    } else { NA }

    sensitivity_long <- rbind(sensitivity_long, data.frame(
      param      = row$param,
      direction  = row$direction,
      chl        = row$chl,
      response   = rv,
      base_val   = base_val,
      test_val   = test_val,
      abs_change = abs_change,
      rel_change = rel_change,
      stringsAsFactors = FALSE
    ))
  }
}

# ── Sensitivity ranking: max |rel_change| across directions and chl levels ──
sensitivity_rank <- sensitivity_long %>%
  group_by(param, response) %>%
  summarise(
    max_abs_rel_change = max(abs(rel_change), na.rm = TRUE),
    mean_abs_rel_change = mean(abs(rel_change), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(response, desc(max_abs_rel_change))

cat("\n--- Top Sensitivities (max |% change| across Chl and direction) ---\n")
print(
  sensitivity_rank %>%
    filter(response %in% c("total_fish_bio", "SL_ratio", "n_fish_alive")) %>%
    head(30),
  n = 30
)

# =============================================================================
# 8. STABILITY SUMMARY
# =============================================================================

cat("\n--- Stability Summary (baseline runs) ---\n")

for (i in seq_len(nrow(baselines))) {
  row <- baselines[i, ]
  cat(sprintf("\nChl = %.2f mg/m³:\n", row$chl))

  for (fn in fish_names) {
    cv_col    <- paste0("CV_", fn)
    trend_col <- paste0("trend_slope_", fn)
    pval_col  <- paste0("trend_pval_", fn)
    pers_col  <- paste0("persists_", fn)
    stdy_col  <- paste0("steady_", fn)

    cv_val    <- row[[cv_col]]
    trend_val <- row[[trend_col]]
    pval_val  <- row[[pval_col]]
    pers_val  <- row[[pers_col]]
    stdy_val  <- row[[stdy_col]]

    status <- if (pers_val == 0) {
      "EXTINCT"
    } else if (stdy_val == 1) {
      "STEADY"
    } else {
      "TRENDING"
    }

    cat(sprintf("  %-12s  CV=%.4f  trend=%.2e  p=%.3f  [%s]\n",
                fn, cv_val, trend_val, pval_val, status))
  }
}

# =============================================================================
# 9. PLOTTING
# =============================================================================

cat("\nGenerating diagnostic plots...\n")

# ── 9a. Tornado plots: % change from baseline for key responses ──
plot_tornado <- function(sens_df, response_var, title_suffix = "") {

  plot_data <- sens_df %>%
    filter(response == response_var) %>%
    group_by(param, chl) %>%
    summarise(
      max_change = rel_change[which.max(abs(rel_change))],
      .groups = "drop"
    ) %>%
    mutate(
      param = factor(param),
      chl_label = paste0("Chl = ", chl, " mg/m³")
    )

  if (nrow(plot_data) == 0) return(NULL)

  ggplot(plot_data, aes(x = max_change, y = reorder(param, abs(max_change)),
                        fill = max_change > 0)) +
    geom_col(show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    facet_wrap(~chl_label) +
    scale_fill_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C")) +
    theme_bw(base_size = 11) +
    labs(
      x = paste("% Change in", response_var),
      y = "Parameter",
      title = paste("OAT Sensitivity:", response_var, title_suffix)
    )
}

# ── 9b. Stability heatmap ──
plot_stability_heatmap <- function(results) {

  # Extract steady state flags for fish groups
  steady_cols <- grep("^steady_", names(results), value = TRUE)
  cv_cols     <- grep("^CV_Fish", names(results), value = TRUE)

  if (length(steady_cols) == 0) return(NULL)

  stability_data <- results %>%
    select(param, direction, chl, all_of(steady_cols)) %>%
    pivot_longer(cols = all_of(steady_cols),
                 names_to = "group", values_to = "steady") %>%
    mutate(
      group = gsub("steady_", "", group),
      label = paste(param, direction, sep = "_"),
      chl_label = paste0("Chl=", chl)
    )

  ggplot(stability_data, aes(x = chl_label, y = label, fill = factor(steady))) +
    geom_tile(colour = "white") +
    facet_wrap(~group) +
    scale_fill_manual(
      values = c("0" = "#E41A1C", "1" = "#4DAF4A", "NA" = "grey80"),
      labels = c("0" = "Unstable/Extinct", "1" = "Steady State"),
      na.value = "grey80",
      name = "Status"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 7)) +
    labs(x = "Chlorophyll", y = "Parameter Perturbation",
         title = "Fish Group Stability Across Perturbations")
}

# Generate and save plots
tryCatch({
  pdf("calibration/oat_sensitivity_plots.pdf", width = 12, height = 8)

  # Tornado plots for key responses
  for (rv in c("total_fish_bio", "SL_ratio",
               paste0("mean_bio_", fish_names))) {
    p <- plot_tornado(sensitivity_long, rv)
    if (!is.null(p)) print(p)
  }

  # CV tornado
  for (fn in fish_names) {
    p <- plot_tornado(sensitivity_long, paste0("CV_", fn),
                      title_suffix = "(biomass variability)")
    if (!is.null(p)) print(p)
  }

  # Stability heatmap
  p_heat <- plot_stability_heatmap(results_df)
  if (!is.null(p_heat)) print(p_heat)

  dev.off()
  cat("  Plots saved to calibration/oat_sensitivity_plots.pdf\n")
}, error = function(e) {
  cat("  WARNING: Plot generation failed:", conditionMessage(e), "\n")
  try(dev.off(), silent = TRUE)
})

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

results_out <- list(
  design           = design,
  results_df       = results_df,
  sensitivity_long = sensitivity_long,
  sensitivity_rank = sensitivity_rank,
  perturbations    = pert,
  config = list(
    chl_levels     = CHL_LEVELS,
    sst            = SST_CONSTANT,
    sim_years      = SIM_YEARS,
    analysis_years = ANALYSIS_YEARS,
    min_R_frac     = MIN_R_FRAC,
    n_cores        = n_cores,
    elapsed_mins   = as.numeric(t_elapsed)
  )
)

saveRDS(results_out, "calibration/oat_sensitivity_results.rds")
cat("  Full results saved to calibration/oat_sensitivity_results.rds\n")

# Save summary CSV
write.csv(sensitivity_rank, "calibration/oat_sensitivity_summary.csv",
          row.names = FALSE)
cat("  Summary saved to calibration/oat_sensitivity_summary.csv\n")

# =============================================================================
# 11. FINAL REPORT
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  SENSITIVITY ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Runs completed: ", n_success, "/", nrow(design), "\n")
cat("Runtime:        ", round(t_elapsed, 1), "minutes\n")
cat("Cores used:     ", n_cores, "\n\n")

cat("Key findings:\n")

# Most influential parameters for total fish biomass
top_bio <- sensitivity_rank %>%
  filter(response == "total_fish_bio") %>%
  head(3)
cat("  Most sensitive (total fish biomass):\n")
for (j in seq_len(nrow(top_bio))) {
  cat(sprintf("    %d. %s (max |%%Δ| = %.1f%%)\n",
              j, top_bio$param[j], top_bio$max_abs_rel_change[j]))
}

# Baseline stability
cat("\n  Baseline stability:\n")
for (i in seq_len(nrow(baselines))) {
  row <- baselines[i, ]
  n_alive  <- row$n_fish_alive
  n_steady <- sum(sapply(paste0("steady_", fish_names),
                         function(x) row[[x]]), na.rm = TRUE)
  cat(sprintf("    Chl=%.2f: %d/%d fish alive, %d/%d steady state\n",
              row$chl, n_alive, length(fish_names),
              n_steady, length(fish_names)))
}

cat("\nNext steps:\n")
cat("  1. Review plots in calibration/oat_sensitivity_plots.pdf\n")
cat("  2. If parameters show large effects, consider wider perturbation ranges\n")
cat("  3. If baseline is unstable, focus on stabilising before differentiation\n")
cat("  4. Load results with: res <- readRDS('calibration/oat_sensitivity_results.rds')\n")

# Clean up parallel workers
plan(sequential)
cat("\nDone.\n")
