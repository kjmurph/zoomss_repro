# =============================================================================
# 02_sensitivity_analysis.R
# One-at-a-time (OAT) sensitivity analysis for calibration parameters
# =============================================================================
#
# Purpose: Vary each candidate parameter while holding all others at defaults.
#          Identifies which parameters most strongly influence the calibration
#          targets (zooplankton proportions, zoo:fish ratio, total biomass).
#          Results inform which parameters to include in the optimisation.
#
# Prerequisites:
#   - Baseline from 01_generate_baseline.R
#   - Revised zoomss package installed (devtools::load_all("."))
#
# Output:
#   - calibration/sensitivity_results.rds   (full OAT results)
#   - calibration/sensitivity_summary.rds   (parameter ranking table)
#   - calibration/sensitivity_plots.pdf     (diagnostic plots)
# =============================================================================

library(zoomss)
library(ggplot2)
library(patchwork)
library(future.apply)

plan(multisession, workers = parallelly::availableCores() - 1)
cat("Using", parallelly::availableCores() - 1, "parallel workers\n")

# ── Load baseline ──
baseline <- readRDS("calibration/baseline_original_zoomss.rds")
chl_levels <- baseline$chl_levels

# ── Define candidate parameters and perturbation ranges ──
# Each parameter is varied across a range while all others are held at defaults.
# The range should span biologically plausible values.

sensitivity_params <- list(
  f_M = list(
    default = 0.50,
    range = seq(0.30, 0.70, by = 0.05),
    description = "Metabolic fraction of assimilated energy"
  ),
  K_growth_zoo_base = list(
    default = 0.35,
    range = seq(0.15, 0.49, by = 0.05),
    description = "Zooplankton growth fraction (base, scaled per group)"
  ),
  K_growth_fish = list(
    default = 0.30,
    range = seq(0.15, 0.45, by = 0.05),
    description = "Fish growth fraction"
  ),
  repro_eff = list(
    default = 0.002,
    range = c(1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2),
    description = "Fish reproductive efficiency (egg-to-recruit survival)"
  ),
  def_low = list(
    default = 0.95,
    range = seq(0.50, 0.95, by = 0.05),
    description = "Defecation fraction for low-Carbon prey"
  ),
  def_high = list(
    default = 0.30,
    range = seq(0.15, 0.45, by = 0.05),
    description = "Defecation fraction for high-Carbon prey"
  )
)

# ── Helper: run model across chl gradient and extract diagnostics ──
run_gradient <- function(Groups, chl_levels, sim_years = 400, avg_years = 100) {
  # Parallelise across chl levels
  results <- future_lapply(chl_levels, function(chl) {
    tryCatch({
      env <- createInputParams(
        time = seq(0, sim_years, by = 0.1),
        sst = 15,
        chl = chl
      )
      mdl <- zoomss_model(input_params = env, Groups = Groups, isave = 10)

      # Average abundance over final avg_years, compute biomass
      avg_N <- averageTimeSeries(mdl, var = "abundance", n_years = avg_years)
      w <- mdl$param$w
      avg_biomass <- sweep(avg_N, 2, w, "*")
      group_biomass <- rowSums(avg_biomass)

      list(group_biomass = group_biomass, success = TRUE)
    }, error = function(e) {
      list(group_biomass = NULL, success = FALSE, error = e$message)
    })
  }, future.seed = TRUE)

  zoo_idx  <- which(Groups$Type == "Zooplankton")
  fish_idx <- which(Groups$Type == "Fish")
  n_chl    <- length(chl_levels)

  zoo_prop <- matrix(NA, n_chl, length(zoo_idx))
  fish_bm  <- matrix(NA, n_chl, length(fish_idx))
  total_bm <- numeric(n_chl)

  for (i in seq_along(chl_levels)) {
    if (!results[[i]]$success) next
    gb <- results[[i]]$group_biomass
    zoo_biomass <- gb[zoo_idx]
    zoo_prop[i, ] <- zoo_biomass / sum(zoo_biomass)
    fish_bm[i, ]  <- gb[fish_idx]
    total_bm[i]   <- sum(gb)
  }

  list(zoo_proportions = zoo_prop, fish_biomass = fish_bm, total_biomass = total_bm)
}

# ── Helper: apply a parameter perturbation to Groups ──
apply_perturbation <- function(Groups_default, param_name, value) {
  Groups <- Groups_default
  zoo_idx  <- which(Groups$Type == "Zooplankton")
  fish_idx <- which(Groups$Type == "Fish")

  switch(param_name,
    f_M = {
      Groups$f_M <- value
    },
    K_growth_zoo_base = {
      default_K_zoo <- Groups$K_growth[zoo_idx]
      K_relative <- default_K_zoo / mean(default_K_zoo)
      Groups$K_growth[zoo_idx] <- pmin(pmax(value * K_relative, 0.05), 0.49)
    },
    K_growth_fish = {
      Groups$K_growth[fish_idx] <- value
    },
    repro_eff = {
      Groups$repro_eff[fish_idx] <- value
    },
    def_low = {
      Groups$def_low <- value
    },
    def_high = {
      Groups$def_high <- value
    }
  )

  # Validate energy budget closure
  R_frac <- 1 - Groups$f_M - Groups$K_growth
  if (any(R_frac < 0)) return(NULL)  # Invalid combination

  Groups
}

# ── Run OAT sensitivity analysis ──
Groups_default <- getGroups()
sensitivity_results <- list()

for (param_name in names(sensitivity_params)) {
  param_info <- sensitivity_params[[param_name]]
  cat("\n=== Sensitivity analysis for:", param_name, "===\n")
  cat("   ", param_info$description, "\n")
  cat("    Range:", paste(param_info$range, collapse = ", "), "\n")

  param_results <- list()
  for (val in param_info$range) {
    cat("  Running", param_name, "=", val, "...")

    Groups_mod <- apply_perturbation(Groups_default, param_name, val)
    if (is.null(Groups_mod)) {
      cat(" SKIPPED (invalid energy budget)\n")
      next
    }

    tryCatch({
      res <- run_gradient(Groups_mod, chl_levels)
      res$param_value <- val
      param_results <- c(param_results, list(res))
      cat(" done\n")
    }, error = function(e) {
      cat(" ERROR:", e$message, "\n")
    })
  }

  sensitivity_results[[param_name]] <- param_results
}

saveRDS(sensitivity_results, "calibration/sensitivity_results.rds")
cat("\nSaved full sensitivity results to calibration/sensitivity_results.rds\n")

# ── Compute sensitivity metrics ──
# For each parameter, compute the range of effect on key outputs:
# 1. Mean absolute change in zoo proportions across gradient
# 2. Mean absolute change in log10 total biomass
# 3. Mean absolute change in log10 zoo:fish ratio

sensitivity_summary <- data.frame(
  parameter = character(),
  delta_zoo_prop = numeric(),   # Sensitivity of zoo proportions
  delta_total_bm = numeric(),   # Sensitivity of total biomass
  delta_zf_ratio = numeric(),   # Sensitivity of zoo:fish ratio
  stringsAsFactors = FALSE
)

for (param_name in names(sensitivity_results)) {
  results_list <- sensitivity_results[[param_name]]
  if (length(results_list) < 2) next

  # Compare each perturbation to the default run
  default_idx <- which(sapply(results_list, function(x) {
    abs(x$param_value - sensitivity_params[[param_name]]$default) < 1e-8
  }))

  if (length(default_idx) == 0) {
    # Use midpoint as reference
    default_idx <- ceiling(length(results_list) / 2)
  }

  ref <- results_list[[default_idx]]
  deltas_prop  <- numeric(length(results_list))
  deltas_bm    <- numeric(length(results_list))
  deltas_ratio <- numeric(length(results_list))

  for (j in seq_along(results_list)) {
    res <- results_list[[j]]
    deltas_prop[j] <- mean(abs(res$zoo_proportions - ref$zoo_proportions), na.rm = TRUE)
    deltas_bm[j]   <- mean(abs(log10(pmax(res$total_biomass, 1e-20)) -
                                log10(pmax(ref$total_biomass, 1e-20))), na.rm = TRUE)

    # Zoo:fish ratio from proportions and biomass
    ref_zoo_bm <- rowSums(ref$zoo_proportions * ref$total_biomass)
    ref_fish_bm <- rowSums(ref$fish_biomass)
    res_zoo_bm <- rowSums(res$zoo_proportions * res$total_biomass)
    res_fish_bm <- rowSums(res$fish_biomass)

    ref_ratio <- ref_zoo_bm / pmax(ref_fish_bm, 1e-20)
    res_ratio <- res_zoo_bm / pmax(res_fish_bm, 1e-20)
    deltas_ratio[j] <- mean(abs(log10(pmax(res_ratio, 1e-20)) -
                                 log10(pmax(ref_ratio, 1e-20))), na.rm = TRUE)
  }

  sensitivity_summary <- rbind(sensitivity_summary, data.frame(
    parameter = param_name,
    delta_zoo_prop = max(deltas_prop),
    delta_total_bm = max(deltas_bm),
    delta_zf_ratio = max(deltas_ratio)
  ))
}

cat("\n=== SENSITIVITY SUMMARY ===\n")
cat("(larger values = higher sensitivity = better calibration parameter)\n\n")
print(sensitivity_summary[order(-sensitivity_summary$delta_zoo_prop), ])

# ── Save summary ──
saveRDS(sensitivity_summary, "calibration/sensitivity_summary.rds")

# ── Generate diagnostic plots ──
cat("\nGenerating sensitivity plots...\n")

plot_list <- list()
for (param_name in names(sensitivity_results)) {
  results_list <- sensitivity_results[[param_name]]
  if (length(results_list) < 2) next

  # Build data frame for plotting
  plot_df <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      param_value = res$param_value,
      total_biomass = res$total_biomass,
      chl = chl_levels
    )
  }))

  p <- ggplot(plot_df, aes(x = factor(chl), y = total_biomass,
                            colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 1) +
    scale_colour_viridis_c() +
    labs(
      title = param_name,
      subtitle = sensitivity_params[[param_name]]$description,
      x = "Chlorophyll (mg/m³)",
      y = "Total Biomass",
      colour = param_name
    ) +
    theme_minimal()

  plot_list[[param_name]] <- p
}

if (length(plot_list) > 0) {
  combined <- wrap_plots(plot_list, ncol = 2)
  ggsave("calibration/sensitivity_plots.pdf", combined,
         width = 14, height = 10)
  cat("Saved plots to calibration/sensitivity_plots.pdf\n")
}

cat("\nSensitivity analysis complete.\n")
cat("Use the summary table above to decide which parameters to include in optimisation.\n")
cat("Parameters with large delta_zoo_prop values are strong candidates.\n")
