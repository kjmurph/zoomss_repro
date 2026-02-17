# =============================================================================
# 05_evaluate_calibration.R
# Evaluate and visualise calibration results
# =============================================================================
#
# Purpose: Compare calibrated revised ZooMSS against the original baseline.
#          Generates diagnostic plots showing how well the calibration matched
#          zooplankton community composition across the productivity gradient.
#
# Prerequisites:
#   - Calibration results from 04_run_calibration.R
#   - Revised zoomss package installed (devtools::load_all("."))
#
# Output:
#   - calibration/calibration_diagnostics.pdf
#   - calibration/calibrated_groups.rds  (Groups with best parameters applied)
# =============================================================================

library(zoomss)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(future.apply)

plan(multisession, workers = parallelly::availableCores() - 1)
cat("Using", parallelly::availableCores() - 1, "parallel workers\n")

# ── Load results ──
results <- readRDS("calibration/calibration_results.rds")
baseline <- results$baseline
chl_levels <- baseline$chl_levels

# ── Pick best result (compare across methods) ──
best_methods <- results$comparison
cat("Method comparison:\n")
print(best_methods)

best_row <- which.min(best_methods$objective)
best_method <- best_methods$method[best_row]
cat("\nBest method:", best_method, "(objective =", best_methods$objective[best_row], ")\n")

# Extract best parameters
if (best_method == "DEoptim" && !is.null(results$deoptim)) {
  best_par <- results$deoptim$optim$bestmem
} else if (best_method == "Nelder-Mead") {
  best_par <- results$nelder_mead$par
} else {
  best_par <- results$lbfgsb$par
}

cat("\nBest parameters:\n")
print(best_par)

# ── Apply best parameters to Groups ──
Groups <- getGroups()
zoo_idx  <- which(Groups$Type == "Zooplankton")
fish_idx <- which(Groups$Type == "Fish")

Groups$f_M <- best_par["f_M"]

K_base_zoo <- best_par["K_growth_zoo_base"]
default_K_zoo <- Groups$K_growth[zoo_idx]
K_relative <- default_K_zoo / mean(default_K_zoo)
Groups$K_growth[zoo_idx] <- pmin(pmax(K_base_zoo * K_relative, 0.05), 0.49)
Groups$K_growth[fish_idx] <- best_par["K_growth_fish"]
Groups$repro_eff[fish_idx] <- best_par["repro_eff"]

# Print the calibrated energy budget
cat("\nCalibrated energy budget:\n")
R_frac <- 1 - Groups$f_M - Groups$K_growth
print(data.frame(
  Species   = Groups$Species,
  Type      = Groups$Type,
  f_M       = Groups$f_M,
  K_growth  = round(Groups$K_growth, 4),
  R_frac    = round(R_frac, 4),
  repro_eff = Groups$repro_eff
))

# Save calibrated Groups for future use
saveRDS(Groups, "calibration/calibrated_groups.rds")

# ── Run calibrated model across gradient (parallelised) ──
cat("\nRunning calibrated model across chlorophyll gradient...\n")

revised_results <- future_lapply(chl_levels, function(chl) {
  env <- createInputParams(
    time = seq(0, 400, by = 0.1),
    sst  = 15,
    chl  = chl
  )
  mdl <- zoomss_model(input_params = env, Groups = Groups, isave = 10)

  # Extract steady-state diagnostics using getBiomass
  Biomass <- getBiomass(mdl, units = "ww")
  time_vec <- mdl$time
  time_idx <- which(time_vec >= max(time_vec) - 100)
  avg_biomass <- apply(Biomass[time_idx, , , drop = FALSE], c(2, 3), mean)
  group_biomass <- rowSums(avg_biomass)

  # Fish reproduction metrics (final 100 years)
  list(
    group_biomass = group_biomass,
    avg_biomass   = avg_biomass,
    mdl           = mdl  # Keep for size spectra
  )
}, future.seed = TRUE)

# ── Extract revised diagnostics ──
n_chl  <- length(chl_levels)
n_zoo  <- length(zoo_idx)
n_fish <- length(fish_idx)

rev_zoo_prop   <- matrix(NA, n_chl, n_zoo)
rev_fish_bm    <- matrix(NA, n_chl, n_fish)
rev_total_bm   <- numeric(n_chl)
rev_zf_ratio   <- numeric(n_chl)

for (i in seq_along(chl_levels)) {
  gb <- revised_results[[i]]$group_biomass
  zoo_bm  <- gb[zoo_idx]
  fish_bm <- gb[fish_idx]

  rev_zoo_prop[i, ] <- zoo_bm / sum(zoo_bm)
  rev_fish_bm[i, ]  <- fish_bm
  rev_total_bm[i]   <- sum(gb)
  rev_zf_ratio[i]   <- sum(zoo_bm) / sum(fish_bm)
}

colnames(rev_zoo_prop)  <- baseline$zoo_species
colnames(rev_fish_bm)   <- baseline$fish_species

# =============================================================================
# DIAGNOSTIC PLOTS
# =============================================================================

cat("Generating diagnostic plots...\n")

# ── Plot 1: Zooplankton proportions — baseline vs calibrated ──
zoo_df <- rbind(
  data.frame(
    chl     = rep(chl_levels, n_zoo),
    species = rep(baseline$zoo_species, each = n_chl),
    prop    = as.vector(baseline$zoo_proportions),
    model   = "Original"
  ),
  data.frame(
    chl     = rep(chl_levels, n_zoo),
    species = rep(baseline$zoo_species, each = n_chl),
    prop    = as.vector(rev_zoo_prop),
    model   = "Calibrated"
  )
)

p1 <- ggplot(zoo_df, aes(x = factor(chl), y = prop, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~model) +
  labs(
    title = "Zooplankton Community Composition",
    subtitle = "Stacked proportions across chlorophyll gradient",
    x     = "Chlorophyll (mg/m³)",
    y     = "Proportion of zooplankton biomass",
    fill  = "Species"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ── Plot 2: Zoo:Fish ratio — baseline vs calibrated ──
ratio_df <- data.frame(
  chl   = rep(chl_levels, 2),
  ratio = c(baseline$zoo_fish_ratio, rev_zf_ratio),
  model = rep(c("Original", "Calibrated"), each = n_chl)
)

p2 <- ggplot(ratio_df, aes(x = chl, y = ratio, colour = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title  = "Zooplankton:Fish Biomass Ratio",
    x      = "Chlorophyll (mg/m³)",
    y      = "Zoo:Fish ratio",
    colour = "Model"
  ) +
  theme_minimal()

# ── Plot 3: Total biomass — baseline vs calibrated ──
bm_df <- data.frame(
  chl    = rep(chl_levels, 2),
  biomass = c(baseline$total_biomass, rev_total_bm),
  model  = rep(c("Original", "Calibrated"), each = n_chl)
)

p3 <- ggplot(bm_df, aes(x = chl, y = biomass, colour = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title  = "Total Community Biomass",
    x      = "Chlorophyll (mg/m³)",
    y      = "Total Biomass (g ww)",
    colour = "Model"
  ) +
  theme_minimal()

# ── Plot 4: Size spectra comparison at 3 representative chl levels ──
# Show averaged spectra at low (0.1), medium (1.0), and high (10.0) chl
chl_plot_idx <- which(chl_levels %in% c(0.1, 1.0, 10.0))

spectra_df <- data.frame()
for (i in chl_plot_idx) {
  mdl <- revised_results[[i]]$mdl
  w   <- mdl$param$w

  avg_N <- averageTimeSeries(mdl, var = "abundance", n_years = 100)
  # Sum across groups for community spectrum
  community_N <- colSums(avg_N)
  community_bm <- community_N * w

  spectra_df <- rbind(spectra_df, data.frame(
    log10_size = log10(w),
    log10_bm   = log10(pmax(community_bm, 1e-30)),
    chl        = paste0("Chl = ", chl_levels[i], " mg/m³")
  ))
}

p4 <- ggplot(spectra_df, aes(x = log10_size, y = log10_bm, colour = chl)) +
  geom_line(linewidth = 0.8) +
  labs(
    title  = "Community Size Spectra (Calibrated Model)",
    x      = expression(log[10]~"Body Size (g)"),
    y      = expression(log[10]~"Biomass"),
    colour = NULL
  ) +
  theme_minimal()

# ── Plot 5: Fish SSB and recruitment (new model only) ──
fish_species <- baseline$fish_species

ssb_df <- data.frame()
rec_df <- data.frame()

for (i in seq_along(chl_levels)) {
  mdl <- revised_results[[i]]$mdl
  if (!is.null(mdl$SSB)) {
    # Average final 100 years of SSB/recruitment
    nsave <- dim(mdl$SSB)[1]
    dt_saved <- mdl$param$isave * mdl$param$dt
    n_steps <- min(nsave, round(100 / dt_saved))
    start_idx <- max(1, nsave - n_steps + 1)

    avg_ssb <- colMeans(mdl$SSB[start_idx:nsave, , drop = FALSE])
    avg_rec <- colMeans(mdl$recruitment[start_idx:nsave, , drop = FALSE])

    ssb_df <- rbind(ssb_df, data.frame(
      chl     = chl_levels[i],
      species = fish_species,
      ssb     = avg_ssb
    ))
    rec_df <- rbind(rec_df, data.frame(
      chl     = chl_levels[i],
      species = fish_species,
      recruitment = avg_rec
    ))
  }
}

p5 <- ggplot(ssb_df, aes(x = chl, y = ssb, colour = species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title  = "Spawning Stock Biomass (Calibrated Model)",
    x      = "Chlorophyll (mg/m³)",
    y      = "SSB",
    colour = "Species"
  ) +
  theme_minimal()

p6 <- ggplot(rec_df, aes(x = chl, y = recruitment, colour = species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title  = "Recruitment (Calibrated Model)",
    x      = "Chlorophyll (mg/m³)",
    y      = "Recruitment",
    colour = "Species"
  ) +
  theme_minimal()

# ── Combine and save ──
combined <- (p1 / (p2 | p3)) / (p4 | (p5 / p6)) +
  plot_annotation(
    title    = "ZooMSS Calibration Diagnostics",
    subtitle = paste("Best method:", best_method,
                     "| Objective:", round(best_methods$objective[best_row], 4))
  )

ggsave("calibration/calibration_diagnostics.pdf", combined,
       width = 16, height = 18)
cat("Saved diagnostic plots to calibration/calibration_diagnostics.pdf\n")

# ── Print summary statistics ──
cat("\n=== CALIBRATION EVALUATION ===\n\n")

cat("Zooplankton proportion RMSE:\n")
prop_rmse <- sqrt(mean((rev_zoo_prop - baseline$zoo_proportions)^2))
cat("  Overall:", round(prop_rmse, 4), "\n")

for (j in seq_along(baseline$zoo_species)) {
  sp_rmse <- sqrt(mean((rev_zoo_prop[, j] - baseline$zoo_proportions[, j])^2))
  cat("  ", baseline$zoo_species[j], ":", round(sp_rmse, 4), "\n")
}

cat("\nZoo:Fish ratio comparison:\n")
print(data.frame(
  chl       = chl_levels,
  original  = round(baseline$zoo_fish_ratio, 2),
  calibrated = round(rev_zf_ratio, 2),
  ratio_diff = round(rev_zf_ratio - baseline$zoo_fish_ratio, 2)
))

cat("\nTotal biomass comparison:\n")
print(data.frame(
  chl            = chl_levels,
  original       = round(baseline$total_biomass, 4),
  calibrated     = round(rev_total_bm, 4),
  log10_diff     = round(log10(rev_total_bm) - log10(baseline$total_biomass), 3)
))

cat("\nCalibration evaluation complete.\n")
