# =============================================================================
# 03_sensitivity_plots.R
# Extended sensitivity diagnostic plots with common y-axes for comparison
# =============================================================================
#
# Purpose: Generate three sets of diagnostic plots from the OAT sensitivity
#          results, each using a common y-axis scale across all parameter panels
#          for consistent visual comparison.
#
# Outputs:
#   1. sensitivity_fish_biomass.pdf   – Total fish biomass vs chl, per parameter
#   2. sensitivity_zoo_biomass.pdf    – Total zoo biomass vs chl, per parameter
#   3. sensitivity_zoo_by_group.pdf   – Per-group zoo biomass vs chl (faceted)
#   4. sensitivity_total_biomass.pdf  – Total biomass vs chl (common y-axis redo)
#
# Prerequisites:
#   - calibration/sensitivity_results.rds from 02_sensitivity_analysis.R
# =============================================================================

library(ggplot2)
library(patchwork)
library(zoomss)

# ── Load results and metadata ──
sensitivity_results <- readRDS("calibration/sensitivity_results.rds")

# Parameter descriptions for subtitles
param_descriptions <- c(
  f_M               = "Metabolic fraction of assimilated energy (fish only)",
  K_growth_zoo_base = "Zooplankton growth fraction (base, scaled per group)",
  K_growth_fish     = "Fish growth fraction",
  repro_eff         = "Fish reproductive efficiency (egg-to-recruit survival)",
  def_low           = "Defecation fraction for low-Carbon prey",
  def_high          = "Defecation fraction for high-Carbon prey"
)

# Group metadata from the package
Groups <- getGroups()
zoo_idx  <- which(Groups$Type == "Zooplankton")
fish_idx <- which(Groups$Type == "Fish")
zoo_names  <- Groups$Species[zoo_idx]
fish_names <- Groups$Species[fish_idx]
zoo_colours  <- Groups$PlotColour[zoo_idx]
fish_colours <- Groups$PlotColour[fish_idx]
names(zoo_colours)  <- zoo_names
names(fish_colours) <- fish_names

# Chl levels
baseline <- readRDS("calibration/baseline_original_zoomss.rds")
chl_levels <- baseline$chl_levels

# ── Common theme ──
theme_sens <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "right"
  )

# =============================================================================
# Helper: build long-format data frames from results
# =============================================================================

build_fish_df <- function(sensitivity_results, chl_levels) {
  rows <- list()
  for (param_name in names(sensitivity_results)) {
    for (res in sensitivity_results[[param_name]]) {
      total_fish <- rowSums(res$fish_biomass)
      rows <- c(rows, list(data.frame(
        parameter   = param_name,
        param_value = res$param_value,
        chl         = chl_levels,
        fish_biomass = total_fish,
        stringsAsFactors = FALSE
      )))
    }
  }
  do.call(rbind, rows)
}

build_zoo_df <- function(sensitivity_results, chl_levels) {
  rows <- list()
  for (param_name in names(sensitivity_results)) {
    for (res in sensitivity_results[[param_name]]) {
      total_zoo <- res$total_biomass - rowSums(res$fish_biomass)
      rows <- c(rows, list(data.frame(
        parameter    = param_name,
        param_value  = res$param_value,
        chl          = chl_levels,
        zoo_biomass  = total_zoo,
        stringsAsFactors = FALSE
      )))
    }
  }
  do.call(rbind, rows)
}

build_total_df <- function(sensitivity_results, chl_levels) {
  rows <- list()
  for (param_name in names(sensitivity_results)) {
    for (res in sensitivity_results[[param_name]]) {
      rows <- c(rows, list(data.frame(
        parameter      = param_name,
        param_value    = res$param_value,
        chl            = chl_levels,
        total_biomass  = res$total_biomass,
        stringsAsFactors = FALSE
      )))
    }
  }
  do.call(rbind, rows)
}

build_zoo_group_df <- function(sensitivity_results, chl_levels, zoo_names) {
  rows <- list()
  for (param_name in names(sensitivity_results)) {
    for (res in sensitivity_results[[param_name]]) {
      # Reconstruct per-group absolute zoo biomass
      total_zoo <- res$total_biomass - rowSums(res$fish_biomass)
      zoo_abs <- sweep(res$zoo_proportions, 1, total_zoo, "*")
      for (j in seq_along(zoo_names)) {
        rows <- c(rows, list(data.frame(
          parameter   = param_name,
          param_value = res$param_value,
          chl         = chl_levels,
          group       = zoo_names[j],
          biomass     = zoo_abs[, j],
          stringsAsFactors = FALSE
        )))
      }
    }
  }
  df <- do.call(rbind, rows)
  df$group <- factor(df$group, levels = zoo_names)
  df
}

build_fish_group_df <- function(sensitivity_results, chl_levels, fish_names) {
  rows <- list()
  for (param_name in names(sensitivity_results)) {
    for (res in sensitivity_results[[param_name]]) {
      for (j in seq_along(fish_names)) {
        rows <- c(rows, list(data.frame(
          parameter   = param_name,
          param_value = res$param_value,
          chl         = chl_levels,
          group       = fish_names[j],
          biomass     = res$fish_biomass[, j],
          stringsAsFactors = FALSE
        )))
      }
    }
  }
  df <- do.call(rbind, rows)
  df$group <- factor(df$group, levels = fish_names)
  df
}

# ── Build data frames ──
cat("Building data frames...\n")
fish_df       <- build_fish_df(sensitivity_results, chl_levels)
zoo_df        <- build_zoo_df(sensitivity_results, chl_levels)
total_df      <- build_total_df(sensitivity_results, chl_levels)
zoo_group_df  <- build_zoo_group_df(sensitivity_results, chl_levels, zoo_names)
fish_group_df <- build_fish_group_df(sensitivity_results, chl_levels, fish_names)

# Parameter display order (by sensitivity ranking)
param_order <- c("f_M", "K_growth_zoo_base", "def_high",
                 "def_low", "repro_eff", "K_growth_fish")
fish_df$parameter       <- factor(fish_df$parameter, levels = param_order)
zoo_df$parameter        <- factor(zoo_df$parameter, levels = param_order)
total_df$parameter      <- factor(total_df$parameter, levels = param_order)
zoo_group_df$parameter  <- factor(zoo_group_df$parameter, levels = param_order)
fish_group_df$parameter <- factor(fish_group_df$parameter, levels = param_order)

# =============================================================================
# 1. Total biomass vs chl (common y-axis) – redo of original
# =============================================================================
cat("Generating total biomass plots...\n")

y_range_total <- range(total_df$total_biomass, na.rm = TRUE)

plot_list_total <- list()
for (p in param_order) {
  sub <- total_df[total_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_total[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = total_biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 1.2) +
    scale_colour_viridis_c() +
    scale_y_continuous(limits = y_range_total) +
    labs(
      title    = p,
      subtitle = param_descriptions[p],
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Total Biomass (g ww/m\u00B2)",
      colour = p
    ) +
    theme_sens
}

if (length(plot_list_total) > 0) {
  combined <- wrap_plots(plot_list_total, ncol = 2)
  ggsave("calibration/sensitivity_total_biomass.pdf", combined,
         width = 14, height = 12)
  cat("  Saved calibration/sensitivity_total_biomass.pdf\n")
}

# =============================================================================
# 2. Fish biomass vs chl (common y-axis)
# =============================================================================
cat("Generating fish biomass plots...\n")

y_range_fish <- range(fish_df$fish_biomass, na.rm = TRUE)

plot_list_fish <- list()
for (p in param_order) {
  sub <- fish_df[fish_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_fish[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = fish_biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 1.2) +
    scale_colour_viridis_c() +
    scale_y_continuous(limits = y_range_fish) +
    labs(
      title    = p,
      subtitle = param_descriptions[p],
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Total Fish Biomass (g ww/m\u00B2)",
      colour = p
    ) +
    theme_sens
}

if (length(plot_list_fish) > 0) {
  combined <- wrap_plots(plot_list_fish, ncol = 2)
  ggsave("calibration/sensitivity_fish_biomass.pdf", combined,
         width = 14, height = 12)
  cat("  Saved calibration/sensitivity_fish_biomass.pdf\n")
}

# =============================================================================
# 3. Zoo biomass vs chl (common y-axis)
# =============================================================================
cat("Generating zooplankton biomass plots...\n")

y_range_zoo <- range(zoo_df$zoo_biomass, na.rm = TRUE)

plot_list_zoo <- list()
for (p in param_order) {
  sub <- zoo_df[zoo_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_zoo[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = zoo_biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 1.2) +
    scale_colour_viridis_c() +
    scale_y_continuous(limits = y_range_zoo) +
    labs(
      title    = p,
      subtitle = param_descriptions[p],
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Total Zooplankton Biomass (g ww/m\u00B2)",
      colour = p
    ) +
    theme_sens
}

if (length(plot_list_zoo) > 0) {
  combined <- wrap_plots(plot_list_zoo, ncol = 2)
  ggsave("calibration/sensitivity_zoo_biomass.pdf", combined,
         width = 14, height = 12)
  cat("  Saved calibration/sensitivity_zoo_biomass.pdf\n")
}

# =============================================================================
# 4. Faceted per-group zooplankton biomass (common y-axis within each param)
# =============================================================================
cat("Generating faceted zooplankton per-group biomass plots...\n")

# Compute global y range across all zoo groups for consistent comparison
y_range_zoo_grp <- range(zoo_group_df$biomass, na.rm = TRUE)

plot_list_zoo_grp <- list()
for (p in param_order) {
  sub <- zoo_group_df[zoo_group_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_zoo_grp[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 0.8) +
    scale_colour_viridis_c() +
    facet_wrap(~ group, scales = "free_y", ncol = 3) +
    labs(
      title    = paste0("Zooplankton per-group biomass — ", p),
      subtitle = param_descriptions[p],
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Biomass (g ww/m\u00B2)",
      colour = p
    ) +
    theme_sens +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1)
    )
}

if (length(plot_list_zoo_grp) > 0) {
  # One page per parameter (9 facets each)
  pdf("calibration/sensitivity_zoo_by_group.pdf", width = 14, height = 12)
  for (p in names(plot_list_zoo_grp)) {
    print(plot_list_zoo_grp[[p]])
  }
  dev.off()
  cat("  Saved calibration/sensitivity_zoo_by_group.pdf\n")
}

# =============================================================================
# 5. Faceted per-group fish biomass (common y-axis within each param)
# =============================================================================
cat("Generating faceted fish per-group biomass plots...\n")

y_range_fish_grp <- range(fish_group_df$biomass, na.rm = TRUE)

plot_list_fish_grp <- list()
for (p in param_order) {
  sub <- fish_group_df[fish_group_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_fish_grp[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 1) +
    scale_colour_viridis_c() +
    facet_wrap(~ group, scales = "free_y", ncol = 3) +
    labs(
      title    = paste0("Fish per-group biomass — ", p),
      subtitle = param_descriptions[p],
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Biomass (g ww/m\u00B2)",
      colour = p
    ) +
    theme_sens +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1)
    )
}

if (length(plot_list_fish_grp) > 0) {
  pdf("calibration/sensitivity_fish_by_group.pdf", width = 14, height = 6)
  for (p in names(plot_list_fish_grp)) {
    print(plot_list_fish_grp[[p]])
  }
  dev.off()
  cat("  Saved calibration/sensitivity_fish_by_group.pdf\n")
}

# =============================================================================
# 6. Faceted per-group zooplankton biomass — COMMON y-axis across all groups
# =============================================================================
cat("Generating faceted zooplankton per-group biomass (common y-axis)...\n")

plot_list_zoo_grp_cy <- list()
for (p in param_order) {
  sub <- zoo_group_df[zoo_group_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_zoo_grp_cy[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 0.8) +
    scale_colour_viridis_c() +
    facet_wrap(~ group, ncol = 3) +
    labs(
      title    = paste0("Zooplankton per-group biomass (common y) \u2014 ", p),
      subtitle = param_descriptions[p],
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Biomass (g ww/m\u00B2)",
      colour = p
    ) +
    theme_sens +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1)
    )
}

if (length(plot_list_zoo_grp_cy) > 0) {
  pdf("calibration/sensitivity_zoo_by_group_common_y.pdf", width = 14, height = 12)
  for (p in names(plot_list_zoo_grp_cy)) {
    print(plot_list_zoo_grp_cy[[p]])
  }
  dev.off()
  cat("  Saved calibration/sensitivity_zoo_by_group_common_y.pdf\n")
}

# =============================================================================
# 7. Faceted per-group fish biomass — COMMON y-axis across all groups
# =============================================================================
cat("Generating faceted fish per-group biomass (common y-axis)...\n")

plot_list_fish_grp_cy <- list()
for (p in param_order) {
  sub <- fish_group_df[fish_group_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_fish_grp_cy[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 1) +
    scale_colour_viridis_c() +
    facet_wrap(~ group, ncol = 3) +
    labs(
      title    = paste0("Fish per-group biomass (common y) \u2014 ", p),
      subtitle = param_descriptions[p],
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Biomass (g ww/m\u00B2)",
      colour = p
    ) +
    theme_sens +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1)
    )
}

if (length(plot_list_fish_grp_cy) > 0) {
  pdf("calibration/sensitivity_fish_by_group_common_y.pdf", width = 14, height = 6)
  for (p in names(plot_list_fish_grp_cy)) {
    print(plot_list_fish_grp_cy[[p]])
  }
  dev.off()
  cat("  Saved calibration/sensitivity_fish_by_group_common_y.pdf\n")
}

# =============================================================================
# 8. Faceted per-group zooplankton — NORMALISED biomass
#    Normalised = biomass / baseline_biomass for each group × chl combination,
#    where baseline is the run with the default parameter value.
# =============================================================================
cat("Generating faceted zooplankton per-group normalised biomass...\n")

# Identify default parameter values (same as in 02_sensitivity_analysis.R)
param_defaults <- c(
  f_M               = 0.50,
  K_growth_zoo_base = 0.411,
  K_growth_fish     = 0.30,
  repro_eff         = 0.001,
  def_low           = 0.95,
  def_high          = 0.30
)

# Build normalised zoo group data frame
zoo_group_norm_df <- zoo_group_df  # copy
zoo_group_norm_df$norm_biomass <- NA_real_

for (p in param_order) {
  default_val <- param_defaults[p]
  for (grp in zoo_names) {
    # Baseline biomass for this group at each chl level
    mask_base <- zoo_group_df$parameter == p &
                 zoo_group_df$group == grp &
                 abs(zoo_group_df$param_value - default_val) < 1e-8
    base_bm <- zoo_group_df$biomass[mask_base]
    if (length(base_bm) != length(chl_levels)) next

    # Normalise all runs for this param × group
    mask_all <- zoo_group_norm_df$parameter == p &
                zoo_group_norm_df$group == grp
    idx <- which(mask_all)
    # Match chl levels to get corresponding baseline
    chl_match <- match(zoo_group_norm_df$chl[idx], chl_levels)
    zoo_group_norm_df$norm_biomass[idx] <-
      zoo_group_norm_df$biomass[idx] / base_bm[chl_match]
  }
}

# Drop rows where normalisation failed (e.g. baseline biomass = 0)
zoo_group_norm_df <- zoo_group_norm_df[
  is.finite(zoo_group_norm_df$norm_biomass), ]

plot_list_zoo_norm <- list()
for (p in param_order) {
  sub <- zoo_group_norm_df[zoo_group_norm_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_zoo_norm[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = norm_biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_colour_viridis_c() +
    facet_wrap(~ group, ncol = 3) +
    labs(
      title    = paste0("Zooplankton normalised biomass \u2014 ", p),
      subtitle = paste0(param_descriptions[p],
                        " (1.0 = default: ", param_defaults[p], ")"),
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Normalised Biomass (fold change from default)",
      colour = p
    ) +
    theme_sens +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1)
    )
}

if (length(plot_list_zoo_norm) > 0) {
  pdf("calibration/sensitivity_zoo_by_group_normalised.pdf",
      width = 14, height = 12)
  for (p in names(plot_list_zoo_norm)) {
    print(plot_list_zoo_norm[[p]])
  }
  dev.off()
  cat("  Saved calibration/sensitivity_zoo_by_group_normalised.pdf\n")
}

# =============================================================================
# 9. Faceted per-group fish — NORMALISED biomass
# =============================================================================
cat("Generating faceted fish per-group normalised biomass...\n")

fish_group_norm_df <- fish_group_df
fish_group_norm_df$norm_biomass <- NA_real_

for (p in param_order) {
  default_val <- param_defaults[p]
  for (grp in fish_names) {
    mask_base <- fish_group_df$parameter == p &
                 fish_group_df$group == grp &
                 abs(fish_group_df$param_value - default_val) < 1e-8
    base_bm <- fish_group_df$biomass[mask_base]
    if (length(base_bm) != length(chl_levels)) next

    mask_all <- fish_group_norm_df$parameter == p &
                fish_group_norm_df$group == grp
    idx <- which(mask_all)
    chl_match <- match(fish_group_norm_df$chl[idx], chl_levels)
    fish_group_norm_df$norm_biomass[idx] <-
      fish_group_norm_df$biomass[idx] / base_bm[chl_match]
  }
}

fish_group_norm_df <- fish_group_norm_df[
  is.finite(fish_group_norm_df$norm_biomass), ]

plot_list_fish_norm <- list()
for (p in param_order) {
  sub <- fish_group_norm_df[fish_group_norm_df$parameter == p, ]
  if (nrow(sub) == 0) next
  plot_list_fish_norm[[p]] <- ggplot(sub,
      aes(x = factor(chl), y = norm_biomass,
          colour = param_value, group = param_value)) +
    geom_line() +
    geom_point(size = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_colour_viridis_c() +
    facet_wrap(~ group, ncol = 3) +
    labs(
      title    = paste0("Fish normalised biomass \u2014 ", p),
      subtitle = paste0(param_descriptions[p],
                        " (1.0 = default: ", param_defaults[p], ")"),
      x = "Chlorophyll (mg/m\u00B3)",
      y = "Normalised Biomass (fold change from default)",
      colour = p
    ) +
    theme_sens +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1)
    )
}

if (length(plot_list_fish_norm) > 0) {
  pdf("calibration/sensitivity_fish_by_group_normalised.pdf",
      width = 14, height = 6)
  for (p in names(plot_list_fish_norm)) {
    print(plot_list_fish_norm[[p]])
  }
  dev.off()
  cat("  Saved calibration/sensitivity_fish_by_group_normalised.pdf\n")
}

cat("\nAll sensitivity plots generated.\n")
