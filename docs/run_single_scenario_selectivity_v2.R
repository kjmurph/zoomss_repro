# =============================================================================
# run_single_scenario.R
#
# Quick interactive sweep: one repro_eff scenario across an effort gradient.
# Plots yield curves, biomass depletion, B/B0, selectivity, and size spectra.
#
# Now includes configurable minimum catch size (w_min_catch) per fish group,
# enforcing realistic fishing selectivity (e.g., no fish below 10g caught).
#
# Uses built-in ZooMSS plotting functions (plotSizeSpectra, plotTimeSeries,
# plotReproduction, plotGrowthComparison) for model diagnostics.
#
# Usage: source from an interactive R session within the zoomss project.
# =============================================================================

library(here)
devtools::load_all(here())
library(parallel)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# =============================================================================
# CONFIGURE YOUR SCENARIO HERE
# =============================================================================

# Reproductive efficiency for each fish group: c(Small, Med, Large)
# repro_eff_vec <- c(1e-5, 1e-5, 1e-6)
repro_eff_vec <- c(1, 1, 1)

# Effort gradient
effort_levels <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                   0.5, 1.0, 1.5, 2.0, 3.0)

# Catchability (q = 1 means F = effort × selectivity)
q_fixed <- 1.0

# ---------------------------------------------------------------------------
# MINIMUM CATCH SIZE (grams) — per fish group: c(Small, Med, Large)
# ---------------------------------------------------------------------------
# This sets the knife-edge selectivity threshold: fish below this mass are
# not caught. Values are in grams (converted to log10 internally).
#
# Rationale:
#   - 10g is a realistic minimum for commercial fishing gear
#   - Protects juvenile/immature size classes from fishing mortality
#   - For large fish, a higher threshold (e.g., near maturity) shelters the
#     sub-adult size classes that sustain recruitment
#
# Set to NA to use the model's default (group W0 — minimum body size)
#
# Suggested configurations:
#   c(10, 10, 10)       — Uniform 10g minimum for all groups
#   c(10, 10, 1000)     — Protect large fish juveniles up to 1kg
#   c(10, 10, 10000)    — Only catch mature large fish (Wmat ~ 10^4 g)
#   c(NA, NA, NA)       — Model defaults (catches from group W0)
# ---------------------------------------------------------------------------
w_min_catch_g <- c(10, 10, 10)

# Simulation settings
sim_time  <- seq(0, 400, by = 0.1)
sst_const <- 15
chl_const <- 1.0
isave     <- 5

# Parallel settings
n_cores <- max(1, detectCores() - 1)

# =============================================================================
# SETUP
# =============================================================================

Groups_base <- getGroups()
fish_rows   <- which(Groups_base$Type == "Fish")
fish_names  <- c("Fish_Small", "Fish_Med", "Fish_Large")
pkg_path    <- here()

# Store original Fmort_W0 for comparison
original_Fmort_W0 <- Groups_base$Fmort_W0[fish_rows]
original_W0       <- Groups_base$W0[fish_rows]
original_Wmat     <- Groups_base$Wmat[fish_rows]

# Convert w_min_catch from grams to log10
w_min_catch_log10 <- rep(NA, length(fish_rows))
for (i in seq_along(w_min_catch_g)) {
  if (!is.na(w_min_catch_g[i]) && w_min_catch_g[i] > 0) {
    w_min_catch_log10[i] <- log10(w_min_catch_g[i])
  }
}

# Report selectivity settings
cat("=== Fishing Selectivity Configuration ===\n")
cat(sprintf("%-12s  %10s  %10s  %10s  %10s  %12s\n",
            "Group", "W0 (log10)", "Wmat(log10)", "Default_W0F", "New_W0F", "w_min(g)"))
for (i in seq_along(fish_rows)) {
  new_w0f <- if (!is.na(w_min_catch_log10[i])) {
    max(original_Fmort_W0[i], w_min_catch_log10[i])
  } else {
    original_Fmort_W0[i]
  }
  w_g <- if (!is.na(w_min_catch_g[i])) sprintf("%.1f", w_min_catch_g[i]) else "default"
  cat(sprintf("%-12s  %10.1f  %10.1f  %10.1f  %10.1f  %12s\n",
              fish_names[i], original_W0[i], original_Wmat[i],
              original_Fmort_W0[i], new_w0f, w_g))
}
cat("\n")

cat("Scenario: repro_eff =", repro_eff_vec, "\n")
cat("Running", length(effort_levels), "effort levels across", n_cores, "cores\n")

# =============================================================================
# RUN SWEEP
# =============================================================================

run_single <- function(repro_eff_vec, effort_val, Groups_base, fish_rows,
                       q_fixed, w_min_catch_log10,
                       sim_time, sst_const, chl_const, isave) {
  Groups_mod <- Groups_base
  Groups_mod$repro_eff[fish_rows] <- repro_eff_vec

  # Apply minimum catch size: set Fmort_W0 to the maximum of default and
  # user-specified threshold (so we never fish below the group's W0 either)
  for (i in seq_along(fish_rows)) {
    if (!is.na(w_min_catch_log10[i])) {
      Groups_mod$Fmort_W0[fish_rows[i]] <- max(
        Groups_mod$Fmort_W0[fish_rows[i]],
        w_min_catch_log10[i]
      )
    }
  }

  Groups_mod <- setFishingParams(Groups_mod,
                                 q_small = q_fixed,
                                 q_med   = q_fixed,
                                 q_large = q_fixed)

  env <- createInputParams(time = sim_time, sst = sst_const, chl = chl_const)
  env <- addFishingEffort(env,
                          effort_small = effort_val,
                          effort_med   = effort_val,
                          effort_large = effort_val)
  zoomss_model(input_params = env, Groups = Groups_mod, isave = isave)
}

t0 <- Sys.time()

cl <- makeCluster(n_cores)
on.exit(stopCluster(cl), add = TRUE)

clusterExport(cl, c("run_single", "repro_eff_vec", "effort_levels",
                     "Groups_base", "fish_rows", "q_fixed", "w_min_catch_log10",
                     "sim_time", "sst_const", "chl_const", "isave", "pkg_path"))
clusterEvalQ(cl, {
  devtools::load_all(pkg_path)
  Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
})

results <- parLapply(cl, seq_along(effort_levels), function(i) {
  run_single(repro_eff_vec, effort_levels[i],
             Groups_base, fish_rows, q_fixed, w_min_catch_log10,
             sim_time, sst_const, chl_const, isave)
})

stopCluster(cl)
on.exit(NULL)

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat("Completed in", elapsed, "minutes\n")

# =============================================================================
# EXTRACT EQUILIBRIUM METRICS (final 200 of 400 years)
# =============================================================================

summary_df <- data.frame()
num_fish   <- length(fish_names)

for (i in seq_along(effort_levels)) {
  mdl       <- results[[i]]
  fish_grps <- mdl$param$fish_grps
  nsave     <- length(mdl$time)
  final_idx <- max(1, floor(nsave * 0.5)):nsave

  for (f in seq_len(num_fish)) {
    fg <- fish_grps[f]

    eq_bm <- mean(apply(
      mdl$abundance[final_idx, fg, , drop = FALSE], 1, function(row) {
        sum(row * mdl$param$w * mdl$param$dx)
      }))

    eq_catch <- if (!is.null(mdl$catch)) {
      mean(mdl$catch[final_idx, f], na.rm = TRUE)
    } else { 0 }

    eq_ssb <- mean(mdl$SSB[final_idx, f], na.rm = TRUE)
    eq_rec <- mean(mdl$recruitment[final_idx, f], na.rm = TRUE)

    summary_df <- rbind(summary_df, data.frame(
      effort       = effort_levels[i],
      fish_group   = fish_names[f],
      eq_catch     = eq_catch,
      eq_biomass   = eq_bm,
      eq_SSB       = eq_ssb,
      eq_recruitment = eq_rec
    ))
  }
}

# Derived metrics
summary_df <- summary_df %>%
  group_by(fish_group) %>%
  mutate(
    B0         = eq_biomass[effort == 0],
    depletion  = eq_biomass / B0,
    yield_norm = eq_catch / max(eq_catch, na.rm = TRUE)
  ) %>%
  ungroup()

# Fish colours from model
fish_colours <- results[[1]]$param$Groups$PlotColour[results[[1]]$param$fish_grps]
names(fish_colours) <- results[[1]]$param$Groups$Species[results[[1]]$param$fish_grps]

# Scenario label for plot titles
w_min_label <- paste0(
  ifelse(is.na(w_min_catch_g), "def", paste0(w_min_catch_g, "g")),
  collapse = ", "
)
scenario_label <- paste0(
  "repro_eff = (", paste(format(repro_eff_vec, scientific = TRUE), collapse = ", "), ")",
  "  |  w_min_catch = (", w_min_label, ")"
)

# Equilibrium averaging period (years) for built-in plot functions
eq_years <- max(results[[1]]$time) / 2

# =============================================================================
# PLOT 0: SELECTIVITY DIAGNOSTIC — what size classes are fished?
# =============================================================================

# Build selectivity profile from a model object
mdl_ref    <- results[[1]]  # Use unfished run for reference
w_log10    <- mdl_ref$param$w_log10
w_grams    <- mdl_ref$param$w
sel_df     <- data.frame()

for (f in seq_len(num_fish)) {
  fg <- mdl_ref$param$fish_grps[f]
  grp_info <- mdl_ref$param$Groups[fg, ]

  # Determine effective Fmort_W0 (with our override applied)
  eff_w0 <- if (!is.na(w_min_catch_log10[f])) {
    max(grp_info$Fmort_W0, w_min_catch_log10[f])
  } else {
    grp_info$Fmort_W0
  }

  # Knife-edge selectivity: 1 where w_log10 is between Fmort_W0 and Fmort_Wmax
  sel <- as.numeric(w_log10 >= eff_w0 & w_log10 <= grp_info$Fmort_Wmax)

  # Also mark the group's size range (W0 to Wmax)
  in_range <- as.numeric(w_log10 >= grp_info$W0 & w_log10 <= grp_info$Wmax)

  sel_df <- rbind(sel_df, data.frame(
    w_log10     = w_log10,
    w_grams     = w_grams,
    selectivity = sel,
    in_range    = in_range,
    fish_group  = fish_names[f],
    stringsAsFactors = FALSE
  ))
}

# Size reference lines
size_refs <- data.frame(
  w_log10 = c(log10(10), log10(100), log10(1000), log10(10000)),
  label   = c("10g", "100g", "1kg", "10kg")
)

# Selectivity plot
p_sel <- ggplot(sel_df %>% filter(in_range > 0),
                aes(x = w_log10, y = selectivity, colour = fish_group)) +
  geom_line(linewidth = 1.2) +
  geom_vline(data = size_refs, aes(xintercept = w_log10),
             linetype = "dotted", colour = "grey50", linewidth = 0.4) +
  geom_text(data = size_refs, aes(x = w_log10, y = 1.08, label = label),
            inherit.aes = FALSE, size = 3, colour = "grey40") +
  scale_colour_manual(values = fish_colours) +
  coord_cartesian(ylim = c(-0.05, 1.15)) +
  theme_bw(base_size = 12) +
  labs(x = expression(log[10] ~ "Body Mass (g)"),
       y = "Selectivity (0 = not fished, 1 = fully fished)",
       colour = "Fish Group",
       title = "Fishing Selectivity by Size",
       subtitle = paste0("Knife-edge at w_min_catch = (", w_min_label, ")"))

# Add maturity markers
wmat_df <- data.frame(
  fish_group = fish_names,
  wmat_log10 = original_Wmat,
  stringsAsFactors = FALSE
)

p_sel <- p_sel +
  geom_point(data = wmat_df, aes(x = wmat_log10, y = 0.5, colour = fish_group),
             shape = 17, size = 3, show.legend = FALSE) +
  annotate("text", x = min(wmat_df$wmat_log10) - 0.3, y = 0.5,
           label = expression(W[mat]), size = 3, colour = "grey30")

# Full range view including unfished sizes
p_sel_full <- ggplot(sel_df, aes(x = w_log10)) +
  geom_ribbon(aes(ymin = 0, ymax = in_range, fill = fish_group), alpha = 0.15) +
  geom_line(aes(y = selectivity, colour = fish_group), linewidth = 1) +
  geom_vline(data = size_refs, aes(xintercept = w_log10),
             linetype = "dotted", colour = "grey50", linewidth = 0.4) +
  geom_text(data = size_refs, aes(y = 1.08, label = label),
            inherit.aes = FALSE, size = 2.5, colour = "grey40") +
  scale_colour_manual(values = fish_colours) +
  scale_fill_manual(values = fish_colours) +
  coord_cartesian(ylim = c(-0.05, 1.15)) +
  theme_bw(base_size = 11) +
  labs(x = expression(log[10] ~ "Body Mass (g)"),
       y = "Selectivity / Presence",
       colour = "Fish Group", fill = "Size Range",
       title = "Full Size Spectrum: Group Range (shaded) vs Fished Range (solid)",
       subtitle = "Shaded = group present; solid line = fished")

# =============================================================================
# PLOTS 1-5: YIELD, BIOMASS, DEPLETION, SSB, RECRUITMENT
# =============================================================================

# 1. Yield-Effort (absolute)
p_yield <- ggplot(summary_df, aes(x = effort, y = eq_catch, colour = fish_group)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = fish_colours) +
  theme_bw(base_size = 12) +
  labs(x = "Effort", y = "Equilibrium Catch (g WW)", colour = "Fish Group",
       title = "Yield-Effort Curve", subtitle = scenario_label)

# 2. Yield-Effort (normalised)
p_yield_norm <- ggplot(summary_df, aes(x = effort, y = yield_norm, colour = fish_group)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = fish_colours) +
  theme_bw(base_size = 12) +
  labs(x = "Effort", y = "Normalised Yield", colour = "Fish Group",
       title = "Normalised Yield-Effort Curve", subtitle = scenario_label)

# 3. Biomass depletion (log scale)
p_biomass <- ggplot(summary_df, aes(x = effort, y = eq_biomass, colour = fish_group)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = fish_colours) +
  scale_y_log10() +
  theme_bw(base_size = 12) +
  labs(x = "Effort", y = "Equilibrium Biomass (log scale)", colour = "Fish Group",
       title = "Biomass Depletion", subtitle = scenario_label)

# 4. B/B0 depletion ratio
p_depletion <- ggplot(summary_df, aes(x = effort, y = depletion, colour = fish_group)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_hline(yintercept = 0.4, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0.2, linetype = "dotted", colour = "red") +
  annotate("text", x = max(effort_levels) * 0.95, y = 0.42,
           label = "B_MSY proxy (0.4)", hjust = 1, size = 3, colour = "grey40") +
  annotate("text", x = max(effort_levels) * 0.95, y = 0.22,
           label = "B_lim (0.2)", hjust = 1, size = 3, colour = "red") +
  scale_colour_manual(values = fish_colours) +
  coord_cartesian(ylim = c(0, 1.05)) +
  theme_bw(base_size = 12) +
  labs(x = "Effort", y = expression(B / B[0]), colour = "Fish Group",
       title = "Stock Depletion", subtitle = scenario_label)

# 5. SSB and Recruitment vs effort
p_ssb <- ggplot(summary_df, aes(x = effort, y = eq_SSB, colour = fish_group)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = fish_colours) +
  scale_y_log10() +
  theme_bw(base_size = 12) +
  labs(x = "Effort", y = "SSB (log scale)", colour = "Fish Group",
       title = "Spawning Stock Biomass")

p_rec <- ggplot(summary_df, aes(x = effort, y = eq_recruitment, colour = fish_group)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = fish_colours) +
  scale_y_log10() +
  theme_bw(base_size = 12) +
  labs(x = "Effort", y = "Recruitment (log scale)", colour = "Fish Group",
       title = "Recruitment")

# =============================================================================
# PLOT 6: SIZE SPECTRA AT KEY EFFORT LEVELS (using built-in plotSizeSpectra)
# =============================================================================
# Use the package's plotSizeSpectra() for each key effort level, averaging over
# the equilibrium period. Faceted comparison across efforts.

spec_efforts <- c(0, 0.2, 0.5, 2.0)
spec_efforts <- spec_efforts[spec_efforts %in% effort_levels]

spec_plots <- lapply(spec_efforts, function(eff) {
  idx <- which(effort_levels == eff)
  plotSizeSpectra(results[[idx]], by = "biomass", n_years = eq_years) +
    labs(title = NULL, subtitle = paste0("Effort = ", eff))
})

p_spectra <- wrap_plots(spec_plots, ncol = 2) +
  plot_annotation(
    title = "Equilibrium Biomass Spectra at Different Effort Levels",
    subtitle = paste0("w_min_catch = (", w_min_label,
                      ") | averaged over final ", eq_years, " years")
  ) +
  plot_layout(guides = "collect")

# =============================================================================
# DISPLAY: SWEEP-LEVEL SUMMARY PLOTS
# =============================================================================

cat("\n=== Sweep Summary Plots ===\n")

# Show selectivity first
print(p_sel + p_sel_full + plot_layout(guides = "collect"))

print(p_yield + p_yield_norm + plot_layout(guides = "collect"))
print(p_biomass + p_depletion + plot_layout(guides = "collect"))
print(p_ssb + p_rec + plot_layout(guides = "collect"))

print(p_spectra)

# =============================================================================
# PLOT 7: MODEL DIAGNOSTICS USING BUILT-IN PLOTTING FUNCTIONS
# =============================================================================
# For key effort levels (unfished, near-MSY, high), show full model diagnostics
# using the package's plotting functions to inspect time series, reproduction,
# growth, and spectra in detail.

cat("\n=== Model Diagnostics (built-in plots) ===\n")

# Select diagnostic effort levels: unfished + near-MSY + high fishing
msy_efforts <- summary_df %>%
  group_by(fish_group) %>%
  summarise(effort_MSY = effort[which.max(eq_catch)], .groups = "drop")
# Use the median MSY effort as the "near-MSY" diagnostic level
median_msy_effort <- median(msy_efforts$effort_MSY)
# Find the closest effort level in our grid
near_msy_effort <- effort_levels[which.min(abs(effort_levels - median_msy_effort))]

diag_efforts <- unique(c(0, near_msy_effort, max(effort_levels)))
diag_efforts <- diag_efforts[diag_efforts %in% effort_levels]

cat("Diagnostic effort levels:", diag_efforts, "\n")
cat("  (unfished, near-MSY at", near_msy_effort, ", max effort)\n\n")

for (eff in diag_efforts) {
  idx <- which(effort_levels == eff)
  mdl <- results[[idx]]
  eff_label <- paste0("Effort = ", eff)

  cat("--- Diagnostics for", eff_label, "---\n")

  # 7a. Biomass time series (all groups, log10)
  p_ts_bm <- plotTimeSeries(mdl, by = "biomass", transform = "log10") +
    labs(title = paste(eff_label, "— Biomass Time Series"))

  # 7b. Biomass community composition (stacked proportional)
  p_ts_fill <- plotTimeSeries(mdl, by = "biomass", type = "fill") +
    labs(title = paste(eff_label, "— Community Composition"))

  print(p_ts_bm + p_ts_fill + plot_layout(guides = "collect"))

  # 7c. Fish-only biomass time series
  p_ts_fish <- plotTimeSeries(mdl, by = "biomass", species = fish_names,
                              transform = "log10") +
    labs(title = paste(eff_label, "— Fish Biomass"))

  # 7d. Fish mortality time series
  p_ts_mort <- plotTimeSeries(mdl, by = "mortality", species = fish_names) +
    labs(title = paste(eff_label, "— Fish Mortality"))

  print(p_ts_fish + p_ts_mort + plot_layout(guides = "collect"))

  # 7e. Reproduction diagnostics (SSB + recruitment time series)
  if ("SSB" %in% names(mdl)) {
    p_ssb_ts <- plotReproduction(mdl, by = "SSB", transform = "log10") +
      labs(title = paste(eff_label, "— SSB"))
    p_rec_ts <- plotReproduction(mdl, by = "recruitment", transform = "log10") +
      labs(title = paste(eff_label, "— Recruitment"))

    print(p_ssb_ts + p_rec_ts + plot_layout(guides = "collect"))
  }

  # 7f. Growth comparison (fish vs zooplankton)
  if ("growth" %in% names(mdl)) {
    p_growth_fish <- plotGrowthComparison(mdl, groups = "fish",
                                          n_years = eq_years) +
      labs(title = paste(eff_label, "— Fish Growth"))
    p_growth_zoo  <- plotGrowthComparison(mdl, groups = "zooplankton",
                                           n_years = eq_years) +
      labs(title = paste(eff_label, "— Zooplankton Growth"))

    print(p_growth_fish + p_growth_zoo + plot_layout(guides = "collect"))
  }

  # 7g. Size spectra (abundance + biomass side by side)
  p_spec_abund <- plotSizeSpectra(mdl, by = "abundance", n_years = eq_years) +
    labs(title = paste(eff_label, "— Abundance Spectrum"))
  p_spec_bm    <- plotSizeSpectra(mdl, by = "biomass", n_years = eq_years) +
    labs(title = paste(eff_label, "— Biomass Spectrum"))

  print(p_spec_abund + p_spec_bm + plot_layout(guides = "collect"))

  # 7h. Reproductive investment by size (if available)
  if ("repro_rate" %in% names(mdl)) {
    p_repro_size <- plotReproSizeSpectra(mdl, n_years = eq_years) +
      labs(title = paste(eff_label, "— Reproductive Investment by Size"))
    print(p_repro_size)
  }
}

# One-off structural diagnostics from unfished model (effort = 0)
mdl_unfished <- results[[which(effort_levels == 0)]]

cat("\n--- Energy Budget & Trophic Structure (unfished) ---\n")

# Energy budget allocation
if ("def_by_prey" %in% names(mdl_unfished)) {
  p_energy  <- plotEnergyBudget(mdl_unfished)
  p_defec   <- plotDefecation(mdl_unfished)
  print(p_energy + p_defec)
}

# Assimilation matrix
if ("assim_by_prey" %in% names(mdl_unfished)) {
  print(plotAssimilationMatrix(mdl_unfished))
}

# =============================================================================
# MSY SUMMARY
# =============================================================================

msy_table <- summary_df %>%
  group_by(fish_group) %>%
  summarise(
    MSY        = max(eq_catch, na.rm = TRUE),
    effort_MSY = effort[which.max(eq_catch)],
    B_at_MSY   = eq_biomass[which.max(eq_catch)],
    B0         = first(B0),
    BdivB0_MSY = eq_biomass[which.max(eq_catch)] / first(B0),
    .groups    = "drop"
  )
cat("\n--- MSY Summary ---\n")
print(msy_table)

# Selectivity summary
cat("\n--- Selectivity Summary ---\n")
sel_summary <- data.frame(
  fish_group    = fish_names,
  W0_log10      = original_W0,
  Wmat_log10    = original_Wmat,
  default_Fmort_W0 = original_Fmort_W0,
  w_min_catch_g = w_min_catch_g,
  effective_Fmort_W0 = sapply(seq_along(fish_rows), function(i) {
    if (!is.na(w_min_catch_log10[i])) {
      max(original_Fmort_W0[i], w_min_catch_log10[i])
    } else {
      original_Fmort_W0[i]
    }
  }),
  fished_fraction = NA
)

# Calculate what fraction of each group's size range is actually fished
for (i in seq_along(fish_rows)) {
  grp <- Groups_base[fish_rows[i], ]
  eff_w0 <- sel_summary$effective_Fmort_W0[i]
  total_range <- grp$Wmax - grp$W0
  fished_range <- grp$Wmax - eff_w0  # could be Fmort_Wmax but usually == Wmax
  sel_summary$fished_fraction[i] <- max(0, fished_range / total_range)
}

print(sel_summary)

cat("\n--- Interpretation ---\n")
cat("fished_fraction: proportion of each group's log10 size range exposed to fishing\n")
cat("  Values < 1.0 mean juveniles/small individuals are sheltered from fishing\n")
cat("  Lower values = more size refuge = greater resilience to fishing pressure\n")
