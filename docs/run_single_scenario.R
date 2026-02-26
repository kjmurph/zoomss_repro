# =============================================================================
# run_single_scenario.R
#
# Quick interactive sweep: one repro_eff scenario across an effort gradient.
# Plots yield curves, biomass depletion, and B/B0 directly.
#
# Usage: source from an interactive R session within the zoomss project.
# =============================================================================

library(here)
devtools::load_all(here())
library(parallel)
library(ggplot2)
library(patchwork)
library(dplyr)

# =============================================================================
# CONFIGURE YOUR SCENARIO HERE
# =============================================================================

# Reproductive efficiency for each fish group: c(Small, Med, Large)
repro_eff_vec <- c(1e-2, 1e-3, 1e-4)

# Effort gradient
effort_levels <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                   0.5, 1.0, 1.5, 2.0, 3.0)

# Catchability (q = 1 means F = effort × selectivity)
q_fixed <- 1.0

# Simulation settings
sim_time  <- seq(0, 400, by = 0.1)
sst_const <- 15
chl_const <- 1.0
isave     <- 5

# Parallel settings
n_cores <- max(1, detectCores() - 1)

# =============================================================================
# RUN SWEEP
# =============================================================================

Groups_base <- getGroups()
fish_rows   <- which(Groups_base$Type == "Fish")
fish_names  <- c("Fish_Small", "Fish_Med", "Fish_Large")
pkg_path    <- here()

cat("Scenario: repro_eff =", repro_eff_vec, "\n")
cat("Running", length(effort_levels), "effort levels across", n_cores, "cores\n")

run_single <- function(repro_eff_vec, effort_val, Groups_base, fish_rows,
                       q_fixed, sim_time, sst_const, chl_const, isave) {
  Groups_mod <- Groups_base
  Groups_mod$repro_eff[fish_rows] <- repro_eff_vec
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
                     "Groups_base", "fish_rows", "q_fixed", "sim_time",
                     "sst_const", "chl_const", "isave", "pkg_path"))
clusterEvalQ(cl, {
  devtools::load_all(pkg_path)
  Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
})

results <- parLapply(cl, seq_along(effort_levels), function(i) {
  run_single(repro_eff_vec, effort_levels[i],
             Groups_base, fish_rows, q_fixed,
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
scenario_label <- paste0("repro_eff = (",
                         paste(format(repro_eff_vec, scientific = TRUE), collapse = ", "), ")")

# =============================================================================
# PLOTS
# =============================================================================

# 1. Yield-Effort (absolute)
p_yield <- ggplot(summary_df, aes(x = effort, y = eq_catch, colour = fish_group)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = fish_colours) +
  theme_bw(base_size = 12) +
  facet_wrap(~fish_group) +
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

# 5. SSB and Recruitment
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

# Display
print(p_yield + p_yield_norm + plot_layout(guides = "collect"))
print(p_biomass + p_depletion + plot_layout(guides = "collect"))
print(p_ssb + p_rec + plot_layout(guides = "collect"))

# MSY summary
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
