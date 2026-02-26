# =============================================================================
# inspect_model_dynamics.R
#
# Diagnostic script to inspect ZooMSS model dynamics for a selected repro_eff
# combination. Runs two models (constant and seasonal forcing) and produces
# detailed time series, stability, and community composition diagnostics.
#
# Usage: source from an interactive R session within the zoomss project.
# =============================================================================

library(here)
devtools::load_all(here())
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# =============================================================================
# CONFIGURE SCENARIO
# =============================================================================

# Reproductive efficiency: c(Small, Med, Large)
repro_eff_vec <- c(1e-5, 1e-5, 1e-8)

# Catchability
q_fixed <- 1.0

# Effort (set to 0 for unfished baseline, or a fishing level to inspect)
effort_val <- 0

# Simulation settings
sim_years <- 400
sim_time  <- seq(0, sim_years, by = 0.1)
sst_const <- 15
chl_const <- 1.0
isave     <- 1  # Save every time step for full dynamics

# Seasonal forcing (sinusoidal SST and Chl with annual cycle)
seasonal_amplitude_sst <- 3     # +/- degrees around sst_const
seasonal_amplitude_chl <- 0.5   # +/- mg/m3 around chl_const
sst_seasonal <- sst_const + seasonal_amplitude_sst * sin(2 * pi * sim_time)
chl_seasonal <- pmax(0.1, chl_const + seasonal_amplitude_chl * sin(2 * pi * sim_time + pi))
# Chl lags SST by half a cycle (bloom after winter mixing)

# =============================================================================
# RUN MODELS
# =============================================================================

Groups_base <- getGroups()
fish_rows   <- which(Groups_base$Type == "Fish")

Groups_mod <- Groups_base
Groups_mod$repro_eff[fish_rows] <- repro_eff_vec
Groups_mod <- setFishingParams(Groups_mod,
                                q_small = q_fixed,
                                q_med   = q_fixed,
                                q_large = q_fixed)

scenario_label <- paste0("repro_eff = (",
                         paste(format(repro_eff_vec, scientific = TRUE), collapse = ", "),
                         "), effort = ", effort_val)

cat("=== Running constant forcing model ===\n")
cat("Scenario:", scenario_label, "\n")
t0 <- Sys.time()

env_const <- createInputParams(time = sim_time, sst = sst_const, chl = chl_const)
env_const <- addFishingEffort(env_const,
                               effort_small = effort_val,
                               effort_med   = effort_val,
                               effort_large = effort_val)
mdl_const <- zoomss_model(input_params = env_const, Groups = Groups_mod, isave = isave)

cat("Constant model:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

cat("\n=== Running seasonal forcing model ===\n")
t1 <- Sys.time()

env_seas <- createInputParams(time = sim_time, sst = sst_seasonal, chl = chl_seasonal)
env_seas <- addFishingEffort(env_seas,
                              effort_small = effort_val,
                              effort_med   = effort_val,
                              effort_large = effort_val)
mdl_seas <- zoomss_model(input_params = env_seas, Groups = Groups_mod, isave = isave)

cat("Seasonal model:", round(difftime(Sys.time(), t1, units = "mins"), 1), "min\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Compute biomass time series for each functional group
get_group_biomass_ts <- function(mdl) {
  ntime   <- length(mdl$time)
  ngroups <- nrow(mdl$param$Groups)
  w       <- mdl$param$w
  dx      <- mdl$param$dx

  bm <- matrix(NA, nrow = ntime, ncol = ngroups)
  for (g in 1:ngroups) {
    bm[, g] <- apply(mdl$abundance[, g, , drop = FALSE], 1, function(row) {
      sum(row * w * dx)
    })
  }
  colnames(bm) <- mdl$param$Groups$Species
  as.data.frame(bm) %>%
    mutate(time = mdl$time) %>%
    pivot_longer(-time, names_to = "group", values_to = "biomass")
}

# Compute biomass by type (Fish, Zooplankton, etc.)
get_type_biomass_ts <- function(mdl) {
  ntime   <- length(mdl$time)
  w       <- mdl$param$w
  dx      <- mdl$param$dx
  types   <- unique(mdl$param$Groups$Type)

  type_bm <- list()
  for (tp in types) {
    grp_idx <- which(mdl$param$Groups$Type == tp)
    bm_vec  <- rep(0, ntime)
    for (g in grp_idx) {
      bm_vec <- bm_vec + apply(mdl$abundance[, g, , drop = FALSE], 1, function(row) {
        sum(row * w * dx)
      })
    }
    type_bm[[tp]] <- bm_vec
  }
  as.data.frame(type_bm) %>%
    mutate(time = mdl$time) %>%
    pivot_longer(-time, names_to = "type", values_to = "biomass")
}

# Rolling CV (window in years, given dt between saved steps)
rolling_cv <- function(x, window_years, dt) {
  win <- round(window_years / dt)
  n   <- length(x)
  cv  <- rep(NA, n)
  for (i in win:n) {
    chunk  <- x[(i - win + 1):i]
    cv[i]  <- sd(chunk) / mean(chunk)
  }
  cv
}

# Get group colours from model
get_group_colours <- function(mdl) {
  cols <- mdl$param$Groups$PlotColour
  names(cols) <- mdl$param$Groups$Species
  cols
}

# =============================================================================
# EXTRACT TIME SERIES
# =============================================================================

cat("\nExtracting biomass time series...\n")

# Per-group biomass
bm_const <- get_group_biomass_ts(mdl_const) %>% mutate(forcing = "Constant")
bm_seas  <- get_group_biomass_ts(mdl_seas)  %>% mutate(forcing = "Seasonal")
bm_all   <- bind_rows(bm_const, bm_seas)

# Per-type biomass
type_const <- get_type_biomass_ts(mdl_const) %>% mutate(forcing = "Constant")
type_seas  <- get_type_biomass_ts(mdl_seas)  %>% mutate(forcing = "Seasonal")
type_all   <- bind_rows(type_const, type_seas)

# Total community biomass
total_const <- type_const %>% group_by(time, forcing) %>%
  summarise(biomass = sum(biomass), .groups = "drop")
total_seas <- type_seas %>% group_by(time, forcing) %>%
  summarise(biomass = sum(biomass), .groups = "drop")
total_all <- bind_rows(total_const, total_seas)

# Group colours
group_colours <- get_group_colours(mdl_const)

# Fish-only subset
fish_species <- mdl_const$param$Groups$Species[fish_rows]
bm_fish <- bm_all %>% filter(group %in% fish_species)

# Zooplankton subset
zoo_rows    <- which(mdl_const$param$Groups$Type == "Zooplankton")
zoo_species <- mdl_const$param$Groups$Species[zoo_rows]
bm_zoo <- bm_all %>% filter(group %in% zoo_species)

# =============================================================================
# PLOT 1: FISH BIOMASS THROUGH TIME
# =============================================================================

cat("Generating plots...\n")

p_fish_ts <- ggplot(bm_fish, aes(x = time, y = biomass, colour = group)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  facet_wrap(~ forcing, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = group_colours) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Biomass (g WW)", colour = "Fish Group",
       title = "Fish Biomass Through Time", subtitle = scenario_label)

# Zoomed to final 50 years
p_fish_zoom <- ggplot(bm_fish %>% filter(time >= sim_years - 50),
                      aes(x = time, y = biomass, colour = group)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ forcing, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = group_colours) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Biomass (g WW)", colour = "Fish Group",
       title = paste0("Fish Biomass — Final 50 Years"),
       subtitle = "Check for steady state vs oscillations")

# =============================================================================
# PLOT 2: ZOOPLANKTON BIOMASS THROUGH TIME
# =============================================================================

p_zoo_ts <- ggplot(bm_zoo, aes(x = time, y = biomass, colour = group)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  facet_wrap(~ forcing, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = group_colours) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Biomass (g WW)", colour = "Group",
       title = "Zooplankton Biomass Through Time", subtitle = scenario_label)

p_zoo_zoom <- ggplot(bm_zoo %>% filter(time >= sim_years - 50),
                     aes(x = time, y = biomass, colour = group)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ forcing, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = group_colours) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Biomass (g WW)", colour = "Group",
       title = "Zooplankton Biomass — Final 50 Years")

# =============================================================================
# PLOT 3: TOTAL COMMUNITY BIOMASS
# =============================================================================

p_total <- ggplot(total_all, aes(x = time, y = biomass, colour = forcing)) +
  geom_line(linewidth = 0.5) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Total Community Biomass (g WW)", colour = "Forcing",
       title = "Total Community Biomass Through Time", subtitle = scenario_label)

# =============================================================================
# PLOT 4: BIOMASS BY TYPE (stacked area)
# =============================================================================

p_type_const <- ggplot(type_const, aes(x = time, y = biomass, fill = type)) +
  geom_area(alpha = 0.7) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Biomass (g WW)", fill = "Type",
       title = "Community Composition — Constant Forcing")

p_type_seas <- ggplot(type_seas, aes(x = time, y = biomass, fill = type)) +
  geom_area(alpha = 0.7) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Biomass (g WW)", fill = "Type",
       title = "Community Composition — Seasonal Forcing")

# =============================================================================
# PLOT 5: BIOMASS PROPORTIONS THROUGH TIME
# =============================================================================

prop_const <- type_const %>%
  group_by(time) %>%
  mutate(proportion = biomass / sum(biomass)) %>%
  ungroup()

prop_seas <- type_seas %>%
  group_by(time) %>%
  mutate(proportion = biomass / sum(biomass)) %>%
  ungroup()

p_prop_const <- ggplot(prop_const, aes(x = time, y = proportion, fill = type)) +
  geom_area() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Proportion of Total Biomass", fill = "Type",
       title = "Biomass Proportions — Constant Forcing")

p_prop_seas <- ggplot(prop_seas, aes(x = time, y = proportion, fill = type)) +
  geom_area() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Proportion of Total Biomass", fill = "Type",
       title = "Biomass Proportions — Seasonal Forcing")

# =============================================================================
# PLOT 6: COEFFICIENT OF VARIATION (rolling 10-year window)
# =============================================================================

dt_save <- mdl_const$time[2] - mdl_const$time[1]

cv_fish <- bm_fish %>%
  filter(forcing == "Constant") %>%
  group_by(group) %>%
  arrange(time) %>%
  mutate(cv = rolling_cv(biomass, window_years = 10, dt = dt_save)) %>%
  ungroup()

cv_fish_seas <- bm_fish %>%
  filter(forcing == "Seasonal") %>%
  group_by(group) %>%
  arrange(time) %>%
  mutate(cv = rolling_cv(biomass, window_years = 10, dt = dt_save)) %>%
  ungroup()

cv_all_fish <- bind_rows(
  cv_fish %>% mutate(forcing = "Constant"),
  cv_fish_seas %>% mutate(forcing = "Seasonal")
)

p_cv <- ggplot(cv_all_fish %>% filter(!is.na(cv)),
               aes(x = time, y = cv, colour = group)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ forcing, ncol = 1) +
  scale_colour_manual(values = group_colours) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "CV (10-year rolling window)",
       colour = "Fish Group",
       title = "Fish Biomass Coefficient of Variation",
       subtitle = "CV → 0 indicates steady state; persistent CV indicates oscillations")

# =============================================================================
# PLOT 7: SEASONAL FORCING INPUTS
# =============================================================================

forcing_df <- data.frame(
  time = sim_time,
  SST  = sst_seasonal,
  Chl  = chl_seasonal
) %>%
  filter(time >= sim_years - 10)  # Show final 10 years

p_forcing <- ggplot(forcing_df %>% pivot_longer(-time),
                    aes(x = time, y = value)) +
  geom_line(colour = "steelblue", linewidth = 0.5) +
  facet_wrap(~ name, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Value",
       title = "Seasonal Forcing — Final 10 Years",
       subtitle = paste0("SST: ", sst_const, " ± ", seasonal_amplitude_sst,
                         "°C, Chl: ", chl_const, " ± ", seasonal_amplitude_chl, " mg/m³"))

# =============================================================================
# PLOT 8: SEASONAL RESPONSE — FINAL 10 YEARS
# =============================================================================

p_fish_seasonal <- ggplot(bm_fish %>%
                            filter(forcing == "Seasonal",
                                   time >= sim_years - 10),
                          aes(x = time, y = biomass, colour = group)) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = group_colours) +
  theme_bw(base_size = 11) +
  labs(x = "Time (years)", y = "Biomass (g WW)", colour = "Fish Group",
       title = "Fish Seasonal Dynamics — Final 10 Years",
       subtitle = "Larger fish should show dampened seasonal response")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

# Equilibrium biomass (final 200 years, constant forcing only)
final_idx <- which(mdl_const$time >= sim_years / 2)

eq_summary <- bm_const %>%
  filter(time >= sim_years / 2) %>%
  group_by(group) %>%
  summarise(
    mean_biomass = mean(biomass),
    sd_biomass   = sd(biomass),
    cv_biomass   = sd(biomass) / mean(biomass),
    min_biomass  = min(biomass),
    max_biomass  = max(biomass),
    .groups      = "drop"
  ) %>%
  arrange(desc(mean_biomass))

cat("\n=== Equilibrium Biomass Summary (Constant, final 200 yr) ===\n")
print(eq_summary, n = 20)

# Community proportions at equilibrium
eq_type <- type_const %>%
  filter(time >= sim_years / 2) %>%
  group_by(type) %>%
  summarise(mean_biomass = mean(biomass), .groups = "drop") %>%
  mutate(proportion = mean_biomass / sum(mean_biomass)) %>%
  arrange(desc(proportion))

cat("\n=== Community Composition at Equilibrium ===\n")
print(eq_type)

# Seasonal CV comparison
seas_cv <- bm_fish %>%
  filter(forcing == "Seasonal", time >= sim_years / 2) %>%
  group_by(group) %>%
  summarise(cv = sd(biomass) / mean(biomass), .groups = "drop")

const_cv <- bm_fish %>%
  filter(forcing == "Constant", time >= sim_years / 2) %>%
  group_by(group) %>%
  summarise(cv = sd(biomass) / mean(biomass), .groups = "drop")

cv_compare <- const_cv %>%
  rename(cv_constant = cv) %>%
  left_join(seas_cv %>% rename(cv_seasonal = cv), by = "group") %>%
  mutate(cv_ratio = cv_seasonal / cv_constant)

cat("\n=== Fish CV Comparison (final 200 yr) ===\n")
cat("cv_ratio >> 1 means seasonal forcing adds substantial variability\n")
cat("cv_constant near 0 means model has reached true steady state\n\n")
print(cv_compare)

# =============================================================================
# DISPLAY PLOTS
# =============================================================================

cat("\n=== Displaying plots ===\n")

print(p_fish_ts)
readline("Press [Enter] for next plot...")

print(p_fish_zoom)
readline("Press [Enter] for next plot...")

print(p_zoo_ts)
readline("Press [Enter] for next plot...")

print(p_zoo_zoom)
readline("Press [Enter] for next plot...")

print(p_total)
readline("Press [Enter] for next plot...")

print(p_type_const / p_type_seas)
readline("Press [Enter] for next plot...")

print(p_prop_const / p_prop_seas)
readline("Press [Enter] for next plot...")

print(p_cv)
readline("Press [Enter] for next plot...")

print(p_forcing)
readline("Press [Enter] for next plot...")

print(p_fish_seasonal)

cat("\n=== Done ===\n")
cat("All model objects (mdl_const, mdl_seas) and summary data remain in memory.\n")
cat("Use summary_df, eq_summary, eq_type, cv_compare for further analysis.\n")
