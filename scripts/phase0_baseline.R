##############################################################################
# Phase 0: Establish Baselines for ZooMSS Energy Budget Revision
# Run this script before making any code changes.
##############################################################################

library(devtools)
load_all(".")

# Create output directory
dir.create("baseline_diagnostics", showWarnings = FALSE)

##############################################################################
# Step 0.1: Verify single grid cell runs
##############################################################################

cat("=== Step 0.1: Single grid cell verification ===\n")

# Moderate conditions: SST = 20°C, chl = 0.5 mg/m³
env_data <- createInputParams(
  time = seq(0, 100, by = 0.1),
  sst = 20,
  chl = 0.5
)

Groups <- getGroups()
cat("Groups loaded:", nrow(Groups), "groups\n")
cat("Species:", paste(Groups$Species, collapse = ", "), "\n")

cat("\nRunning model (100 years, dt=0.1, isave=2)...\n")
t0 <- Sys.time()
mdl <- zoomss_model(input_params = env_data, Groups = Groups, isave = 2)
t1 <- Sys.time()
cat("Model completed in", round(difftime(t1, t0, units = "secs"), 1), "seconds\n")

# Check output structure
cat("\nModel output names:", paste(names(mdl), collapse = ", "), "\n")
cat("Abundance dims:", paste(dim(mdl$abundance), collapse = " x "), "\n")
cat("Growth dims:", paste(dim(mdl$growth), collapse = " x "), "\n")
cat("Mortality dims:", paste(dim(mdl$mortality), collapse = " x "), "\n")
cat("Time range:", range(mdl$time), "\n")
cat("Weight grid length:", length(mdl$param$w), "\n")
cat("Weight range (log10):", range(log10(mdl$param$w)), "\n")

# Check equilibrium abundance
N_eq <- averageTimeSeries(mdl, var = "abundance", n_years = 10)
cat("\nEquilibrium total abundance per group:\n")
for (i in 1:nrow(N_eq)) {
  cat(sprintf("  %s: %.3e\n", Groups$Species[i], sum(N_eq[i, ])))
}

cat("\n=== Step 0.1 PASSED ===\n")

##############################################################################
# Step 0.2: Run baseline simulations across environmental gradients
# SST: 5, 10, 15, 20, 25, 30°C
# Chl: 0.05, 0.15, 0.5, 1.5, 5.0 mg/m³
# = 30 grid cells total
##############################################################################

cat("\n=== Step 0.2: Environmental gradient simulations ===\n")

sst_gradient <- c(5, 10, 15, 20, 25, 30)
chl_gradient <- c(0.05, 0.15, 0.5, 1.5, 5.0)

# Store results
baseline_results <- list()
baseline_summary <- data.frame()

total_runs <- length(sst_gradient) * length(chl_gradient)
run_count <- 0

for (sst_val in sst_gradient) {
  for (chl_val in chl_gradient) {
    run_count <- run_count + 1
    run_id <- paste0("sst", sst_val, "_chl", chl_val)
    cat(sprintf("\n[%d/%d] Running SST=%.0f°C, chl=%.2f mg/m³...\n",
                run_count, total_runs, sst_val, chl_val))

    env <- createInputParams(
      time = seq(0, 100, by = 0.1),
      sst = sst_val,
      chl = chl_val
    )

    t0 <- Sys.time()
    result <- tryCatch(
      zoomss_model(input_params = env, Groups = Groups, isave = 2),
      error = function(e) {
        cat("  ERROR:", e$message, "\n")
        return(NULL)
      }
    )
    t1 <- Sys.time()

    if (!is.null(result)) {
      elapsed <- round(difftime(t1, t0, units = "secs"), 1)
      cat(sprintf("  Completed in %s seconds\n", elapsed))

      # Extract equilibrium values (average final 10 years)
      N_eq <- averageTimeSeries(result, var = "abundance", n_years = 10)
      gg_eq <- averageTimeSeries(result, var = "growth", n_years = 10)
      Z_eq <- averageTimeSeries(result, var = "mortality", n_years = 10)

      # Compute biomass at equilibrium
      biomass_ww <- getBiomass(result, units = "ww")
      # Average last 10 years of biomass
      n_save <- dim(biomass_ww)[1]
      dt_save <- diff(result$time[1:2])
      n_steps_10yr <- min(round(10 / dt_save), n_save)
      biomass_eq <- apply(biomass_ww[(n_save - n_steps_10yr + 1):n_save, , , drop = FALSE],
                          c(2, 3), mean)
      total_biomass <- rowSums(biomass_eq)

      # Store complete data
      baseline_results[[run_id]] <- list(
        sst = sst_val,
        chl = chl_val,
        N = N_eq,
        gg = gg_eq,
        Z = Z_eq,
        biomass = biomass_eq,
        total_biomass = total_biomass,
        w = result$param$w,
        time = result$time
      )

      # Summary row
      row <- data.frame(
        run_id = run_id,
        sst = sst_val,
        chl = chl_val,
        elapsed_sec = as.numeric(elapsed),
        stringsAsFactors = FALSE
      )
      for (sp in Groups$Species) {
        row[[paste0("biomass_", sp)]] <- total_biomass[which(Groups$Species == sp)]
      }
      baseline_summary <- rbind(baseline_summary, row)
    }
  }
}

cat("\n=== Step 0.2 completed:", nrow(baseline_summary), "of", total_runs, "runs successful ===\n")

##############################################################################
# Step 0.3: Generate and save baseline diagnostics
##############################################################################

cat("\n=== Step 0.3: Generating baseline diagnostics ===\n")

library(ggplot2)

# ---- 0.3a: Size spectra for low/moderate/high productivity ----
cat("Generating size spectra plots...\n")

representative_runs <- c(
  "sst20_chl0.05",  # oligotrophic
  "sst20_chl0.5",   # moderate
  "sst20_chl5"      # eutrophic
)

pdf("baseline_diagnostics/baseline_size_spectra.pdf", width = 14, height = 10)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (rid in representative_runs) {
  if (rid %in% names(baseline_results)) {
    res <- baseline_results[[rid]]
    w_log10 <- log10(res$w)

    plot(NA, xlim = range(w_log10), ylim = c(-5, 25),
         xlab = "log10 body mass (g)", ylab = "log10 abundance",
         main = sprintf("SST=%.0f°C, chl=%.2f", res$sst, res$chl))

    for (i in 1:nrow(res$N)) {
      valid <- res$N[i, ] > 0
      if (any(valid)) {
        lines(w_log10[valid], log10(res$N[i, valid]),
              col = Groups$PlotColour[i], lwd = 2)
      }
    }
    legend("topright", Groups$Species, col = Groups$PlotColour,
           lwd = 2, cex = 0.5, ncol = 2)
  }
}
dev.off()
cat("  Saved: baseline_size_spectra.pdf\n")

# ---- 0.3b: Growth rate profiles ----
cat("Generating growth rate profiles...\n")

pdf("baseline_diagnostics/baseline_growth_rates.pdf", width = 16, height = 12)
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
# Use moderate conditions
res <- baseline_results[["sst20_chl0.5"]]
w_log10 <- log10(res$w)
for (i in 1:nrow(res$gg)) {
  valid <- res$gg[i, ] > 0
  if (any(valid)) {
    plot(w_log10[valid], res$gg[i, valid], type = "l",
         col = Groups$PlotColour[i], lwd = 2,
         xlab = "log10 body mass (g)", ylab = "Growth rate (g/g/yr)",
         main = Groups$Species[i])
  } else {
    plot(1, 1, type = "n", main = paste(Groups$Species[i], "(no growth)"))
  }
}
dev.off()
cat("  Saved: baseline_growth_rates.pdf\n")

# ---- 0.3c: Biomass by group across gradient ----
cat("Generating biomass gradient plots...\n")

# Build a tidy biomass data frame
biomass_df <- data.frame()
for (rid in names(baseline_results)) {
  res <- baseline_results[[rid]]
  for (i in seq_along(Groups$Species)) {
    biomass_df <- rbind(biomass_df, data.frame(
      sst = res$sst,
      chl = res$chl,
      species = Groups$Species[i],
      biomass = res$total_biomass[i],
      stringsAsFactors = FALSE
    ))
  }
}

# Biomass vs SST (at moderate chl)
p_sst <- ggplot(biomass_df[biomass_df$chl == 0.5, ],
                aes(x = sst, y = log10(biomass + 1e-30), colour = species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = setNames(Groups$PlotColour, Groups$Species)) +
  labs(x = "SST (°C)", y = "log10 Biomass (gww/m²)",
       title = "Biomass vs SST (chl = 0.5 mg/m³)") +
  theme_bw()

ggsave("baseline_diagnostics/baseline_biomass_vs_sst.pdf", p_sst, width = 12, height = 7)
cat("  Saved: baseline_biomass_vs_sst.pdf\n")

# Biomass vs Chl (at moderate SST)
p_chl <- ggplot(biomass_df[biomass_df$sst == 20, ],
                aes(x = log10(chl), y = log10(biomass + 1e-30), colour = species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = setNames(Groups$PlotColour, Groups$Species)) +
  labs(x = "log10 Chlorophyll (mg/m³)", y = "log10 Biomass (gww/m²)",
       title = "Biomass vs Chlorophyll (SST = 20°C)") +
  theme_bw()

ggsave("baseline_diagnostics/baseline_biomass_vs_chl.pdf", p_chl, width = 12, height = 7)
cat("  Saved: baseline_biomass_vs_chl.pdf\n")

# ---- 0.3d: Proportional biomass along gradients ----
cat("Generating proportional biomass plots...\n")

# Compute proportion
biomass_df_prop <- biomass_df %>%
  dplyr::group_by(sst, chl) %>%
  dplyr::mutate(total = sum(biomass),
                proportion = biomass / total) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(species = factor(species, levels = Groups$Species))

# Proportion vs SST (at moderate chl)
p_prop_sst <- ggplot(biomass_df_prop[biomass_df_prop$chl == 0.5, ],
                     aes(x = sst, y = proportion, fill = species)) +
  geom_area(position = "fill") +
  scale_fill_manual(values = setNames(Groups$PlotColour, Groups$Species)) +
  labs(x = "SST (°C)", y = "Biomass proportion",
       title = "Community composition vs SST (chl = 0.5 mg/m³)") +
  theme_bw()

ggsave("baseline_diagnostics/baseline_proportion_vs_sst.pdf", p_prop_sst, width = 12, height = 7)
cat("  Saved: baseline_proportion_vs_sst.pdf\n")

# Proportion vs Chl (at moderate SST)
p_prop_chl <- ggplot(biomass_df_prop[biomass_df_prop$sst == 20, ],
                     aes(x = log10(chl), y = proportion, fill = species)) +
  geom_area(position = "fill") +
  scale_fill_manual(values = setNames(Groups$PlotColour, Groups$Species)) +
  labs(x = "log10 Chlorophyll (mg/m³)", y = "Biomass proportion",
       title = "Community composition vs Chlorophyll (SST = 20°C)") +
  theme_bw()

ggsave("baseline_diagnostics/baseline_proportion_vs_chl.pdf", p_prop_chl, width = 12, height = 7)
cat("  Saved: baseline_proportion_vs_chl.pdf\n")

# ---- 0.3e: Predation mortality profiles ----
cat("Generating predation mortality profiles...\n")

pdf("baseline_diagnostics/baseline_predation_mortality.pdf", width = 16, height = 12)
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
res <- baseline_results[["sst20_chl0.5"]]
w_log10 <- log10(res$w)
for (i in 1:nrow(res$Z)) {
  valid <- res$Z[i, ] > 0
  if (any(valid)) {
    plot(w_log10[valid], res$Z[i, valid], type = "l",
         col = Groups$PlotColour[i], lwd = 2,
         xlab = "log10 body mass (g)", ylab = "Mortality (1/yr)",
         main = Groups$Species[i])
  } else {
    plot(1, 1, type = "n", main = paste(Groups$Species[i], "(no mortality)"))
  }
}
dev.off()
cat("  Saved: baseline_predation_mortality.pdf\n")

# ---- 0.3f: Save numerical reference data ----
cat("Saving numerical baseline reference...\n")

# Use the moderate conditions as the primary reference
ref <- baseline_results[["sst20_chl0.5"]]

baseline_data <- list(
  commit = "d07b5ab",
  date = Sys.Date(),
  # Primary reference (moderate conditions)
  gg = ref$gg,
  N = ref$N,
  Z = ref$Z,
  biomass = ref$biomass,
  total_biomass = ref$total_biomass,
  w = ref$w,
  # Full gradient results
  all_results = baseline_results,
  summary = baseline_summary,
  # Environmental conditions
  sst_gradient = sst_gradient,
  chl_gradient = chl_gradient,
  # Group info
  Groups = Groups
)

saveRDS(baseline_data, "baseline_diagnostics/baseline_reference.rds")
cat("  Saved: baseline_reference.rds\n")

# Also save the biomass summary as CSV for easy inspection
write.csv(baseline_summary, "baseline_diagnostics/baseline_biomass_summary.csv",
          row.names = FALSE)
cat("  Saved: baseline_biomass_summary.csv\n")

cat("\n=== Step 0.3 completed ===\n")

##############################################################################
# Step 0.4: Document the baseline
##############################################################################
cat("\n=== Phase 0 COMPLETE ===\n")
cat("All baseline diagnostics saved to baseline_diagnostics/\n")
cat("Proceed to Phase 1 only after reviewing the baseline outputs.\n")
