# ZooMSS Reproduction Visualization Suite
# Creates comprehensive plots using ZooMSS built-in plotting functions

library(devtools)
load_all()
library(ggplot2)
library(patchwork)

cat("=== ZOOMSS REPRODUCTION VISUALIZATION SUITE ===\n")

# Load results from long-term simulation
if(file.exists("long_term_reproduction_results.RData")) {
  load("long_term_reproduction_results.RData")
  cat("✅ Long-term results loaded successfully\n")
} else {
  cat("❌ No long-term results found. Run comprehensive_reproduction_test.R first\n")
  stop("Missing results file")
}

cat("Creating comprehensive visualization plots...\n\n")

# === 1. SIZE SPECTRA ANALYSIS ===
cat("1. Creating size spectra plots...\n")

# Abundance size spectra
p1_abundance <- plotSizeSpectra(results, by = "abundance", n_years = 10) +
  labs(title = "Final Abundance Size Spectra (Last 10 Years)",
       subtitle = "Log-log plot showing abundance vs body weight") +
  theme_minimal()

# Biomass size spectra  
p1_biomass <- plotSizeSpectra(results, by = "biomass", n_years = 10) +
  labs(title = "Final Biomass Size Spectra (Last 10 Years)",
       subtitle = "Log-log plot showing biomass vs body weight") +
  theme_minimal()

# Growth size spectra
p1_growth <- plotSizeSpectra(results, by = "growth", n_years = 10) +
  labs(title = "Final Growth Rate Size Spectra (Last 10 Years)",
       subtitle = "Log-log plot showing growth rates vs body weight") +
  theme_minimal()

# Mortality size spectra
p1_mortality <- plotSizeSpectra(results, by = "mortality", n_years = 10) +
  labs(title = "Final Mortality Rate Size Spectra (Last 10 Years)", 
       subtitle = "Log-log plot showing mortality rates vs body weight") +
  theme_minimal()

# Combine size spectra plots
size_spectra_combined <- (p1_abundance + p1_biomass) / (p1_growth + p1_mortality)

ggsave("plots_size_spectra_combined.png", size_spectra_combined, 
       width = 16, height = 12, dpi = 300)
cat("✅ Size spectra plots saved to 'plots_size_spectra_combined.png'\n")

# === 2. TIME SERIES ANALYSIS ===
cat("2. Creating time series plots...\n")

# Abundance time series (log scale)
p2_abundance <- plotTimeSeries(results, by = "abundance", transform = "log10") +
  labs(title = "Abundance Time Series (200 Years)",
       subtitle = "Log10 abundance over time for all functional groups") +
  theme_minimal()

# Biomass time series (log scale)
p2_biomass <- plotTimeSeries(results, by = "biomass", transform = "log10") +
  labs(title = "Biomass Time Series (200 Years)",
       subtitle = "Log10 biomass over time for all functional groups") +
  theme_minimal()

# Growth rate time series
p2_growth <- plotTimeSeries(results, by = "growth") +
  labs(title = "Growth Rate Time Series (200 Years)",
       subtitle = "Average growth rates over time") +
  theme_minimal()

# Mortality time series  
p2_mortality <- plotTimeSeries(results, by = "mortality") +
  labs(title = "Mortality Rate Time Series (200 Years)",
       subtitle = "Average mortality rates over time") +
  theme_minimal()

# Save individual time series plots
ggsave("plots_abundance_timeseries.png", p2_abundance, width = 12, height = 8, dpi = 300)
ggsave("plots_biomass_timeseries.png", p2_biomass, width = 12, height = 8, dpi = 300)
ggsave("plots_growth_timeseries.png", p2_growth, width = 12, height = 8, dpi = 300)
ggsave("plots_mortality_timeseries.png", p2_mortality, width = 12, height = 8, dpi = 300)

cat("✅ Time series plots saved\n")

# === 3. FISH-SPECIFIC ANALYSIS ===
cat("3. Creating fish-specific plots...\n")

fish_names <- Groups$Species[Groups$Type == "Fish"]

# Fish abundance time series
p3_fish_abundance <- plotTimeSeries(results, by = "abundance", 
                                   species = fish_names, transform = "log10") +
  labs(title = "Fish Abundance Time Series (200 Years)", 
       subtitle = "Log10 abundance for fish groups only") +
  theme_minimal()

# Fish biomass time series
p3_fish_biomass <- plotTimeSeries(results, by = "biomass",
                                 species = fish_names, transform = "log10") +
  labs(title = "Fish Biomass Time Series (200 Years)",
       subtitle = "Log10 biomass for fish groups only") +
  theme_minimal()

# Fish biomass proportional (stacked)
p3_fish_prop <- plotTimeSeries(results, by = "biomass", type = "fill",
                              species = fish_names) +
  labs(title = "Fish Biomass Proportions (200 Years)",
       subtitle = "Proportional biomass composition of fish community") +
  theme_minimal()

# Save fish-specific plots
ggsave("plots_fish_abundance.png", p3_fish_abundance, width = 12, height = 8, dpi = 300)
ggsave("plots_fish_biomass.png", p3_fish_biomass, width = 12, height = 8, dpi = 300)
ggsave("plots_fish_proportions.png", p3_fish_prop, width = 12, height = 8, dpi = 300)

cat("✅ Fish-specific plots saved\n")

# === 4. ZOOPLANKTON REFERENCE ANALYSIS ===
cat("4. Creating zooplankton reference plots...\n")

zoo_names <- Groups$Species[Groups$Type == "Zooplankton"]

# Zooplankton abundance time series
p4_zoo_abundance <- plotTimeSeries(results, by = "abundance",
                                  species = zoo_names, transform = "log10") +
  labs(title = "Zooplankton Abundance Time Series (200 Years)",
       subtitle = "Log10 abundance for zooplankton groups (reference dynamics)") +
  theme_minimal()

# Zooplankton biomass time series
p4_zoo_biomass <- plotTimeSeries(results, by = "biomass", 
                                species = zoo_names, transform = "log10") +
  labs(title = "Zooplankton Biomass Time Series (200 Years)",
       subtitle = "Log10 biomass for zooplankton groups (reference dynamics)") +
  theme_minimal()

# Save zooplankton plots
ggsave("plots_zooplankton_abundance.png", p4_zoo_abundance, width = 12, height = 8, dpi = 300)
ggsave("plots_zooplankton_biomass.png", p4_zoo_biomass, width = 12, height = 8, dpi = 300)

cat("✅ Zooplankton reference plots saved\n")

# === 5. REPRODUCTION-SPECIFIC ANALYSIS ===
cat("5. Creating reproduction-specific plots...\n")

if(!is.null(results$SSB)) {
  # Create SSB time series plot
  ssb_data <- data.frame(
    time = rep(results$time, length(fish_names)),
    SSB = as.vector(results$SSB),
    Fish_Group = rep(fish_names, each = length(results$time))
  )
  
  p5_ssb <- ggplot(ssb_data, aes(x = time, y = SSB, color = Fish_Group)) +
    geom_line(size = 1) +
    scale_y_log10() +
    labs(title = "Spawning Stock Biomass Time Series (200 Years)",
         subtitle = "Log10 SSB showing reproductive capacity over time",
         x = "Time (years)", y = "SSB (log10 scale)", color = "Fish Group") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("plots_ssb_timeseries.png", p5_ssb, width = 12, height = 8, dpi = 300)
  cat("✅ SSB time series plot saved\n")
  
  # SSB statistics
  cat("\nSSB Statistics Summary:\n")
  for(i in seq_along(fish_names)) {
    ssb_group <- results$SSB[i, ]
    ssb_group <- ssb_group[ssb_group > 0]  # Remove zeros
    if(length(ssb_group) > 0) {
      cat(sprintf("%s SSB: mean = %.2e, min = %.2e, max = %.2e, CV = %.2f\n",
                  fish_names[i], mean(ssb_group), min(ssb_group), max(ssb_group),
                  sd(ssb_group)/mean(ssb_group)))
    }
  }
} else {
  cat("❌ No SSB data available - reproduction may not be enabled\n")
}

# === 6. EQUILIBRATION ANALYSIS ===
cat("\n6. Creating equilibration analysis...\n")

# Calculate rolling means for stability assessment
n_saves <- length(results$time)
window_size <- min(20, floor(n_saves/10))  # 20 saves or 10% of simulation

# Fish biomass equilibration
fish_indices <- which(Groups$Type == "Fish")
fish_biomass_total <- matrix(NA, nrow = n_saves, ncol = length(fish_indices))

for(i in seq_along(fish_indices)) {
  grp_idx <- fish_indices[i]
  for(t in 1:n_saves) {
    fish_biomass_total[t, i] <- sum(results$biomass[t, grp_idx, ], na.rm = TRUE)
  }
}

# Create equilibration plot
equil_data <- data.frame(
  time = rep(results$time, length(fish_names)),
  biomass = as.vector(fish_biomass_total),
  Fish_Group = rep(fish_names, each = n_saves)
)

p6_equil <- ggplot(equil_data, aes(x = time, y = biomass, color = Fish_Group)) +
  geom_line(size = 1) +
  scale_y_log10() +
  labs(title = "Fish Biomass Equilibration (200 Years)",
       subtitle = "Total biomass per fish group - assess model stability",
       x = "Time (years)", y = "Total Biomass (log10 g/m³)", color = "Fish Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("plots_fish_equilibration.png", p6_equil, width = 12, height = 8, dpi = 300)
cat("✅ Equilibration plot saved\n")

# === SUMMARY ===
cat("\n=== VISUALIZATION COMPLETE ===\n")
cat("All plots saved to current directory with 'plots_' prefix\n")
cat("\nFiles created:\n")
cat("- plots_size_spectra_combined.png (4-panel size spectra)\n")
cat("- plots_abundance_timeseries.png (all groups abundance)\n") 
cat("- plots_biomass_timeseries.png (all groups biomass)\n")
cat("- plots_growth_timeseries.png (growth rates)\n")
cat("- plots_mortality_timeseries.png (mortality rates)\n")
cat("- plots_fish_abundance.png (fish abundance only)\n")
cat("- plots_fish_biomass.png (fish biomass only)\n")
cat("- plots_fish_proportions.png (fish community composition)\n")
cat("- plots_zooplankton_abundance.png (zooplankton reference)\n")
cat("- plots_zooplankton_biomass.png (zooplankton reference)\n")
if(!is.null(results$SSB)) {
  cat("- plots_ssb_timeseries.png (spawning stock biomass)\n")
}
cat("- plots_fish_equilibration.png (stability assessment)\n")

cat("\nRecommendations for analysis:\n")
cat("1. Check equilibration plots for stability (last 25% should be relatively stable)\n")
cat("2. Compare fish dynamics to zooplankton reference patterns\n") 
cat("3. Examine size spectra for realistic slopes and distributions\n")
cat("4. Check SSB values for biological realism\n")
cat("5. Assess whether 200 years was sufficient for equilibration\n")