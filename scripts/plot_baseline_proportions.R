# Generate and save baseline biomass proportion plot for the chl gradient
# Uses the already-computed baseline from calibration/baseline_original_zoomss.rds

library(ggplot2)
library(dplyr)
library(tidyr)

bl <- readRDS("calibration/baseline_original_zoomss.rds")

# ── Group-level biomass from baseline ──
# Extract per-group biomass across chl levels
species <- bl$zoo_species
chl_levels <- bl$chl_levels
zoo_idx <- bl$zoo_idx
fish_idx <- bl$fish_idx
fish_species <- bl$fish_species

# Build long-form data frame of all group biomass
all_species <- c(species, fish_species)
n_groups <- length(all_species)
n_chl <- length(chl_levels)

bm_matrix <- matrix(NA, n_chl, n_groups)
for (i in seq_len(n_chl)) {
  bm_matrix[i, ] <- bl$raw_results[[i]]$group_biomass
}
colnames(bm_matrix) <- all_species

cat("=== Baseline Group Biomass ===\n")
print(round(bm_matrix, 4))

cat("\n=== Total Biomass ===\n")
cat(paste(chl_levels, ":", round(bl$total_biomass, 2)), sep = "\n")

# Flag: chl levels where fish blow up (Fish_Large > 1000x median)
fish_large_bm <- bm_matrix[, "Fish_Large"]
median_fl <- median(fish_large_bm[fish_large_bm < 100])
unstable <- fish_large_bm > 100
cat("\nUnstable chl levels (Fish_Large blow-up):\n")
cat(paste("chl =", chl_levels[unstable], ": Fish_Large =", round(fish_large_bm[unstable], 1)), sep = "\n")

# ── Biomass proportion plot (all groups) ──
bm_df <- as.data.frame(bm_matrix) %>%
  mutate(chl = chl_levels) %>%
  pivot_longer(-chl, names_to = "species", values_to = "biomass")

# Compute proportion within each chl level
bm_df <- bm_df %>%
  group_by(chl) %>%
  mutate(proportion = biomass / sum(biomass)) %>%
  ungroup()

# Set factor order to match Groups table order
bm_df$species <- factor(bm_df$species, levels = all_species)

# Colour palette from corrected GroupInputs
corrected <- utils::read.csv("data-raw/GroupInputs.csv", stringsAsFactors = FALSE)
all_colors <- setNames(corrected$PlotColour, corrected$Species)

p_all <- ggplot(bm_df, aes(x = factor(chl), y = proportion, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = all_colors) +
  labs(
    title = "Baseline (Original ZooMSS): Biomass Proportions Across Chlorophyll Gradient",
    x = "Chlorophyll (mg/m³)",
    y = "Biomass Proportion",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("baseline_diagnostics/baseline_biomass_proportions.png",
       p_all, width = 12, height = 7, dpi = 150)
cat("\nSaved: baseline_diagnostics/baseline_biomass_proportions.png\n")

# ── Zooplankton-only proportions (excludes fish blow-up effect) ──
zoo_bm_matrix <- bm_matrix[, species]
zoo_long <- as.data.frame(zoo_bm_matrix) %>%
  mutate(chl = chl_levels) %>%
  pivot_longer(-chl, names_to = "species", values_to = "biomass") %>%
  group_by(chl) %>%
  mutate(proportion = biomass / sum(biomass)) %>%
  ungroup()

zoo_long$species <- factor(zoo_long$species, levels = species)
zoo_colors <- all_colors[species]

p_zoo <- ggplot(zoo_long, aes(x = factor(chl), y = proportion, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = zoo_colors) +
  labs(
    title = "Baseline (Original ZooMSS): Zooplankton Biomass Proportions",
    x = "Chlorophyll (mg/m³)",
    y = "Zooplankton Biomass Proportion",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("baseline_diagnostics/baseline_zoo_proportions.png",
       p_zoo, width = 12, height = 7, dpi = 150)
cat("Saved: baseline_diagnostics/baseline_zoo_proportions.png\n")

# ── Log-scale total biomass to show the instability clearly ──
total_df <- data.frame(
  chl = chl_levels,
  total_bm = bl$total_biomass,
  stable = !unstable
)

p_total <- ggplot(total_df, aes(x = chl, y = total_bm)) +
  geom_line(linewidth = 1) +
  geom_point(aes(colour = stable), size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "red"),
                      labels = c("TRUE" = "Stable", "FALSE" = "Anomalous")) +
  labs(
    title = "Baseline Total Biomass Across Chlorophyll Gradient",
    subtitle = "Red points: anomalous Fish_Large blow-up at chl \u2248 1.5\u20133.0",
    x = "Chlorophyll (mg/m³)",
    y = "Total Biomass (g ww, log scale)",
    colour = "Status"
  ) +
  theme_minimal()

ggsave("baseline_diagnostics/baseline_total_biomass.png",
       p_total, width = 10, height = 6, dpi = 150)
cat("Saved: baseline_diagnostics/baseline_total_biomass.png\n")
