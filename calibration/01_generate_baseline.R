# =============================================================================
# 01_generate_baseline.R
# Generate baseline outputs from the ORIGINAL ZooMSS package (MathMarEcol/zoomss)
# =============================================================================
#
# Purpose: Run the original (pre-energy-budget) ZooMSS across a chlorophyll
#          gradient at constant temperature and save steady-state diagnostics.
#          These outputs serve as the calibration target for the revised model.
#
# Prerequisites:
#   - Install the original zoomss package (uncomment the line below):
#     remotes::install_github("MathMarEcol/zoomss")
#   - The original package should be loaded INSTEAD of the revised version.
#
# Output: calibration/baseline_original_zoomss.rds
# =============================================================================

library(zoomss)
library(future.apply)

# ── Set up parallel execution ──
plan(multisession, workers = parallelly::availableCores() - 1)
cat("Using", parallelly::availableCores() - 1, "parallel workers\n")

# ── Define chlorophyll gradient ──
chl_levels <- c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)
sst_constant <- 15  # Constant temperature (°C)
sim_years <- 400    # Minimum for oscillating steady state
avg_years <- 100    # Average over final 100 years

cat("Running original ZooMSS across", length(chl_levels), "chlorophyll levels\n")
cat("Simulation:", sim_years, "years | Averaging final:", avg_years, "years\n\n")

# ── Load default groups from original package ──
Groups <- getGroups()
zoo_idx <- which(Groups$Type == "Zooplankton")
fish_idx <- which(Groups$Type == "Fish")
n_zoo <- length(zoo_idx)
n_fish <- length(fish_idx)

cat("Zooplankton groups (", n_zoo, "):", Groups$Species[zoo_idx], "\n")
cat("Fish groups (", n_fish, "):", Groups$Species[fish_idx], "\n\n")

# ── Run model across chlorophyll gradient (parallelised) ──
results <- future_lapply(chl_levels, function(chl) {
  cat("  Starting chl =", chl, "mg/m³\n")

  env <- createInputParams(
    time = seq(0, sim_years, by = 0.1),
    sst  = sst_constant,
    chl  = chl
  )

  mdl <- zoomss_model(input_params = env, Groups = Groups, isave = 10)

  # Extract time-averaged abundance over final avg_years
  avg_N <- averageTimeSeries(mdl, var = "N", n_years = avg_years)
  # avg_N is a 2D matrix: groups x size_classes


  # Compute biomass: abundance * weight at each size class
  w <- mdl$param$w
  avg_biomass <- sweep(avg_N, 2, w, "*")  # groups x size_classes

  # Sum across sizes to get total biomass per group
  group_biomass <- rowSums(avg_biomass)

  list(
    group_biomass = group_biomass,
    avg_N         = avg_N,
    avg_biomass   = avg_biomass,
    chl           = chl
  )
}, future.seed = TRUE)

# ── Extract diagnostics into matrices ──
zoo_prop_matrix   <- matrix(NA, nrow = length(chl_levels), ncol = n_zoo)
fish_biomass_matrix <- matrix(NA, nrow = length(chl_levels), ncol = n_fish)
total_biomass_vec <- numeric(length(chl_levels))
zoo_fish_ratio_vec <- numeric(length(chl_levels))

for (i in seq_along(chl_levels)) {
  gb <- results[[i]]$group_biomass

  zoo_bm  <- gb[zoo_idx]
  fish_bm <- gb[fish_idx]

  zoo_prop_matrix[i, ]   <- zoo_bm / sum(zoo_bm)
  fish_biomass_matrix[i, ] <- fish_bm
  total_biomass_vec[i]   <- sum(gb)
  zoo_fish_ratio_vec[i]  <- sum(zoo_bm) / sum(fish_bm)
}

# Add column names for clarity
colnames(zoo_prop_matrix)   <- Groups$Species[zoo_idx]
colnames(fish_biomass_matrix) <- Groups$Species[fish_idx]

# ── Build baseline object ──
baseline <- list(
  chl_levels      = chl_levels,
  sst             = sst_constant,
  sim_years       = sim_years,
  avg_years       = avg_years,
  zoo_idx         = zoo_idx,
  fish_idx        = fish_idx,
  zoo_species     = Groups$Species[zoo_idx],
  fish_species    = Groups$Species[fish_idx],
  zoo_proportions = zoo_prop_matrix,    # matrix: n_chl x n_zoo_groups
  fish_biomass    = fish_biomass_matrix, # matrix: n_chl x n_fish_groups
  total_biomass   = total_biomass_vec,  # vector: n_chl
  zoo_fish_ratio  = zoo_fish_ratio_vec, # vector: n_chl
  raw_results     = results             # full results for further analysis
)

# ── Save ──
dir.create("calibration", showWarnings = FALSE, recursive = TRUE)
saveRDS(baseline, "calibration/baseline_original_zoomss.rds")

cat("\n=== Baseline Summary ===\n")
cat("Saved to: calibration/baseline_original_zoomss.rds\n\n")

cat("Total biomass across gradient:\n")
print(data.frame(
  chl = chl_levels,
  total_biomass = round(total_biomass_vec, 4),
  zoo_fish_ratio = round(zoo_fish_ratio_vec, 2)
))

cat("\nZooplankton proportions:\n")
print(round(zoo_prop_matrix, 3))

cat("\nDone.\n")
