# =============================================================================
# 01_generate_baseline.R
# Generate baseline outputs from the ORIGINAL ZooMSS package (MathMarEcol/zoomss)
# =============================================================================
#
# Purpose: Run the original (pre-energy-budget) ZooMSS across a chlorophyll
#          gradient at constant temperature and save steady-state diagnostics.
#          These outputs serve as the calibration target for the revised model.
#
# IMPORTANT: The original zoomss package ships with several parameter errors
#   in its built-in GroupInputs.rda that must be corrected:
#     - Carbon content: Flagellates (0.14->0.15), Larvaceans (0.01->0.02),
#       OmniCopepods (0.10->0.12), CarnCopepods (0.10->0.12),
#       Euphausiids (0.11->0.12), Salps (0.01->0.02)
#     - Fish_Large maximum size: Wmax (7->6) and Fmort_Wmax (7->6)
#     - Flagellates minimum size: W0 (-12 -> -10.7)
#     - CarnCopepods feeding kernel width: FeedWidth (0.36->0.40)
#   This script loads the original model's default Groups, then patches in
#   the corrected values from our revised data-raw/GroupInputs.csv. This
#   ensures the original model code runs with the correct biological parameters.
#
# Prerequisites:
#   - Install the original zoomss package (uncomment the line below):
#     remotes::install_github("MathMarEcol/zoomss")
#   - The original package should be loaded INSTEAD of the revised version.
#   - The corrected GroupInputs.csv must be available at data-raw/GroupInputs.csv
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

# ── Load original Groups and patch with corrected values ──
# The original package's GroupInputs has incorrect Carbon values and other
# parameter errors. We load its default structure (which the original model
# code expects, including GrossGEscale, Repro, etc.) then overwrite shared
# columns with the corrected values from our revised CSV.

Groups <- getGroups()  # Original package's default Groups

# Load corrected parameter values from the revised CSV
corrected <- utils::read.csv("data-raw/GroupInputs.csv", stringsAsFactors = FALSE)

# Verify species match (order must be identical)
stopifnot(
  "Species mismatch between original and corrected Groups" =
    all(Groups$Species == corrected$Species)
)

# Identify columns that exist in BOTH the original and corrected data frames.
# These are the shared biological parameters that need the corrected values.
# New energy-budget-specific columns (def_high, def_low, f_M, K_growth, etc.)
# only exist in the corrected CSV and are irrelevant to the original model.
shared_cols <- intersect(names(Groups), names(corrected))
patched_cols <- character(0)

for (col in shared_cols) {
  if (!identical(Groups[[col]], corrected[[col]])) {
    cat("  Patching column '", col, "': original -> corrected\n", sep = "")
    Groups[[col]] <- corrected[[col]]
    patched_cols <- c(patched_cols, col)
  }
}

if (length(patched_cols) == 0) {
  cat("  No columns needed patching (original and corrected values match).\n")
} else {
  cat("  Patched", length(patched_cols), "column(s):", paste(patched_cols, collapse = ", "), "\n")
}

cat("\nCorrected key values:\n")
print(data.frame(
  Species  = Groups$Species,
  Carbon   = Groups$Carbon,
  W0       = Groups$W0,
  Wmax     = Groups$Wmax,
  FeedWidth = Groups$FeedWidth,
  Fmort_Wmax = Groups$Fmort_Wmax
))

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

  # Compute wet-weight biomass using getBiomass (time x groups x size)
  Biomass <- getBiomass(mdl, units = "ww")

  # Time-average final avg_years using model time vector
  time_vec <- mdl$time
  start_time <- max(0, max(time_vec) - avg_years)
  time_idx <- which(time_vec >= start_time)

  # Average biomass over time, then sum over size for each group
  avg_biomass <- apply(Biomass[time_idx, , , drop = FALSE], c(2, 3), mean)
  group_biomass <- rowSums(avg_biomass)

  list(
    group_biomass = group_biomass,
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
