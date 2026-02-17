# Diagnose the chl=2 anomaly in the baseline
library(zoomss)

# Load corrected groups (same as baseline script)
Groups <- getGroups()
corrected <- utils::read.csv("data-raw/GroupInputs.csv", stringsAsFactors = FALSE)
shared_cols <- intersect(names(Groups), names(corrected))
for (col in shared_cols) {
  if (!identical(Groups[[col]], corrected[[col]])) {
    Groups[[col]] <- corrected[[col]]
  }
}

cat("=== Comparing getBiomass() vs manual sweep at chl = 1.0 and chl = 2.0 ===\n\n")

for (chl_val in c(1.0, 2.0)) {
  cat("--- chl =", chl_val, "---\n")
  env <- createInputParams(time = seq(0, 400, by = 0.1), sst = 15, chl = chl_val)
  mdl <- zoomss_model(input_params = env, Groups = Groups, isave = 10)

  # Method 1: getBiomass() on full 3D output, then time-average
  Biomass <- getBiomass(mdl, units = "ww")
  cat("Biomass dim:", dim(Biomass), "\n")
  time_vec <- mdl$time
  start_time <- max(0, max(time_vec) - 100)
  time_idx <- which(time_vec >= start_time)
  avg_bm_3d <- apply(Biomass[time_idx, , , drop = FALSE], c(2, 3), mean)
  group_bm_method1 <- rowSums(avg_bm_3d)
  names(group_bm_method1) <- Groups$Species

  # Method 2: averageTimeSeries on abundance, then sweep by w
  avg_N <- averageTimeSeries(mdl, var = "abundance", n_years = 100)
  cat("avg_N dim:", dim(avg_N), "\n")
  w <- mdl$param$w
  cat("w length:", length(w), "\n")
  manual_bm <- sweep(avg_N, 2, w, "*")
  group_bm_method2 <- rowSums(manual_bm)
  names(group_bm_method2) <- Groups$Species

  cat("\nMethod 1 (getBiomass time-avg):\n")
  print(round(group_bm_method1, 4))
  cat("\nMethod 2 (averageTimeSeries + sweep):\n")
  print(round(group_bm_method2, 4))

  # Check for anomalous values
  cat("\nMax abundance per group (avg_N):\n")
  print(round(apply(avg_N, 1, max), 6))

  cat("\nTotal biomass Method 1:", round(sum(group_bm_method1), 4), "\n")
  cat("Total biomass Method 2:", round(sum(group_bm_method2), 4), "\n\n")
}

# Also examine the saved baseline
cat("=== Saved baseline summary ===\n")
bl <- readRDS("calibration/baseline_original_zoomss.rds")
cat("Total biomass per chl level:\n")
print(data.frame(
  chl = bl$chl_levels,
  total_bm = round(bl$total_biomass, 4),
  zoo_fish = round(bl$zoo_fish_ratio, 4)
))

cat("\nZoo proportions:\n")
print(round(bl$zoo_proportions, 4))
