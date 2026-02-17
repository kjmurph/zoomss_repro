# Quick check: chl=1.5 and chl=3.0 to bracket the instability
library(zoomss)

Groups <- getGroups()
corrected <- utils::read.csv("data-raw/GroupInputs.csv", stringsAsFactors = FALSE)
shared_cols <- intersect(names(Groups), names(corrected))
for (col in shared_cols) {
  if (!identical(Groups[[col]], corrected[[col]])) Groups[[col]] <- corrected[[col]]
}

for (chl_val in c(1.5, 3.0)) {
  cat("--- chl =", chl_val, "---\n")
  env <- createInputParams(time = seq(0, 400, by = 0.1), sst = 15, chl = chl_val)
  mdl <- zoomss_model(input_params = env, Groups = Groups, isave = 10)
  Biomass <- getBiomass(mdl, units = "ww")
  time_vec <- mdl$time
  time_idx <- which(time_vec >= max(time_vec) - 100)
  avg_bm <- apply(Biomass[time_idx, , , drop = FALSE], c(2, 3), mean)
  group_bm <- rowSums(avg_bm)
  names(group_bm) <- Groups$Species
  cat("Total:", round(sum(group_bm), 2), "\n")
  print(round(group_bm, 4))
  cat("\n")
}
