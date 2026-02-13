# Phase 1 Verification Script
library(devtools)
load_all(".")

cat("=== Phase 1 Verification ===\n\n")

# Step 1: Verify package loads and validation works
cat("--- Test 1: Package loads and validation ---\n")
Groups <- getGroups()
cat("Groups loaded:", ncol(Groups), "columns,", nrow(Groups), "rows\n")
cat("New columns present:", all(c("AssimCategory", "Kappa", "MetabConst", "MetabExp", "StarvSens") %in% names(Groups)), "\n")
validateGroups(Groups)

# Step 2: Verify assim_eff values
cat("\n--- Test 2: Verify assim_eff values ---\n")
env_data <- createInputParams(
  time = seq(0, 10, by = 0.1),
  sst = 20,
  chl = 0.5
)
param <- zoomss_params(Groups, env_data, isave = 2)
model <- zoomss_setup(param)

cat("assim_eff[1,1] (Flagellates as prey):", model$assim_eff[1, 1], "\n")
cat("  Expected: 0.75 * 0.14 =", 0.75 * 0.14, "\n")
cat("  Old value would be: 2.5 * 0.14 =", 2.5 * 0.14, "\n")

cat("assim_eff[4,1] (OmniCopepods as prey):", model$assim_eff[4, 1], "\n")
cat("  Expected: 0.70 * 0.10 =", 0.70 * 0.10, "\n")

cat("assim_eff[10,1] (Fish_Small as prey):", model$assim_eff[10, 1], "\n")
cat("  Expected: 0.80 * 0.10 =", 0.80 * 0.10, "\n")

cat("assim_eff[9,1] (Jellyfish as prey):", model$assim_eff[9, 1], "\n")
cat("  Expected: 0.50 * 0.005 =", 0.50 * 0.005, "\n")

# Step 3: Verify carbon_i, kappa, metab_cost
cat("\n--- Test 3: Verify new model arrays ---\n")
cat("carbon_i:", model$carbon_i, "\n")
cat("kappa[1,1] (Flagellates):", model$kappa[1, 1], "  Expected: 0.7\n")
cat("kappa[8,1] (Salps):", model$kappa[8, 1], "  Expected: 0.5\n")
cat("kappa[10,1] (Fish_Small):", model$kappa[10, 1], "  Expected: 0.7 (placeholder)\n")
cat("max(metab_cost):", max(model$metab_cost), "  Expected: 0 (Phase 1)\n")
cat("max(starv_mort):", max(model$starv_mort), "  Expected: 0\n")

# Step 4: Run model with MetabConst=0 (current default)
cat("\n--- Test 4: Full model run (10 years) ---\n")
t0 <- Sys.time()
mdl <- zoomss_model(input_params = env_data, Groups = Groups, isave = 2)
t1 <- Sys.time()
cat("Model completed in", round(difftime(t1, t0, units = "secs"), 1), "seconds\n")

# Check outputs are non-zero
N_eq <- averageTimeSeries(mdl, var = "abundance", n_years = 5)
gg_eq <- averageTimeSeries(mdl, var = "growth", n_years = 5)
cat("\nEquilibrium abundance per group:\n")
for (i in 1:nrow(N_eq)) {
  cat(sprintf("  %s: %.3e (max growth: %.3e)\n",
              Groups$Species[i], sum(N_eq[i, ]), max(gg_eq[i, ])))
}

# Step 5: Verify growth rates are reasonable
cat("\n--- Test 5: Growth rate sanity checks ---\n")
cat("All groups have positive growth:", all(apply(gg_eq, 1, max) > 0), "\n")
cat("No NaN in growth:", !any(is.nan(gg_eq)), "\n")
cat("No Inf in growth:", !any(is.infinite(gg_eq)), "\n")

# Step 6: Run with Kappa=1 to compare against old model
cat("\n--- Test 6: Run with Kappa=1, MetabConst=0 (comparable to old model) ---\n")
Groups_k1 <- Groups
Groups_k1$Kappa[Groups_k1$Type == "Zooplankton"] <- 1.0
Groups_k1$Kappa[Groups_k1$Type == "Fish"] <- NA  # will be set to 0.7 internally

mdl_k1 <- zoomss_model(
  input_params = createInputParams(time = seq(0, 10, by = 0.1), sst = 20, chl = 0.5),
  Groups = Groups_k1,
  isave = 2
)
N_k1 <- averageTimeSeries(mdl_k1, var = "abundance", n_years = 5)
cat("With Kappa=1, equilibrium total abundance per group:\n")
for (i in 1:nrow(N_k1)) {
  cat(sprintf("  %s: %.3e\n", Groups$Species[i], sum(N_k1[i, ])))
}

cat("\n=== Phase 1 Verification COMPLETE ===\n")
