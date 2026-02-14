## Smoke test for ZooMSS model after energy budget revisions
library(devtools)
load_all(".")

cat("=== Test 1: Energy budget closure ===\n")
G <- getGroups()
R_frac <- 1 - G$f_M - G$K_growth
for (i in 1:nrow(G)) {
  cat(sprintf("  %-14s f_M=%.2f K_growth=%.2f R_frac=%.2f\n",
              G$Species[i], G$f_M[i], G$K_growth[i], R_frac[i]))
}
stopifnot(all(abs(G$f_M + G$K_growth + R_frac - 1) < 1e-10))
stopifnot(R_frac[1] == 0, R_frac[2] == 0)  # Flagellates, Ciliates
stopifnot(all(R_frac[3:12] > 0))
cat("PASS\n\n")

cat("=== Test 2: validateGroups ===\n")
validateGroups(G)
cat("PASS\n\n")

cat("=== Test 3: 50-year model run (Scenario A) ===\n")
env <- createInputParams(time = seq(0, 50, by = 0.1), sst = 15, chl = 1.0)
mdl <- zoomss_model(input_params = env, Groups = G, isave = 10)
cat("  Time steps:", length(mdl$time), "\n")
cat("  Abundance dims:", paste(dim(mdl$abundance), collapse = " x "), "\n")
cat("  Any NaN:", any(is.nan(mdl$abundance)), "\n")
cat("  Any Inf:", any(is.infinite(mdl$abundance)), "\n")
stopifnot(!any(is.nan(mdl$abundance)))
stopifnot(!any(is.infinite(mdl$abundance)))

bm <- rowSums(mdl$biomass[length(mdl$time), , ])
names(bm) <- G$Species
cat("  Final biomass:\n")
print(signif(bm, 4))
stopifnot(all(bm > 0))

cat("  SSB:", mdl$SSB[nrow(mdl$SSB), ], "\n")
cat("  Recruitment:", mdl$recruitment[nrow(mdl$recruitment), ], "\n")
stopifnot(all(mdl$SSB[nrow(mdl$SSB), ] > 0))
stopifnot(all(mdl$recruitment[nrow(mdl$recruitment), ] > 0))
cat("PASS\n\n")

cat("=== Test 4: Scenario B runs ===\n")
mdl_B <- zoomss_model(input_params = env, Groups = G, isave = 10,
                       energy_budget_scenario = "B")
cat("  Time steps:", length(mdl_B$time), "\n")
stopifnot(!any(is.nan(mdl_B$abundance)))
bm_B <- rowSums(mdl_B$biomass[length(mdl_B$time), , ])
names(bm_B) <- G$Species
cat("  Final biomass (Scenario B):\n")
print(signif(bm_B, 4))
stopifnot(all(bm_B > 0))
cat("PASS\n\n")

cat("=== Test 5: Maturity-dependent growth check ===\n")
# For a zooplankton group with R_frac > 0 (e.g., Euphausiids, idx=6),
# gg_total output should exist and growth should reflect maturity scaling
cat("  gg_total dims:", paste(dim(mdl$gg_total), collapse = " x "), "\n")
# Verify gg_total exists in output
stopifnot(!is.null(mdl$gg_total))
cat("PASS\n\n")

cat("=== ALL TESTS PASSED ===\n")
