# =============================================================================
# run_repro_eff_sweep.R
#
# Runs the reproductive efficiency sensitivity analysis model sweeps in
# parallel and saves results as .rds files for use by the companion
# repro_eff_sensitivity.Rmd document.
#
# Usage:
#   Rscript run_repro_eff_sweep.R
#   # or source("run_repro_eff_sweep.R") from an interactive R session
#
# Output files (saved to cache_dir):
#   results_uniform_repro.rds   — Part 1: uniform repro_eff sweep
#   results_group_repro.rds     — Part 2: group-specific repro_eff sweep
#   sweep_config.rds            — configuration used (for .Rmd validation)
# =============================================================================

devtools::load_all()
library(parallel)

# Path to the zoomss package source (adjust if your project root differs)
pkg_path <- here::here()
cat("Package path:", pkg_path, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

# Output directory for cached results (in vignettes/ alongside the .Rmd)
cache_dir <- here::here("vignettes", "repro_eff_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Number of cores — leave one free for system stability
n_cores <- max(1, detectCores() - 1)
cat("Using", n_cores, "cores\n")

# Force single-threaded BLAS to avoid contention across parallel workers
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)

# Shared parameters
effort_levels    <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0, 1.5, 2.0, 3.0)
q_fixed          <- 1.0
sim_time         <- seq(0, 400, by = 0.1)
sst_const        <- 15
chl_const        <- 1.0
isave            <- 5

Groups_base <- getGroups()
fish_rows   <- which(Groups_base$Type == "Fish")
fish_names  <- c("Fish_Small", "Fish_Med", "Fish_Large")

# ─────────────────────────────────────────────────────────────────────────────
# Helper: run a single model for one repro_eff vector + one effort level
# ─────────────────────────────────────────────────────────────────────────────

run_single <- function(repro_eff_vec, effort_val, Groups_base, fish_rows,
                       q_fixed, sim_time, sst_const, chl_const, isave) {

  Groups_mod <- Groups_base
  Groups_mod$repro_eff[fish_rows] <- repro_eff_vec

  Groups_mod <- setFishingParams(Groups_mod,
                                 q_small = q_fixed,
                                 q_med   = q_fixed,
                                 q_large = q_fixed)

  env <- createInputParams(time = sim_time, sst = sst_const, chl = chl_const)
  env <- addFishingEffort(env,
                          effort_small = effort_val,
                          effort_med   = effort_val,
                          effort_large = effort_val)

  zoomss_model(input_params = env, Groups = Groups_mod, isave = isave)
}


# =============================================================================
# PART 1: Uniform repro_eff sweep
# =============================================================================

repro_eff_values <- c(1, 0.1, 0.01)

# Build job list: all combinations of repro_eff × effort
uniform_jobs <- expand.grid(
  re_idx  = seq_along(repro_eff_values),
  eff_idx = seq_along(effort_levels)
)

cat("\n=== Part 1: Uniform repro_eff ===\n")
cat("Running", nrow(uniform_jobs), "model combinations across", n_cores, "cores\n")

t1 <- Sys.time()

# Create PSOCK cluster (works on Windows, macOS, and Linux)
cl <- makeCluster(n_cores)
on.exit(stopCluster(cl), add = TRUE)

# Export objects and load zoomss on each worker
clusterExport(cl, c("run_single", "uniform_jobs", "repro_eff_values",
                     "effort_levels", "Groups_base", "fish_rows",
                     "q_fixed", "sim_time", "sst_const", "chl_const", "isave",
                     "pkg_path"))
clusterEvalQ(cl, {
  devtools::load_all(pkg_path)
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
})

uniform_results_flat <- parLapply(cl, seq_len(nrow(uniform_jobs)), function(j) {
  re_idx  <- uniform_jobs$re_idx[j]
  eff_idx <- uniform_jobs$eff_idx[j]

  repro_vec <- rep(repro_eff_values[re_idx], length(fish_rows))

  run_single(repro_vec, effort_levels[eff_idx],
             Groups_base, fish_rows, q_fixed,
             sim_time, sst_const, chl_const, isave)
})

# Reshape flat list back to nested list: results_repro[[re_idx]][[eff_idx]]
results_repro <- vector("list", length(repro_eff_values))
for (i in seq_along(results_repro)) results_repro[[i]] <- vector("list", length(effort_levels))

for (j in seq_len(nrow(uniform_jobs))) {
  re_idx  <- uniform_jobs$re_idx[j]
  eff_idx <- uniform_jobs$eff_idx[j]
  results_repro[[re_idx]][[eff_idx]] <- uniform_results_flat[[j]]
}

t2 <- Sys.time()
cat("Part 1 completed in", round(difftime(t2, t1, units = "mins"), 1), "minutes\n")

saveRDS(results_repro, file = file.path(cache_dir, "results_uniform_repro_v2.rds"))
cat("Saved:", file.path(cache_dir, "results_uniform_repro_v2.rds"), "\n")

rm(uniform_results_flat)  # Free memory before Part 2


# =============================================================================
# PART 2: Group-specific repro_eff sweep
# =============================================================================

# repro_scenarios <- list(
#   "Uniform default"     = c(1e-3,  1e-3,  1e-3),
#   "Uniform low"         = c(1e-5,  1e-5,  1e-5),
#   "Gradient 10x"        = c(1e-3,  1e-4,  1e-5),
#   "Gradient 10x (L=1e-8)"  = c(1e-3,  1e-4,  1e-8),
#   "Gradient 100x"       = c(1e-2,  1e-4,  1e-6),
#   "Gradient 100x (L=1e-8)" = c(1e-2,  1e-4,  1e-8),
#   "Large fish fragile"  = c(1e-3,  1e-3,  1e-6),
#   "Small fish dominant" = c(1e-2,  1e-5,  1e-5)
# )

repro_scenarios <- list(
  "Large fish fragile"  = c(0.01,  0.01,  0.001),
  "Small fish dominant" = c(0.1,  0.01,  0.01)
)

scenario_names <- names(repro_scenarios)

# Build job list
group_jobs <- expand.grid(
  sc_idx  = seq_along(repro_scenarios),
  eff_idx = seq_along(effort_levels)
)

cat("\n=== Part 2: Group-specific repro_eff ===\n")
cat("Running", nrow(group_jobs), "model combinations across", n_cores, "cores\n")

t3 <- Sys.time()

# Export additional Part 2 objects to existing cluster
clusterExport(cl, c("group_jobs", "repro_scenarios"))

group_results_flat <- parLapply(cl, seq_len(nrow(group_jobs)), function(j) {
  sc_idx  <- group_jobs$sc_idx[j]
  eff_idx <- group_jobs$eff_idx[j]

  repro_vec <- repro_scenarios[[sc_idx]]

  run_single(repro_vec, effort_levels[eff_idx],
             Groups_base, fish_rows, q_fixed,
             sim_time, sst_const, chl_const, isave)
})

# Reshape to nested list: results_group_repro[[sc_idx]][[eff_idx]]
results_group_repro <- vector("list", length(repro_scenarios))
for (i in seq_along(results_group_repro)) results_group_repro[[i]] <- vector("list", length(effort_levels))

for (j in seq_len(nrow(group_jobs))) {
  sc_idx  <- group_jobs$sc_idx[j]
  eff_idx <- group_jobs$eff_idx[j]
  results_group_repro[[sc_idx]][[eff_idx]] <- group_results_flat[[j]]
}

t4 <- Sys.time()
cat("Part 2 completed in", round(difftime(t4, t3, units = "mins"), 1), "minutes\n")

saveRDS(results_group_repro, file = file.path(cache_dir, "results_group_repro_v2.rds"))
cat("Saved:", file.path(cache_dir, "results_group_repro_v2.rds"), "\n")

# Shut down the cluster
stopCluster(cl)


# =============================================================================
# Save configuration for .Rmd validation
# =============================================================================

sweep_config <- list(
  effort_levels    = effort_levels,
  repro_eff_values = repro_eff_values,
  repro_scenarios  = repro_scenarios,
  scenario_names   = scenario_names,
  q_fixed          = q_fixed,
  fish_names       = fish_names,
  fish_rows        = fish_rows,
  sim_time_range   = range(sim_time),
  sst_const        = sst_const,
  chl_const        = chl_const,
  isave            = isave,
  n_cores_used     = n_cores,
  timestamp        = Sys.time(),
  r_version        = R.version.string
)

saveRDS(sweep_config, file = file.path(cache_dir, "sweep_config.rds"))
cat("\nSaved:", file.path(cache_dir, "sweep_config.rds"), "\n")

cat("\n=== All sweeps complete ===\n")
cat("Total time:", round(difftime(t4, t1, units = "mins"), 1), "minutes\n")
cat("Cache directory:", normalizePath(cache_dir), "\n")
cat("Files:\n")
cat("  ", list.files(cache_dir, full.names = TRUE), sep = "\n  ")
cat("\n\nYou can now knit repro_eff_sensitivity.Rmd\n")
