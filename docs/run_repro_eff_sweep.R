# =============================================================================
# run_repro_eff_sweep.R
#
# Reproductive efficiency × fishing yield sweep for the dual-pathway ZooMSS.
# Explores fish co-existence under different repro_eff and fishing scenarios.
#
# Usage:
#   Rscript run_repro_eff_sweep.R
#   # or source("run_repro_eff_sweep.R") from an interactive R session
#
# Output files (saved to cache_dir):
#   results_uniform_repro.rds   — Part 1: uniform repro_eff × effort sweep
#   results_group_repro.rds     — Part 2: group-specific repro_eff × effort
#   results_selective_fishing.rds — Part 3: selective fishing × repro_eff
#   sweep_config.rds            — configuration used (for .Rmd validation)
# =============================================================================

devtools::load_all()
library(parallel)

pkg_path <- here::here()
cat("Package path:", pkg_path, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

cache_dir <- here::here("vignettes", "cache", "repro_eff_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- max(1, detectCores() - 2)
cat("Using", n_cores, "cores\n")

# Force single-threaded BLAS to avoid contention across parallel workers
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)

# Shared simulation parameters
n_years   <- 300
dt_step   <- 0.1
sim_time  <- seq(0, n_years, by = dt_step)
sst_const <- 15
chl_const <- 0.3
isave     <- 2

# Effort sweep — finer resolution at low effort where yield curves peak
effort_levels <- c(0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3,
                   0.4, 0.5, 0.75, 1.0, 1.5, 2.0)

# Default catchability
q_fixed <- 1.0

# Base groups
Groups_base <- getGroups()
fish_rows   <- which(Groups_base$Type == "Fish")
fish_names  <- Groups_base$Species[fish_rows]
num_fish    <- length(fish_rows)

# Dynamic effort column names (dual-pathway convention)
effort_col_names <- paste0("effort_", tolower(gsub(" ", "_", fish_names)))
names(effort_col_names) <- fish_names
cat("Fish groups:", paste(fish_names, collapse = ", "), "\n")
cat("Effort columns:", paste(effort_col_names, collapse = ", "), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build input_params with effort columns using dual-pathway naming
# ─────────────────────────────────────────────────────────────────────────────

build_env <- function(sim_time, sst, chl, effort_vec, effort_col_names) {
  # effort_vec: named or positional vector of effort per fish group
  env <- createInputParams(time = sim_time, sst = sst, chl = chl)

  # Add effort columns with dual-pathway naming convention
  for (i in seq_along(effort_col_names)) {
    env[[effort_col_names[i]]] <- rep(effort_vec[i], nrow(env))
  }
  return(env)
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: run a single model for one repro_eff vector + effort vector
# ─────────────────────────────────────────────────────────────────────────────

run_single <- function(repro_eff_vec, effort_vec, Groups_base, fish_rows,
                       q_vec, sim_time, sst_const, chl_const, isave,
                       effort_col_names) {

  Groups_mod <- Groups_base
  Groups_mod$repro_eff[fish_rows] <- repro_eff_vec

  # Set catchability (q) per fish group
  Groups_mod <- setFishingParams(Groups_mod,
                                 q_small = q_vec[1],
                                 q_med   = q_vec[2],
                                 q_large = q_vec[3])

  env <- build_env(sim_time, sst_const, chl_const, effort_vec, effort_col_names)

  zoomss_model(input_params = env, Groups = Groups_mod, isave = isave)
}


# =============================================================================
# PART 1: Uniform repro_eff × uniform effort sweep
# =============================================================================
# All fish get the same repro_eff, all fish get the same effort.
# This gives yield curves under different reproductive efficiency regimes.

repro_eff_values <- c(1, 0.5, 0.1)

uniform_jobs <- expand.grid(
  re_idx  = seq_along(repro_eff_values),
  eff_idx = seq_along(effort_levels)
)

cat("\n=== Part 1: Uniform repro_eff × effort yield curves ===\n")
cat("repro_eff values:", paste(repro_eff_values, collapse = ", "), "\n")
cat("effort levels:", length(effort_levels), "values from",
    min(effort_levels), "to", max(effort_levels), "\n")
cat("Running", nrow(uniform_jobs), "model combinations across", n_cores, "cores\n")

t1 <- Sys.time()

cl <- makeCluster(n_cores)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

clusterExport(cl, c("run_single", "build_env", "uniform_jobs",
                     "repro_eff_values", "effort_levels",
                     "Groups_base", "fish_rows", "num_fish",
                     "q_fixed", "sim_time", "sst_const", "chl_const", "isave",
                     "effort_col_names", "pkg_path"))
clusterEvalQ(cl, {
  devtools::load_all(pkg_path, quiet = TRUE)
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
})

uniform_results_flat <- parLapply(cl, seq_len(nrow(uniform_jobs)), function(j) {
  re_idx  <- uniform_jobs$re_idx[j]
  eff_idx <- uniform_jobs$eff_idx[j]

  repro_vec  <- rep(repro_eff_values[re_idx], num_fish)
  effort_vec <- rep(effort_levels[eff_idx], num_fish)
  q_vec      <- rep(q_fixed, num_fish)

  run_single(repro_vec, effort_vec, Groups_base, fish_rows, q_vec,
             sim_time, sst_const, chl_const, isave, effort_col_names)
})

# Reshape to nested list: results[[re_idx]][[eff_idx]]
results_repro <- vector("list", length(repro_eff_values))
for (i in seq_along(results_repro)) {
  results_repro[[i]] <- vector("list", length(effort_levels))
}

for (j in seq_len(nrow(uniform_jobs))) {
  re_idx  <- uniform_jobs$re_idx[j]
  eff_idx <- uniform_jobs$eff_idx[j]
  results_repro[[re_idx]][[eff_idx]] <- uniform_results_flat[[j]]
}

t2 <- Sys.time()
cat("Part 1 completed in", round(difftime(t2, t1, units = "mins"), 1), "minutes\n")

saveRDS(results_repro, file = file.path(cache_dir, "results_uniform_repro.rds"))
cat("Saved:", file.path(cache_dir, "results_uniform_repro.rds"), "\n")

rm(uniform_results_flat)
gc()


# =============================================================================
# PART 2: Group-specific repro_eff × uniform effort
# =============================================================================
# Key hypothesis: co-existence may require size-dependent repro_eff because
# larger fish produce more eggs per unit biomass but have lower survival,
# while smaller fish compete more directly with zooplankton.
#
# Scenarios explore gradients, inversions, and extreme contrasts.

repro_scenarios <- list(
  # --- Baselines ---
  "Uniform 1.0"            = c(1,      1,      1),
  "Uniform 0.1"            = c(0.1,    0.1,    0.1),
  "Uniform 0.01"           = c(0.01,   0.01,   0.01),

  # --- Size-dependent gradients (large fish lower) ---
  # Rationale: larger fish have higher fecundity but lower larval survival
  "Gradient 10x"           = c(1,      0.1,    0.01),
  "Gradient 100x"          = c(1,      0.01,   0.0001),
  "Mild gradient"          = c(1,      0.5,    0.1),
  "Steep gradient"         = c(1,      0.05,   0.001)

  # # --- Inverted gradients (small fish lower) ---
  # # Tests whether small fish are the ones being competitively excluded
  # "Inverted 10x"           = c(0.01,   0.1,    1),
  # "Inverted mild"          = c(0.1,    0.5,    1),
  #
  # # --- Single group dominant ---
  # # Identifies which groups can persist alone
  # "Small only"             = c(1,      0.001,  0.001),
  # "Med only"               = c(0.001,  1,      0.001),
  # "Large only"             = c(0.001,  0.001,  1),
  #
  # # --- Pairwise co-existence ---
  # "Small+Med"              = c(1,      1,      0.001),
  # "Small+Large"            = c(1,      0.001,  1),
  # "Med+Large"              = c(0.001,  1,      1)
)

scenario_names <- names(repro_scenarios)

group_jobs <- expand.grid(
  sc_idx  = seq_along(repro_scenarios),
  eff_idx = seq_along(effort_levels)
)

cat("\n=== Part 2: Group-specific repro_eff × effort ===\n")
cat("Scenarios:", length(repro_scenarios), "\n")
for (i in seq_along(repro_scenarios)) {
  cat(sprintf("  %-25s [%s]\n", scenario_names[i],
              paste(repro_scenarios[[i]], collapse = ", ")))
}
cat("Running", nrow(group_jobs), "model combinations across", n_cores, "cores\n")

t3 <- Sys.time()

clusterExport(cl, c("group_jobs", "repro_scenarios"))

group_results_flat <- parLapply(cl, seq_len(nrow(group_jobs)), function(j) {
  sc_idx  <- group_jobs$sc_idx[j]
  eff_idx <- group_jobs$eff_idx[j]

  repro_vec  <- repro_scenarios[[sc_idx]]
  effort_vec <- rep(effort_levels[eff_idx], num_fish)
  q_vec      <- rep(q_fixed, num_fish)

  run_single(repro_vec, effort_vec, Groups_base, fish_rows, q_vec,
             sim_time, sst_const, chl_const, isave, effort_col_names)
})

# Reshape to nested list: results[[sc_idx]][[eff_idx]]
results_group_repro <- vector("list", length(repro_scenarios))
for (i in seq_along(results_group_repro)) {
  results_group_repro[[i]] <- vector("list", length(effort_levels))
}

for (j in seq_len(nrow(group_jobs))) {
  sc_idx  <- group_jobs$sc_idx[j]
  eff_idx <- group_jobs$eff_idx[j]
  results_group_repro[[sc_idx]][[eff_idx]] <- group_results_flat[[j]]
}

t4 <- Sys.time()
cat("Part 2 completed in", round(difftime(t4, t3, units = "mins"), 1), "minutes\n")

saveRDS(results_group_repro, file = file.path(cache_dir, "results_group_repro.rds"))
cat("Saved:", file.path(cache_dir, "results_group_repro.rds"), "\n")

rm(group_results_flat)
gc()


# =============================================================================
# PART 3: Selective fishing — different effort per group
# =============================================================================
# Explores whether selective fishing pressure can promote co-existence.
# Uses repro_eff = 1 (baseline) with varying effort patterns.
# Each "fishing strategy" defines relative effort per group.

fishing_strategies <- list(
  # Name                     = c(Small, Med, Large)
  "Uniform"                  = c(1,    1,    1),
  # "Small only"               = c(1,    0,    0),
  # "Med only"                 = c(0,    1,    0),
  # "Large only"               = c(0,    0,    1),
  # "Balanced (no small)"      = c(0,    1,    1),
  # "Balanced (no large)"      = c(1,    1,    0),
  "Progressive (more large)" = c(0.25, 0.5,  1),
  "Inverse (more small)"     = c(1,    0.5,  0.25)
)

# Effort multiplier levels (applied to the strategy weights)
effort_multipliers <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0)

# repro_eff scenarios to combine with selective fishing
selective_repro_effs <- list(
  "Uniform 1.0"     = c(1, 1, 1),
  "Gradient 10x"    = c(1, 0.1, 0.01),
  "Mild gradient"   = c(1, 0.5, 0.1)
)

selective_jobs <- expand.grid(
  strat_idx  = seq_along(fishing_strategies),
  mult_idx   = seq_along(effort_multipliers),
  re_idx     = seq_along(selective_repro_effs)
)

cat("\n=== Part 3: Selective fishing × repro_eff ===\n")
cat("Fishing strategies:", length(fishing_strategies), "\n")
cat("Effort multipliers:", length(effort_multipliers), "\n")
cat("repro_eff scenarios:", length(selective_repro_effs), "\n")
cat("Running", nrow(selective_jobs), "model combinations across", n_cores, "cores\n")

t5 <- Sys.time()

clusterExport(cl, c("selective_jobs", "fishing_strategies",
                     "effort_multipliers", "selective_repro_effs"))

selective_results_flat <- parLapply(cl, seq_len(nrow(selective_jobs)), function(j) {
  strat_idx <- selective_jobs$strat_idx[j]
  mult_idx  <- selective_jobs$mult_idx[j]
  re_idx    <- selective_jobs$re_idx[j]

  strategy   <- fishing_strategies[[strat_idx]]
  multiplier <- effort_multipliers[mult_idx]
  effort_vec <- strategy * multiplier
  repro_vec  <- selective_repro_effs[[re_idx]]
  q_vec      <- rep(q_fixed, num_fish)

  run_single(repro_vec, effort_vec, Groups_base, fish_rows, q_vec,
             sim_time, sst_const, chl_const, isave, effort_col_names)
})

# Reshape to nested list: results[[strat_idx]][[mult_idx]][[re_idx]]
results_selective <- vector("list", length(fishing_strategies))
for (i in seq_along(results_selective)) {
  results_selective[[i]] <- vector("list", length(effort_multipliers))
  for (k in seq_along(results_selective[[i]])) {
    results_selective[[i]][[k]] <- vector("list", length(selective_repro_effs))
  }
}

for (j in seq_len(nrow(selective_jobs))) {
  strat_idx <- selective_jobs$strat_idx[j]
  mult_idx  <- selective_jobs$mult_idx[j]
  re_idx    <- selective_jobs$re_idx[j]
  results_selective[[strat_idx]][[mult_idx]][[re_idx]] <- selective_results_flat[[j]]
}

t6 <- Sys.time()
cat("Part 3 completed in", round(difftime(t6, t5, units = "mins"), 1), "minutes\n")

saveRDS(results_selective, file = file.path(cache_dir, "results_selective_fishing.rds"))
cat("Saved:", file.path(cache_dir, "results_selective_fishing.rds"), "\n")

rm(selective_results_flat)
gc()


# =============================================================================
# Shut down cluster
# =============================================================================
stopCluster(cl)


# =============================================================================
# Save configuration for .Rmd validation
# =============================================================================

sweep_config <- list(
  # Simulation
  n_years          = n_years,
  dt_step          = dt_step,
  sim_time_range   = range(sim_time),
  sst_const        = sst_const,
  chl_const        = chl_const,
  isave            = isave,

  # Fish
  fish_names       = fish_names,
  fish_rows        = fish_rows,
  num_fish         = num_fish,
  effort_col_names = effort_col_names,
  q_fixed          = q_fixed,

  # Part 1
  effort_levels    = effort_levels,
  repro_eff_values = repro_eff_values,

  # Part 2
  repro_scenarios  = repro_scenarios,
  scenario_names   = scenario_names,

  # Part 3
  fishing_strategies    = fishing_strategies,
  effort_multipliers    = effort_multipliers,
  selective_repro_effs  = selective_repro_effs,

  # Meta
  n_cores_used     = n_cores,
  timestamp        = Sys.time(),
  r_version        = R.version.string
)

saveRDS(sweep_config, file = file.path(cache_dir, "sweep_config.rds"))
cat("\nSaved:", file.path(cache_dir, "sweep_config.rds"), "\n")


# =============================================================================
# Summary
# =============================================================================

cat("\n==================================================\n")
cat("All sweeps complete\n")
cat("==================================================\n")
cat("Part 1 (uniform):   ", length(repro_eff_values), "× ",
    length(effort_levels), "=", length(repro_eff_values) * length(effort_levels),
    "models\n")
cat("Part 2 (group):     ", length(repro_scenarios), "× ",
    length(effort_levels), "=", length(repro_scenarios) * length(effort_levels),
    "models\n")
cat("Part 3 (selective): ", length(fishing_strategies), "× ",
    length(effort_multipliers), "× ", length(selective_repro_effs), "=",
    nrow(selective_jobs), "models\n")
cat("Total models:       ",
    nrow(uniform_jobs) + nrow(group_jobs) + nrow(selective_jobs), "\n")
cat("Total time:         ", round(difftime(t6, t1, units = "mins"), 1), "minutes\n")
cat("Cache directory:    ", normalizePath(cache_dir), "\n")
cat("Files:\n")
cat(paste0("  ", list.files(cache_dir, full.names = FALSE)), sep = "\n")
cat("\n")
