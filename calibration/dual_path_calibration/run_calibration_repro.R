#!/usr/bin/env Rscript
# =============================================================================
# Run Fish Reproduction Calibration
# =============================================================================
#
# Standalone script to run the calibration pipeline.
# Designed for both local execution and HPC submission.
#
# Usage:
#   Rscript vignettes/run_calibration_repro.R                    # defaults
#   Rscript vignettes/run_calibration_repro.R --n_samples 1000   # more samples
#   Rscript vignettes/run_calibration_repro.R --phase benchmark  # benchmark only
#
# Requires: zoomss package (devtools::load_all()), lhs, future, furrr
# =============================================================================

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default, type = "character") {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  val <- args[idx + 1]
  switch(type,
    numeric  = as.numeric(val),
    integer  = as.integer(val),
    logical  = as.logical(val),
    val
  )
}

# Configuration
n_samples        <- parse_arg("--n_samples",        500L,    "integer")
n_workers        <- parse_arg("--n_workers",         14L,    "integer")
sst              <- parse_arg("--sst",               15,     "numeric")
screening_years  <- parse_arg("--screening_years",   100L,   "integer")
refinement_years <- parse_arg("--refinement_years",  200L,   "integer")
top_n_refine     <- parse_arg("--top_n_refine",       5L,    "integer")
seed             <- parse_arg("--seed",               42L,   "integer")
cache_dir        <- parse_arg("--cache_dir",         "calibration_repro_cache")
phase            <- parse_arg("--phase",             "all")   # all, benchmark, lhs, refine

cat("=============================================================\n")
cat("ZooMSS Fish Reproduction Calibration\n")
cat("=============================================================\n")
cat(sprintf("  n_samples:        %d\n", n_samples))
cat(sprintf("  n_workers:        %d\n", n_workers))
cat(sprintf("  sst:              %.1f\n", sst))
cat(sprintf("  screening_years:  %d\n", screening_years))
cat(sprintf("  refinement_years: %d\n", refinement_years))
cat(sprintf("  top_n_refine:     %d\n", top_n_refine))
cat(sprintf("  seed:             %d\n", seed))
cat(sprintf("  cache_dir:        %s\n", cache_dir))
cat(sprintf("  phase:            %s\n", phase))
cat("=============================================================\n\n")

# --- Load package ---
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(zoomss)
}

# Check dependencies
required_pkgs <- c("lhs", "future", "furrr")
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c('",
       paste(missing, collapse = "', '"), "'))")
}

# --- Setup ---
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

log10_chl_seq <- seq(-1.7, 0.5, by = 0.1)
chl_levels <- 10^log10_chl_seq

t_start <- Sys.time()

# =============================================================================
# Phase 1: Benchmark
# =============================================================================
if (phase %in% c("all", "benchmark")) {
  cat("\n=== Phase 1: Generating Legacy Benchmark ===\n")

  benchmark <- generate_legacy_benchmark(
    chl_levels = chl_levels,
    sst = sst,
    n_workers = n_workers,
    cache_dir = file.path(cache_dir, "benchmark")
  )

  cat(sprintf("Benchmark complete: %d chl levels, SST = %.0f\n",
              length(chl_levels), sst))
  cat(sprintf("Zoo proportions range: [%.4f, %.4f]\n",
              min(benchmark$zoo_proportions, na.rm = TRUE),
              max(benchmark$zoo_proportions, na.rm = TRUE)))

  if (phase == "benchmark") {
    cat("\nBenchmark-only run complete.\n")
    cat(sprintf("Total time: %.1f minutes\n", difftime(Sys.time(), t_start, units = "mins")))
    quit(save = "no", status = 0)
  }
}

# =============================================================================
# Phase 2: LHS Exploration
# =============================================================================
if (phase %in% c("all", "lhs")) {
  cat("\n=== Phase 2: LHS Parameter Exploration ===\n")

  # Load benchmark if running lhs-only
  if (phase == "lhs") {
    benchmark_file <- file.path(cache_dir, "benchmark",
                                sprintf("benchmark_sst%.0f.rds", sst))
    if (!file.exists(benchmark_file)) {
      stop("Benchmark not found. Run with --phase benchmark first.")
    }
    benchmark <- readRDS(benchmark_file)
  }

  # Generate samples
  cat(sprintf("Generating %d LHS samples (14 dimensions)...\n", n_samples))
  lhs_samples <- generate_lhs_samples(n_samples = n_samples, seed = seed)

  # Verify energy constraint
  n_valid <- sum(sapply(seq_len(nrow(lhs_samples)), function(i) {
    check_energy_constraint(lhs_samples[i, ])
  }))
  cat(sprintf("Energy constraint check: %d/%d valid\n", n_valid, nrow(lhs_samples)))

  # Save LHS design for reproducibility
  saveRDS(lhs_samples, file.path(cache_dir, "lhs_design.rds"))

  # Run exploration
  lhs_results <- run_lhs_exploration(
    lhs_samples = lhs_samples,
    benchmark = benchmark,
    n_years = screening_years,
    n_workers = n_workers,
    cache_dir = file.path(cache_dir, "lhs")
  )

  # Summary statistics
  cat("\n--- LHS Results Summary ---\n")
  cat(sprintf("Total evaluated:    %d\n", nrow(lhs_results)))
  cat(sprintf("Score range:        [%.4f, %.4f]\n",
              min(lhs_results$score), max(lhs_results$score)))
  cat(sprintf("Full coexistence:   %d (%.1f%%)\n",
              sum(lhs_results$coexistence <= 0.01),
              100 * mean(lhs_results$coexistence <= 0.01)))
  cat(sprintf("Good zoo comp:      %d (%.1f%%)\n",
              sum(lhs_results$zoo_comp <= 0.3),
              100 * mean(lhs_results$zoo_comp <= 0.3)))

  # Filter candidates
  cat("\n=== Phase 2b: Filtering Candidates ===\n")
  candidates <- filter_lhs_candidates(lhs_results, top_n = top_n_refine * 2)

  if (nrow(candidates) == 0) {
    cat("\nWARNING: No candidates passed filters!\n")
    cat("Consider relaxing constraints or increasing n_samples.\n")

    # Relaxed filter
    candidates_relaxed <- filter_lhs_candidates(
      lhs_results,
      max_coexistence = 0.1,
      max_zoo_comp = 0.5,
      top_n = top_n_refine * 2
    )
    cat(sprintf("Relaxed filter: %d candidates\n", nrow(candidates_relaxed)))
    candidates <- candidates_relaxed
  }

  saveRDS(list(lhs_results = lhs_results, candidates = candidates),
          file.path(cache_dir, "lhs_summary.rds"))

  if (phase == "lhs") {
    cat("\nLHS-only run complete.\n")
    cat(sprintf("Total time: %.1f minutes\n",
                difftime(Sys.time(), t_start, units = "mins")))
    quit(save = "no", status = 0)
  }
}

# =============================================================================
# Phase 3: Refinement
# =============================================================================
if (phase %in% c("all", "refine")) {
  cat("\n=== Phase 3: Refining Top Candidates ===\n")

  # Load if running refine-only
  if (phase == "refine") {
    benchmark_file <- file.path(cache_dir, "benchmark",
                                sprintf("benchmark_sst%.0f.rds", sst))
    summary_file <- file.path(cache_dir, "lhs_summary.rds")
    if (!file.exists(benchmark_file) || !file.exists(summary_file)) {
      stop("Benchmark or LHS results not found. Run earlier phases first.")
    }
    benchmark <- readRDS(benchmark_file)
    lhs_data <- readRDS(summary_file)
    candidates <- lhs_data$candidates
    lhs_results <- lhs_data$lhs_results
    lhs_samples <- readRDS(file.path(cache_dir, "lhs_design.rds"))
  }

  n_refine <- min(top_n_refine, nrow(candidates))
  param_names <- names(generate_lhs_samples(1))  # Get param names
  refined <- list()

  for (k in seq_len(n_refine)) {
    cat(sprintf("\n--- Refining candidate %d/%d (LHS score: %.4f) ---\n",
                k, n_refine, candidates$score[k]))

    par_init <- as.numeric(candidates[k, param_names])
    names(par_init) <- param_names

    refined[[k]] <- refine_candidate(
      par_init = par_init,
      benchmark = benchmark,
      n_years = refinement_years
    )

    cat(sprintf("  Refined score: %.4f (convergence: %d)\n",
                refined[[k]]$score, refined[[k]]$convergence))
  }

  # Select best
  best_idx <- which.min(sapply(refined, function(x) x$score))
  best <- refined[[best_idx]]

  cat("\n=== Best Calibrated Parameters ===\n")
  for (nm in names(best$par)) {
    cat(sprintf("  %-12s = %.4f\n", nm, best$par[nm]))
  }
  cat(sprintf("\n  R_frac = %.4f\n", 1 - best$par["f_M"] - best$par["K_growth"]))
  cat(sprintf("  Best score: %.4f\n", best$score))

  # Detailed metric breakdown
  cat("\n  Metric breakdown:\n")
  for (nm in names(best$details$metric_scores)) {
    cat(sprintf("    %-15s = %.4f\n", nm, best$details$metric_scores[nm]))
  }

  # Save final results
  calibration <- list(
    benchmark = benchmark,
    lhs_results = lhs_results,
    candidates = candidates,
    refined = refined,
    best = best,
    param_space = repro_param_space(),
    settings = list(
      n_samples = n_samples, sst = sst,
      screening_years = screening_years,
      refinement_years = refinement_years, seed = seed
    )
  )

  results_file <- file.path(cache_dir, "calibration_results.rds")
  saveRDS(calibration, results_file)
  cat(sprintf("\nResults saved to: %s\n", results_file))
}

# --- Final timing ---
t_end <- Sys.time()
cat(sprintf("\n=== Total elapsed time: %.1f minutes ===\n",
            difftime(t_end, t_start, units = "mins")))
