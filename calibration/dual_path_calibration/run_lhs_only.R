#!/usr/bin/env Rscript
# =============================================================================
# Run LHS Exploration Only (Phases 1-2)
# =============================================================================
#
# Runs benchmark generation + LHS exploration without refinement.
# All simulations use seasonal environmental forcing via createEnviroData().
#
# Run lengths:
#   Benchmark:  400 years (seasonal SST/chl), assess final 100 years
#   LHS screen: 300 years (seasonal SST/chl), assess final 100 years
#
# Outputs:
#   calibration_repro_cache/benchmark/  — legacy benchmark runs
#   calibration_repro_cache/lhs/        — LHS results with checkpoints
#
# The output .rds files can then be transferred to a VM or picked up
# by the full pipeline or refinement step later.
#
# Usage:
#   Rscript run_lhs_only.R [n_samples] [n_workers]
#   Rscript run_lhs_only.R 500 14      # default
#   Rscript run_lhs_only.R 200 8       # lighter run
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
n_samples  <- if (length(args) >= 1) as.integer(args[1]) else 500L
n_workers  <- if (length(args) >= 2) as.integer(args[2]) else 14L

cat("=============================================================\n")
cat("ZooMSS Fish Reproduction Calibration — LHS Exploration Only\n")
cat(sprintf("  Samples:  %d\n", n_samples))
cat(sprintf("  Workers:  %d\n", n_workers))
cat(sprintf("  Forcing:  Seasonal (SST amp=4, chl amp=0.5)\n"))
cat(sprintf("  Benchmark: 400yr run, assess final 100yr\n"))
cat(sprintf("  LHS:       300yr run, assess final 100yr\n"))
cat(sprintf("  Started:  %s\n", Sys.time()))
cat("=============================================================\n\n")

# --- Setup ---
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(zoomss)
}

source(file.path("R", "zoomss_calibration_repro.R"))

cache_dir <- "calibration_repro_cache"
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

sst <- 15
seed <- 42L
sst_amplitude <- 4
chl_amplitude <- 0.5

# --- Phase 1: Benchmark ---
cat("=== Phase 1: Legacy Benchmark (400yr, seasonal) ===\n")
t1 <- Sys.time()

log10_chl_seq <- seq(-1.7, 0.5, by = 0.1)
chl_levels <- 10^log10_chl_seq
cat(sprintf("  %d chlorophyll levels (log10 chl: %.1f to %.1f)\n",
            length(chl_levels), min(log10_chl_seq), max(log10_chl_seq)))

benchmark <- generate_legacy_benchmark(
  chl_levels = chl_levels,
  sst = sst,
  n_years = 400,
  dt = 0.1,
  sst_amplitude = sst_amplitude,
  chl_amplitude = chl_amplitude,
  assess_years = 100,
  cache_dir = file.path(cache_dir, "benchmark"),
  n_workers = n_workers,
  force_rerun = FALSE
)

t1_elapsed <- difftime(Sys.time(), t1, units = "mins")
cat(sprintf("  Benchmark complete: %.1f minutes\n\n", as.numeric(t1_elapsed)))

# --- Phase 2: LHS ---
cat("=== Phase 2: LHS Parameter Exploration (300yr, seasonal) ===\n")
t2 <- Sys.time()

lhs_samples <- generate_lhs_samples(n_samples = n_samples, seed = seed)
cat(sprintf("  Generated %d LHS samples (14 dimensions)\n", nrow(lhs_samples)))

# 5 representative chl levels for screening
target_log10 <- c(-1.5, -1.0, -0.5, 0.0, 0.4)
chl_indices <- sapply(target_log10, function(t) {
  which.min(abs(log10(chl_levels) - t))
})
chl_indices <- unique(chl_indices)
cat(sprintf("  Evaluating at %d chl levels: %s\n",
            length(chl_indices),
            paste(sprintf("%.2f", log10(chl_levels[chl_indices])), collapse = ", ")))

# Estimate runtime (300yr seasonal runs are ~3x longer than 100yr constant)
est_per_sample_sec <- 5 * 120  # 5 chl levels × ~120s per 300yr seasonal run
est_total_min <- (n_samples / n_workers) * est_per_sample_sec / 60
cat(sprintf("  Estimated runtime: %.0f–%.0f hours (depends on hardware)\n",
            est_total_min * 0.5 / 60, est_total_min * 1.5 / 60))

lhs_results <- run_lhs_exploration(
  lhs_samples = lhs_samples,
  benchmark = benchmark,
  chl_indices = chl_indices,
  n_years = 300,
  n_workers = n_workers,
  cache_dir = file.path(cache_dir, "lhs"),
  batch_size = 50
)

t2_elapsed <- difftime(Sys.time(), t2, units = "mins")
cat(sprintf("\n  LHS complete: %.1f minutes\n", as.numeric(t2_elapsed)))

# --- Summary ---
cat("\n=== LHS Results Summary ===\n")
cat(sprintf("  Total samples:     %d\n", nrow(lhs_results)))
cat(sprintf("  Score range:       %.4f – %.4f\n",
            min(lhs_results$score), max(lhs_results$score)))
cat(sprintf("  Median score:      %.4f\n", median(lhs_results$score)))

# Quick filter preview
n_coexist <- sum(lhs_results$coexistence <= 0.01)
n_zoo     <- sum(lhs_results$zoo_comp <= 0.3)
n_both    <- sum(lhs_results$coexistence <= 0.01 & lhs_results$zoo_comp <= 0.3)
cat(sprintf("  Coexistence pass:  %d / %d\n", n_coexist, nrow(lhs_results)))
cat(sprintf("  Zoo comp pass:     %d / %d\n", n_zoo, nrow(lhs_results)))
cat(sprintf("  Both pass:         %d / %d\n", n_both, nrow(lhs_results)))

if (n_both > 0) {
  top5 <- head(lhs_results[lhs_results$coexistence <= 0.01 &
                            lhs_results$zoo_comp <= 0.3, ], 5)
  top5 <- top5[order(top5$score), ]
  cat("\n  Top 5 candidates:\n")
  for (i in seq_len(nrow(top5))) {
    cat(sprintf("    #%d  score=%.4f  coex=%.4f  zoo=%.4f  stab=%.4f  spec=%.4f\n",
                top5$sample_id[i], top5$score[i], top5$coexistence[i],
                top5$zoo_comp[i], top5$stability[i], top5$spectrum[i]))
  }
}

total_elapsed <- difftime(Sys.time(), t1, units = "hours")
cat(sprintf("\n  Total elapsed: %.1f hours\n", as.numeric(total_elapsed)))

results_file <- file.path(cache_dir, "lhs", "lhs_results.rds")
cat(sprintf("  Results saved:  %s\n", results_file))
cat(sprintf("  Benchmark at:   %s\n",
            file.path(cache_dir, "benchmark",
                      sprintf("benchmark_sst%.0f.rds", sst))))

cat("\nNext step: run refinement on top candidates (400yr runs).\n")
cat("=============================================================\n")
