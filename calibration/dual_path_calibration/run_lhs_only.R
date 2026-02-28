#!/usr/bin/env Rscript
# =============================================================================
# Run LHS Exploration Only (Phases 1-2)
# =============================================================================
#
# All 8 parameter types are group-specific (24 dimensions):
#   PPMR, FeedWidth, K_growth, f_M, repro_eff, Wmat, ZSpre, ZSexp
#   each x 3 fish groups (Small, Medium, Large)
#
# Constraints: per-group R_frac >= 0.15, Wmat_S <= Wmat_M <= Wmat_L
#
# Usage:
#   source("run_lhs_only.R")          # in RStudio
#   Rscript run_lhs_only.R 500 14     # from terminal
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
n_samples  <- if (length(args) >= 1) as.integer(args[1]) else 500L
n_workers  <- if (length(args) >= 2) as.integer(args[2]) else 14L

cat("=============================================================\n")
cat("ZooMSS Fish Reproduction Calibration - LHS Exploration Only\n")
cat(sprintf("  Samples:    %d\n", n_samples))
cat(sprintf("  Workers:    %d\n", n_workers))
cat(sprintf("  Parameters: 24 (all group-specific)\n"))
cat(sprintf("  Constraint: per-group R_frac >= 0.15\n"))
cat(sprintf("  Forcing:    Seasonal (SST amp=4, chl amp=0.5)\n"))
cat(sprintf("  Benchmark:  400yr, assess final 100yr\n"))
cat(sprintf("  LHS:        300yr, assess final 100yr\n"))
cat(sprintf("  Started:    %s\n", Sys.time()))
cat("=============================================================\n\n")

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
  chl_levels = chl_levels, sst = sst,
  n_years = 400, dt = 0.1,
  sst_amplitude = sst_amplitude, chl_amplitude = chl_amplitude,
  assess_years = 100,
  cache_dir = file.path(cache_dir, "benchmark"),
  n_workers = n_workers, force_rerun = FALSE
)

t1_elapsed <- difftime(Sys.time(), t1, units = "mins")
cat(sprintf("  Benchmark complete: %.1f minutes\n\n", as.numeric(t1_elapsed)))

# --- Phase 2: LHS ---
cat("=== Phase 2: LHS Parameter Exploration (300yr, seasonal, 24 dims) ===\n")
t2 <- Sys.time()

lhs_samples <- generate_lhs_samples(n_samples = n_samples, seed = seed)
cat(sprintf("  Generated %d LHS samples (%d dimensions)\n",
            nrow(lhs_samples), ncol(lhs_samples)))

target_log10 <- c(-1.5, -1.0, -0.5, 0.0, 0.4)
chl_indices <- sapply(target_log10, function(t) {
  which.min(abs(log10(chl_levels) - t))
})
chl_indices <- unique(chl_indices)
cat(sprintf("  Evaluating at %d chl levels: %s\n",
            length(chl_indices),
            paste(sprintf("%.2f", log10(chl_levels[chl_indices])), collapse = ", ")))

est_per_sample_sec <- 5 * 120
est_total_min <- (n_samples / n_workers) * est_per_sample_sec / 60
cat(sprintf("  Estimated runtime: %.0f-%.0f hours\n",
            est_total_min * 0.5 / 60, est_total_min * 1.5 / 60))

lhs_results <- run_lhs_exploration(
  lhs_samples = lhs_samples, benchmark = benchmark,
  chl_indices = chl_indices, n_years = 300,
  n_workers = n_workers,
  cache_dir = file.path(cache_dir, "lhs"),
  batch_size = 50
)

t2_elapsed <- difftime(Sys.time(), t2, units = "mins")
cat(sprintf("\n  LHS complete: %.1f minutes\n", as.numeric(t2_elapsed)))

# --- Summary ---
cat("\n=== LHS Results Summary ===\n")
cat(sprintf("  Total samples:     %d\n", nrow(lhs_results)))
cat(sprintf("  Score range:       %.4f - %.4f\n",
            min(lhs_results$score), max(lhs_results$score)))
cat(sprintf("  Median score:      %.4f\n", median(lhs_results$score)))

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
cat(sprintf("  Results saved:  %s\n", file.path(cache_dir, "lhs", "lhs_results.rds")))
cat("\nNext step: run refinement on top candidates (400yr runs).\n")
cat("=============================================================\n")
