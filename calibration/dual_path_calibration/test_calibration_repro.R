#!/usr/bin/env Rscript
# =============================================================================
# Test Run: Fish Reproduction Calibration Pipeline (24 dimensions)
# =============================================================================

cat("=============================================================\n")
cat("ZooMSS Fish Reproduction Calibration - TEST RUN (24-dim)\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat("=============================================================\n\n")

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(zoomss)
}

calib_source <- file.path("R", "zoomss_calibration_repro.R")
if (!file.exists(calib_source)) {
  calib_source <- "zoomss_calibration_repro.R"
  if (!file.exists(calib_source)) stop("Cannot find zoomss_calibration_repro.R")
}
source(calib_source)

required_pkgs <- c("lhs", "future", "furrr")
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))

test_cache_dir <- file.path(tempdir(), "zoomss_calib_test")
if (dir.exists(test_cache_dir)) unlink(test_cache_dir, recursive = TRUE)
dir.create(test_cache_dir, recursive = TRUE)

test_n_workers  <- 1L
test_sst        <- 15
test_n_years_bm <- 30
test_n_years_sc <- 20
test_n_lhs      <- 5L
test_seed       <- 42L
test_dt         <- 0.1
test_assess_yr  <- 10
test_sst_amp    <- 4
test_chl_amp    <- 0.5
test_chl_levels <- 10^c(-1.0, 0.0)

pass_count <- 0; fail_count <- 0; test_log <- character(0)
report <- function(test_name, passed, detail = "") {
  status <- if (passed) "PASS" else "FAIL"
  msg <- sprintf("[%s] %s %s", status, test_name, detail)
  cat(msg, "\n")
  test_log <<- c(test_log, msg)
  if (passed) pass_count <<- pass_count + 1 else fail_count <<- fail_count + 1
}


# =============================================================================
# Test 1: Parameter space (24 dimensions, all group-specific)
# =============================================================================
cat("\n--- Test 1: Parameter Space Definition ---\n")
tryCatch({
  ps <- repro_param_space()
  report("param_space has 24 parameters", nrow(ps) == 24, sprintf("(got %d)", nrow(ps)))
  report("all group-specific (no shared)", all(!is.na(ps$group_idx)))
  report("8 param types x 3 groups",
         length(unique(sub("_[SML]$", "", ps$name))) == 8)
  report("all lower < upper", all(ps$lower < ps$upper))
  # Check PPMR, FeedWidth, f_M, K_growth, repro_eff each have 3 variants
  for (ptype in c("PPMR", "FeedWidth", "f_M", "K_growth", "repro_eff")) {
    n <- sum(grepl(paste0("^", ptype, "_"), ps$name))
    report(sprintf("%s has 3 group variants", ptype), n == 3)
  }
}, error = function(e) report("param_space", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 2: Energy constraint (per-group, R_frac >= 0.15)
# =============================================================================
cat("\n--- Test 2: Per-Group Energy Constraint ---\n")
tryCatch({
  # All groups valid
  par_ok <- c(f_M_S = 0.40, K_growth_S = 0.30,
              f_M_M = 0.35, K_growth_M = 0.35,
              f_M_L = 0.30, K_growth_L = 0.40)
  report("all groups valid passes", check_energy_constraint(par_ok))

  # One group invalid (L: R_frac = 1 - 0.60 - 0.40 = 0.00)
  par_bad <- c(f_M_S = 0.40, K_growth_S = 0.30,
               f_M_M = 0.35, K_growth_M = 0.35,
               f_M_L = 0.60, K_growth_L = 0.40)
  report("one group invalid fails", !check_energy_constraint(par_bad))

  # Edge case: R_frac = 0.15 exactly (FP: 1 - 0.45 - 0.40 = 0.14999...)
  par_edge <- c(f_M_S = 0.40, K_growth_S = 0.30,
                f_M_M = 0.35, K_growth_M = 0.35,
                f_M_L = 0.45, K_growth_L = 0.40)
  report("edge case R_frac=0.15 passes with tolerance",
         check_energy_constraint(par_edge))

  # Below threshold: R_frac = 0.14
  par_low <- c(f_M_S = 0.40, K_growth_S = 0.30,
               f_M_M = 0.35, K_growth_M = 0.35,
               f_M_L = 0.46, K_growth_L = 0.40)
  report("R_frac=0.14 fails", !check_energy_constraint(par_low))
}, error = function(e) report("energy constraint", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 3: Wmat ordering
# =============================================================================
cat("\n--- Test 3: Wmat Ordering Constraint ---\n")
tryCatch({
  report("ordered passes", check_wmat_ordering(c(Wmat_S=-1, Wmat_M=1, Wmat_L=3)))
  report("equal passes",   check_wmat_ordering(c(Wmat_S=1, Wmat_M=1, Wmat_L=1)))
  report("reversed fails", !check_wmat_ordering(c(Wmat_S=3, Wmat_M=1, Wmat_L=-1)))
  report("S>M fails",      !check_wmat_ordering(c(Wmat_S=2, Wmat_M=1, Wmat_L=3)))
}, error = function(e) report("Wmat ordering", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 4: Rounding
# =============================================================================
cat("\n--- Test 4: Parameter Rounding ---\n")
tryCatch({
  par_raw <- c(PPMR_S = 145.7, PPMR_M = 200.3, PPMR_L = 312.8,
               FeedWidth_S = 1.756, K_growth_M = 0.2173, f_M_L = 0.3114)
  par_r <- round_params(par_raw)
  report("PPMR_S -> integer", par_r["PPMR_S"] == 146)
  report("PPMR_L -> integer", par_r["PPMR_L"] == 313)
  report("FeedWidth_S -> 2dp", par_r["FeedWidth_S"] == 1.76)
  report("K_growth_M -> 2dp", par_r["K_growth_M"] == 0.22)
}, error = function(e) report("rounding", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 5: Seasonal environment helper
# =============================================================================
cat("\n--- Test 5: Seasonal Environment ---\n")
tryCatch({
  ip <- create_seasonal_input(n_years = 10, dt = 0.1,
                              base_sst = 15, base_chl = 1.0,
                              sst_amplitude = 4, chl_amplitude = 0.5)
  report("returns list with time/sst/chl", is.list(ip) && all(c("time","sst","chl") %in% names(ip)))
  report("same lengths", length(ip$time) == length(ip$sst) && length(ip$time) == length(ip$chl))
  report("sst varies", sd(ip$sst) > 0)
  report("chl always positive", all(ip$chl > 0))
}, error = function(e) report("seasonal env", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 6: LHS generation (24 dims)
# =============================================================================
cat("\n--- Test 6: LHS Sample Generation ---\n")
tryCatch({
  lhs <- generate_lhs_samples(n_samples = test_n_lhs, seed = test_seed)
  ps <- repro_param_space()
  report("correct dims (n x 24)", nrow(lhs) == test_n_lhs && ncol(lhs) == 24,
         sprintf("(%d x %d)", nrow(lhs), ncol(lhs)))
  report("column names match", all(names(lhs) == ps$name))

  all_energy <- all(sapply(seq_len(nrow(lhs)), function(i) check_energy_constraint(lhs[i,])))
  report("all satisfy energy constraint", all_energy)

  all_wmat <- all(sapply(seq_len(nrow(lhs)), function(i) check_wmat_ordering(lhs[i,])))
  report("all satisfy Wmat ordering", all_wmat)

  in_bounds <- TRUE
  for (i in seq_len(ncol(lhs))) {
    if (any(lhs[,i] < ps$lower[i] - 1e-10) || any(lhs[,i] > ps$upper[i] + 1e-10)) {
      in_bounds <- FALSE; break
    }
  }
  report("all within bounds", in_bounds)
}, error = function(e) report("LHS generation", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 7: apply_repro_params (group-specific)
# =============================================================================
cat("\n--- Test 7: apply_repro_params ---\n")
tryCatch({
  Groups <- getGroups()
  fish_idx <- which(Groups$Type == "Fish")
  par <- as.numeric(lhs[1,]); names(par) <- names(lhs)
  Groups_mod <- apply_repro_params(par, Groups)
  par_r <- round_params(par)

  report("PPMR differs across fish groups",
         !all(Groups_mod$PPMR[fish_idx] == Groups_mod$PPMR[fish_idx[1]]) ||
           par_r["PPMR_S"] == par_r["PPMR_M"],  # may be equal by chance
         "(may match by chance)")
  report("FeedWidth group-specific",
         Groups_mod$FeedWidth[fish_idx[1]] == par_r["FeedWidth_S"] &&
           Groups_mod$FeedWidth[fish_idx[2]] == par_r["FeedWidth_M"])
  report("f_M group-specific",
         Groups_mod$f_M[fish_idx[1]] == par_r["f_M_S"] &&
           Groups_mod$f_M[fish_idx[3]] == par_r["f_M_L"])
  report("K_growth group-specific",
         Groups_mod$K_growth[fish_idx[2]] == par_r["K_growth_M"])
  report("repro_eff group-specific",
         Groups_mod$repro_eff[fish_idx[1]] == par_r["repro_eff_S"])
  report("repro_on set to 1", all(Groups_mod$repro_on[fish_idx] == 1L))

  zoo_rows <- which(Groups$Type == "Zooplankton")
  report("zooplankton unchanged",
         all(sapply(c("Species","Type","W0","Wmax","GrossGEscale"), function(col)
           identical(Groups_mod[[col]][zoo_rows], Groups[[col]][zoo_rows]))))
}, error = function(e) report("apply_repro_params", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 8: Single model run
# =============================================================================
cat("\n--- Test 8: Single Model Run ---\n")
tryCatch({
  ip <- create_seasonal_input(n_years = test_n_years_bm, dt = test_dt,
                              base_sst = test_sst, base_chl = test_chl_levels[1],
                              sst_amplitude = test_sst_amp, chl_amplitude = test_chl_amp)
  mdl <- zoomss_model(input_params = ip, Groups = Groups_mod, isave = 2)
  report("biomass is 3D", length(dim(mdl$biomass)) == 3)
  avg <- averageTimeSeries(mdl, var = "biomass", n_years = test_assess_yr)
  fish_bm <- rowSums(avg)[mdl$param$fish_grps]
  report("fish biomass > 0", all(fish_bm > 0),
         sprintf("(%s)", paste(sprintf("%.2e", fish_bm), collapse=", ")))
}, error = function(e) report("model run", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 9: Benchmark generation
# =============================================================================
cat("\n--- Test 9: Benchmark Generation ---\n")
benchmark <- NULL
tryCatch({
  benchmark <- generate_legacy_benchmark(
    chl_levels = test_chl_levels, sst = test_sst,
    n_years = test_n_years_bm, dt = test_dt,
    sst_amplitude = test_sst_amp, chl_amplitude = test_chl_amp,
    assess_years = test_assess_yr,
    cache_dir = file.path(test_cache_dir, "benchmark"),
    n_workers = test_n_workers, force_rerun = TRUE)
  report("benchmark is a list", is.list(benchmark))
  report("zoo proportions sum to ~1",
         all(abs(rowSums(benchmark$zoo_proportions, na.rm=TRUE) - 1) < 0.01))
  report("stores seasonal params",
         !is.null(benchmark$sst_amplitude) && !is.null(benchmark$chl_amplitude))
}, error = function(e) report("benchmark", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 10: Objective function
# =============================================================================
cat("\n--- Test 10: Objective Function ---\n")
tryCatch({
  if (is.null(benchmark)) stop("Benchmark not available")
  par <- as.numeric(lhs[1,]); names(par) <- names(lhs)

  result <- repro_objective(par = par, benchmark = benchmark,
    chl_indices = 1:2, n_years = test_n_years_sc, dt = test_dt,
    assess_years = test_assess_yr, return_details = TRUE)
  report("returns list with score", is.list(result) && is.numeric(result$score))
  report("5 metrics in [0,1]",
         length(result$metric_scores) == 5 &&
           all(result$metric_scores >= 0 & result$metric_scores <= 1))

  # Energy violation (one group)
  par_bad <- par; par_bad["f_M_L"] <- 0.70; par_bad["K_growth_L"] <- 0.45
  report("energy violation -> 1e6",
         repro_objective(par=par_bad, benchmark=benchmark,
                         chl_indices=1, n_years=test_n_years_sc,
                         assess_years=test_assess_yr) == 1e6)

  # Wmat violation
  par_wmat <- par; par_wmat["Wmat_S"] <- 3.0; par_wmat["Wmat_L"] <- -1.0
  report("Wmat violation -> 1e6",
         repro_objective(par=par_wmat, benchmark=benchmark,
                         chl_indices=1, n_years=test_n_years_sc,
                         assess_years=test_assess_yr) == 1e6)
}, error = function(e) report("objective", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 11: LHS exploration
# =============================================================================
cat("\n--- Test 11: LHS Exploration ---\n")
lhs_results <- NULL
tryCatch({
  if (is.null(benchmark)) stop("Benchmark not available")
  lhs_results <- run_lhs_exploration(
    lhs_samples = lhs, benchmark = benchmark,
    chl_indices = 1:2, n_years = test_n_years_sc,
    n_workers = test_n_workers,
    cache_dir = file.path(test_cache_dir, "lhs"), batch_size = 3)
  report("correct rows", nrow(lhs_results) == test_n_lhs)
  report("all scores finite", all(is.finite(lhs_results$score)))
}, error = function(e) report("LHS exploration", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 12: Filtering
# =============================================================================
cat("\n--- Test 12: Filtering ---\n")
candidates <- NULL
tryCatch({
  if (is.null(lhs_results)) stop("LHS results not available")
  candidates <- filter_lhs_candidates(lhs_results, max_coexistence=1, max_zoo_comp=1, top_n=2)
  report("filter works", is.data.frame(candidates) && nrow(candidates) <= 2)
  report("sorted by score", nrow(candidates) <= 1 || all(diff(candidates$score) >= 0))
}, error = function(e) report("filtering", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 13: Refinement (interface)
# =============================================================================
cat("\n--- Test 13: Refinement (interface) ---\n")
tryCatch({
  report("refine_candidate exists", is.function(refine_candidate))
  if (!is.null(candidates) && nrow(candidates) > 0) {
    param_names <- names(lhs)
    par_init <- as.numeric(candidates[1, param_names]); names(par_init) <- param_names
    details <- repro_objective(par=par_init, benchmark=benchmark,
                               chl_indices=1, n_years=2, dt=test_dt,
                               assess_years=1, return_details=TRUE)
    report("objective on candidate works", is.numeric(details$score))
  }
}, error = function(e) report("refinement", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Test 14: Yield curves
# =============================================================================
cat("\n--- Test 14: Yield Curves ---\n")
tryCatch({
  par <- as.numeric(lhs[1,]); names(par) <- names(lhs)
  yd <- yield_curve_validation(par=par, fmort_levels=c(0,0.5,1),
    chl=1, sst=test_sst, n_years=test_n_years_sc, dt=test_dt,
    sst_amplitude=test_sst_amp, chl_amplitude=test_chl_amp,
    assess_years=test_assess_yr)
  report("returns data.frame", is.data.frame(yd))
  report("9 rows (3F x 3fish)", nrow(yd) == 9, sprintf("(got %d)", nrow(yd)))
  report("yield at F=0 is 0", all(yd$Yield[yd$Fmort==0] == 0, na.rm=TRUE))
}, error = function(e) report("yield curves", FALSE, paste("ERROR:", e$message)))


# =============================================================================
# Summary
# =============================================================================
cat("\n=============================================================\n")
cat(sprintf("TEST SUMMARY: %d PASSED, %d FAILED (of %d)\n",
            pass_count, fail_count, pass_count + fail_count))
cat("=============================================================\n")
if (fail_count > 0) {
  cat("\nFailed tests:\n")
  for (f in grep("^\\[FAIL\\]", test_log, value=TRUE)) cat("  ", f, "\n")
}
writeLines(test_log, file.path(test_cache_dir, "test_log.txt"))
cat(sprintf("Completed: %s\n", Sys.time()))
unlink(test_cache_dir, recursive = TRUE)

if (fail_count > 0) {
  cat("\n*** TESTS FAILED ***\n")
  if (!interactive()) quit(save="no", status=1)
} else {
  cat("\n*** ALL TESTS PASSED ***\n")
  if (!interactive()) quit(save="no", status=0)
}
