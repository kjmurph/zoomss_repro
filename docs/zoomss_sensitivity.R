# =============================================================================
# ZooMSS Dual-Pathway Sensitivity Analysis
# =============================================================================
#
# Workflow:
#   setup_parallel()                                 # detect cores, pick method
#   test_run()                                       # serial smoke test + timing
#   test_parallel()                                  # parallel smoke test
#   run_morris_stage(outfile = "morris.rds")         # overnight on laptop
#   morris_res <- readRDS("morris.rds")
#   plot_morris(morris_res)
#   get_top_params_morris(morris_res)
#   run_lhs_stage(n_lhs = 500, outfile = "lhs.rds") # on VM
#   resume_lhs_stage("lhs.rds")                     # if crashed
#   run_sobol_stage(lhs_file = "lhs.rds")
#   run_optimise_stage(lhs_file = "lhs.rds")
#
# =============================================================================

# devtools::load_all("path/to/zoomss_repro")   # ← load your package first
library(lhs)
library(sensitivity)
library(dplyr)
library(ggplot2)

# =============================================================================
# RUN CONSTANTS  — adjust to match your setup
# =============================================================================

SA_DT         <- 0.1    # years; time step for sensitivity runs
SA_ISAVE      <- 10     # save every N time steps (50 yr / 0.1 dt = 500 steps → 50 saved)
SA_T_RUN      <- 50     # years per run
SA_SST        <- 10     # °C; held constant during sensitivity runs
SA_T_STABLE_F <- 0.4    # fraction of saved steps used for stability assessment

# Chlorophyll gradient for zooplankton composition checks (mg m^-3)
SA_CHL_GRADIENT <- 10^seq(-1.5, 1.5, length.out = 7)

# A group is considered extinct if mean biomass in the stable window falls below:
SA_BIOMASS_THRESHOLD <- 1e-8   # g m^-2 (integrated over size; adjust to your units)

# =============================================================================
# PARALLEL BACKEND
# =============================================================================
# Call setup_parallel() once at session start.
# Mac/Linux: fork-based mclapply — no cluster overhead, objects shared for free
# Windows:   PSOCK cluster — must export objects explicitly (handled automatically)
# VM/HPC:    set n_cores manually; fork is safe on all Linux VMs

setup_parallel <- function(n_cores = NULL, verbose = TRUE) {
  os        <- .Platform$OS.type   # "unix" or "windows"
  available <- parallel::detectCores(logical = FALSE)  # physical cores only

  if (is.null(n_cores)) n_cores <- max(1L, available - 1L)
  n_cores <- as.integer(min(n_cores, available))
  method  <- if (os == "windows") "psock" else "fork"

  if (verbose) {
    message(sprintf("[parallel] OS: %s | method: %s | cores: %d / %d physical",
                    os, method, n_cores, available))
    if (method == "psock")
      message("[parallel] PSOCK cluster: objects/packages exported automatically via par_lapply().")
  }

  cfg <- list(n_cores = n_cores, method = method)
  assign("PAR_CONFIG", cfg, envir = .GlobalEnv)
  invisible(cfg)
}

# Internal parallel map.  Handles fork vs PSOCK automatically.
# extra_export: names of additional .GlobalEnv objects PSOCK workers need.
par_lapply <- function(X, FUN,
                       extra_export = character(0),
                       config       = get0("PAR_CONFIG", envir = .GlobalEnv)) {

  if (is.null(config) || config$n_cores == 1L)
    return(lapply(X, FUN))

  if (config$method == "fork")
    return(parallel::mclapply(X, FUN,
                              mc.cores      = config$n_cores,
                              mc.preschedule = TRUE))

  # PSOCK (Windows)
  cl <- parallel::makeCluster(config$n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  always_export <- c(
    # Parameter space constants
    "PARAM_NAMES", "N_PARAMS", "LOWER", "UPPER", "FISH_GROUPS",
    # Run constants
    "SA_DT", "SA_ISAVE", "SA_T_RUN", "SA_SST",
    "SA_T_STABLE_F", "SA_CHL_GRADIENT", "SA_BIOMASS_THRESHOLD",
    # Helper functions used inside run_one
    "unit_to_params", "patch_group_inputs", "run_one", "run_model_once"
  )
  parallel::clusterExport(cl, unique(c(always_export, extra_export)),
                          envir = .GlobalEnv)

  parallel::clusterEvalQ(cl, {
    # ---- load your package on each worker ----
    # devtools::load_all("path/to/zoomss_repro")
    library(dplyr)
  })

  parallel::parLapply(cl, X, FUN)
}

# =============================================================================
# PARAMETER SPACE  (24 parameters: 8 per fish group × 3 groups)
# =============================================================================

FISH_GROUPS <- c("Fish_Small", "Fish_Med", "Fish_Large")

param_bounds <- list(
  Fish_Small = list(
    Wmat          = c(-2.0,  1.5),   # log10(g); W0=-3, Wmax=2
    log_PPMR      = c( 1.0,  3.0),   # log10(PPMR); 10–1000
    FeedWidth     = c( 0.5,  2.5),   # feeding kernel width (decades)
    f_M           = c( 0.25, 0.65),  # metabolic fraction
    K_frac        = c( 0.2,  0.8),   # K_growth = K_frac * (1 - f_M)
    log_repro_eff = c(-2.0,  1.0),   # log10(repro_eff)
    log_ZSpre     = c(-2.5,  0.0),   # log10(ZSpre)
    ZSexp         = c( 0.1,  0.8)
  ),
  Fish_Med = list(
    Wmat          = c( 0.5,  3.5),   # W0=-3, Wmax=4
    log_PPMR      = c( 1.0,  3.0),
    FeedWidth     = c( 0.5,  2.5),
    f_M           = c( 0.25, 0.65),
    K_frac        = c( 0.2,  0.8),
    log_repro_eff = c(-2.0,  1.0),
    log_ZSpre     = c(-2.5,  0.0),
    ZSexp         = c( 0.1,  0.8)
  ),
  Fish_Large = list(
    Wmat          = c( 2.5,  5.5),   # W0=-3, Wmax=6
    log_PPMR      = c( 1.0,  3.0),
    FeedWidth     = c( 0.5,  2.5),
    f_M           = c( 0.25, 0.65),
    K_frac        = c( 0.2,  0.8),
    log_repro_eff = c(-2.0,  1.0),
    log_ZSpre     = c(-2.5,  0.0),
    ZSexp         = c( 0.1,  0.8)
  )
)

PARAM_NAMES <- unlist(lapply(FISH_GROUPS, function(g)
  paste0(g, "_", names(param_bounds[[g]]))))
N_PARAMS    <- length(PARAM_NAMES)   # 24

LOWER <- unlist(lapply(FISH_GROUPS, function(g) sapply(param_bounds[[g]], `[`, 1)))
UPPER <- unlist(lapply(FISH_GROUPS, function(g) sapply(param_bounds[[g]], `[`, 2)))
names(LOWER) <- PARAM_NAMES
names(UPPER) <- PARAM_NAMES

# =============================================================================
# PARAMETER TRANSLATION
# =============================================================================
# Maps [0,1]^24 hypercube sample → named list of actual parameter values.
# K_growth reparameterised via K_frac so R_frac = 1 - f_M - K_growth > 0 always.

unit_to_params <- function(x_unit) {
  x <- LOWER + x_unit * (UPPER - LOWER)
  names(x) <- PARAM_NAMES
  lapply(FISH_GROUPS, function(g) {
    f_M    <- x[[paste0(g, "_f_M")]]
    K_frac <- x[[paste0(g, "_K_frac")]]
    list(
      Wmat      = x[[paste0(g, "_Wmat")]],
      PPMR      = 10^x[[paste0(g, "_log_PPMR")]],
      FeedWidth = x[[paste0(g, "_FeedWidth")]],
      f_M       = f_M,
      K_growth  = K_frac * (1 - f_M),   # R_frac = 1 - f_M - K_growth > 0
      repro_eff = 10^x[[paste0(g, "_log_repro_eff")]],
      ZSpre     = 10^x[[paste0(g, "_log_ZSpre")]],
      ZSexp     = x[[paste0(g, "_ZSexp")]]
    )
  }) |> setNames(FISH_GROUPS)
}

# Patch a copy of base_groups with the translated parameters
patch_group_inputs <- function(base_groups, param_list) {
  gi <- base_groups
  for (g in FISH_GROUPS) {
    idx <- which(gi$Species == g)
    p   <- param_list[[g]]
    gi$Wmat[idx]      <- p$Wmat
    gi$PPMR[idx]      <- p$PPMR
    gi$FeedWidth[idx] <- p$FeedWidth
    gi$f_M[idx]       <- p$f_M
    gi$K_growth[idx]  <- p$K_growth
    gi$repro_eff[idx] <- p$repro_eff
    gi$ZSpre[idx]     <- p$ZSpre
    gi$ZSexp[idx]     <- p$ZSexp
  }
  gi
}

# =============================================================================
# MODEL WRAPPER
# =============================================================================
# Single call to the actual ZooMSS API.
# Returns [saved_steps × groups] biomass matrix, or NULL on error.
#
# suppressMessages + capture.output silence the progress bar and cat() calls
# from zoomss_params / zoomss_run, which would flood the console and cause
# issues inside parallel workers.

run_model_once <- function(groups, chl, sst = SA_SST,
                           t_run = SA_T_RUN, dt = SA_DT, isave = SA_ISAVE) {
  time_vec    <- seq(0, t_run, by = dt)
  input_params <- createInputParams(time_vec,
                                    sst = rep(sst, length(time_vec)),
                                    chl = rep(chl, length(time_vec)))
  result <- suppressMessages(
    capture.output(
      mdl <- zoomss_model(input_params, Groups = groups, isave = isave)
    )
  )
  # reduceSize returns [saved_steps × groups] by summing across size bins
  reduceSize(mdl$biomass)
}

# =============================================================================
# SINGLE PARAMETER SET: RUN + COMPUTE METRICS
# =============================================================================

run_one <- function(x_unit, base_groups,
                    chl       = median(SA_CHL_GRADIENT),
                    sst       = SA_SST,
                    t_run     = SA_T_RUN,
                    chl_sweep = FALSE) {
  tryCatch({
    gi     <- patch_group_inputs(base_groups, unit_to_params(x_unit))
    bio_ts <- run_model_once(gi, chl = chl, sst = sst, t_run = t_run)
    # bio_ts: [n_saved × n_groups], rows = time, cols = groups in gi$Species order

    fish_idx <- which(gi$Type == "Fish")
    zoo_idx  <- which(gi$Type == "Zooplankton")
    n_saved  <- nrow(bio_ts)
    t0       <- max(1L, floor(n_saved * (1 - SA_T_STABLE_F)))
    bio_fin  <- bio_ts[t0:n_saved, , drop = FALSE]

    # (a) Coexistence: all fish groups above threshold throughout stable window
    fish_survive <- all(apply(bio_fin[, fish_idx, drop = FALSE], 2,
                              function(b) all(b > SA_BIOMASS_THRESHOLD)))

    # (b) Stability: CV in (0.005, 2.0) — oscillating but not exploding
    fish_cv <- apply(bio_fin[, fish_idx, drop = FALSE], 2,
                     function(b) sd(b) / (mean(b) + 1e-30))
    fish_stable <- fish_survive &&
      all(fish_cv > 0.005) && all(fish_cv < 2.0)

    # (c) Shannon entropy of mean fish biomass — evenness proxy
    fish_mean    <- colMeans(bio_fin[, fish_idx, drop = FALSE])
    fish_prop    <- fish_mean / (sum(fish_mean) + 1e-30)
    fish_entropy <- -sum(fish_prop * log(fish_prop + 1e-30))

    # (d) Zooplankton chl-gradient MAE vs reference (optional — ~7× slower)
    zoo_mae_chl <- NA_real_
    if (chl_sweep && fish_survive) {
      zoo_ref <- get0("ZOO_REF_COMP", envir = .GlobalEnv)
      if (!is.null(zoo_ref)) {
        zoo_comp <- do.call(rbind, lapply(SA_CHL_GRADIENT, function(cv) {
          bt <- run_model_once(gi, chl = cv, sst = sst, t_run = t_run)
          bio_end <- bt[nrow(bt), zoo_idx]
          bio_end / (sum(bio_end) + 1e-30)
        }))
        zoo_mae_chl <- mean(abs(zoo_comp - zoo_ref))
      }
    }

    list(
      fish_survive  = as.integer(fish_survive),
      fish_stable   = as.integer(fish_stable),
      fish_entropy  = fish_entropy,
      fish_cv_small = fish_cv[1],
      fish_cv_med   = fish_cv[2],
      fish_cv_large = fish_cv[3],
      zoo_mae_chl   = zoo_mae_chl
    )

  }, error = function(e) {
    message("[run_one] ", conditionMessage(e))
    NULL
  })
}

# Collapse a list of run_one() outputs to a data frame (NAs for failed runs)
results_to_df <- function(res_list) {
  null_row <- list(fish_survive  = NA_integer_, fish_stable   = NA_integer_,
                   fish_entropy  = NA_real_,    fish_cv_small = NA_real_,
                   fish_cv_med   = NA_real_,    fish_cv_large = NA_real_,
                   zoo_mae_chl   = NA_real_)
  rows <- lapply(res_list, function(r) if (is.null(r)) null_row else r)
  as.data.frame(do.call(rbind, lapply(rows, unlist)))
}

# =============================================================================
# TEST RUN  (serial — accurate timing, checks output structure)
# =============================================================================
# Always run this before any batch.  Uses a short t_run so it completes quickly,
# then projects wall-time to the full SA_T_RUN.

test_run <- function(base_groups = NULL,
                     n_test      = 5,
                     t_run       = 5,     # short for speed
                     verbose     = TRUE) {

  if (is.null(base_groups)) { data(GroupInputs); base_groups <- GroupInputs }

  message("\n", strrep("=", 62))
  message("ZooMSS SENSITIVITY — SERIAL TEST  (n=", n_test, ", t=", t_run, " yr)")
  message(strrep("=", 62))

  cfg <- get0("PAR_CONFIG", envir = .GlobalEnv)
  if (is.null(cfg)) message("[test] No parallel config — call setup_parallel() first.")
  else message("[test] Parallel config: ", cfg$method, " | cores: ", cfg$n_cores)

  # Check package is loaded
  if (!exists("zoomss_model")) {
    stop("[test] zoomss_model() not found. Load your package first:\n",
         "  devtools::load_all('path/to/zoomss_repro')")
  }
  message("[test] zoomss_model() found OK")

  # Check default run with unmodified GroupInputs
  message("[test] Running default-parameter smoke test...")
  t_default <- proc.time()["elapsed"]
  tryCatch({
    time_vec     <- seq(0, t_run, by = SA_DT)
    ip           <- createInputParams(time_vec,
                                      rep(SA_SST, length(time_vec)),
                                      rep(median(SA_CHL_GRADIENT), length(time_vec)))
    default_mdl  <- suppressMessages(capture.output(
      zoomss_model(ip, Groups = base_groups, isave = SA_ISAVE)
    ))
    bt           <- reduceSize(default_mdl$biomass)
    message(sprintf("[test] Default run OK: biomass matrix %dx%d, t=%.1fs",
                    nrow(bt), ncol(bt),
                    proc.time()["elapsed"] - t_default))
  }, error = function(e) {
    stop("[test] Default run FAILED: ", conditionMessage(e),
         "\nCheck SA_DT, SA_ISAVE, SA_CHL_GRADIENT constants match your model.")
  })

  # Randomised parameter sets (run serially for clean timing)
  set.seed(42)
  X <- lhs::randomLHS(n_test, N_PARAMS)
  colnames(X) <- PARAM_NAMES

  times   <- numeric(n_test)
  results <- vector("list", n_test)

  for (i in seq_len(n_test)) {
    t0           <- proc.time()["elapsed"]
    results[[i]] <- run_one(X[i, ], base_groups,
                            chl   = median(SA_CHL_GRADIENT),
                            t_run = t_run)
    times[i] <- proc.time()["elapsed"] - t0
    status <- if (is.null(results[[i]])) {
      "ERROR — check run_one() / model API"
    } else {
      sprintf("survive=%d  stable=%d  entropy=%.3f  cv_med=%.3f",
              results[[i]]$fish_survive, results[[i]]$fish_stable,
              results[[i]]$fish_entropy, results[[i]]$fish_cv_med)
    }
    if (verbose) message(sprintf("  run %d/%d  %5.1fs  %s", i, n_test, times[i], status))
  }

  ok     <- sum(!sapply(results, is.null))
  mean_t <- mean(times)
  # Project from short t_run to full SA_T_RUN (time scales ~linearly with n_steps)
  est_full <- mean_t * (SA_T_RUN / t_run)

  message(sprintf("\n  Completed: %d/%d | Mean per run (t=%dyr): %.1fs",
                  ok, n_test, t_run, mean_t))
  message(sprintf("  Projected per run at t=%dyr:  ~%.0fs  (~%.1f min)",
                  SA_T_RUN, est_full, est_full / 60))

  n_cores <- if (is.null(cfg)) 1L else cfg$n_cores
  message("\n  Wall-time estimates (projected, parallel):")
  for (n_runs in c(250, 500, 1000, 6000)) {
    wall <- est_full * n_runs / n_cores / 3600
    message(sprintf("    %5d runs × %2d cores → ~%.1f hrs", n_runs, n_cores, wall))
  }

  if (ok < n_test)
    message("\n  WARNING: ", n_test - ok, " run(s) returned NULL. ",
            "Inspect errors above before proceeding.")
  else
    message("\n  All runs completed successfully. Ready for parallel test.")

  invisible(list(results = results, times = times, est_full_s = est_full))
}

# =============================================================================
# PARALLEL TEST
# =============================================================================
# Must be run after test_run() passes.  Runs the same n_test parameter sets
# in parallel so you can compare:
#   - results match the serial run  (catches export / serialisation errors)
#   - no crashes or silent failures on workers
#   - actual parallel speedup is reasonable
#
# Common failure modes this catches:
#   PSOCK: object not found on worker  → fix always_export in par_lapply()
#   PSOCK: package not loaded          → fix clusterEvalQ block in par_lapply()
#   fork:  worker segfault from        → reduce n_cores, check for C-level
#          large object copying          globals in your package
#   both:  progress bar output floods  → confirms capture.output() is working
#          workers / corrupts results
#   both:  RNG collision (all workers  → set per-worker seeds (handled here)
#          produce identical samples)

test_parallel <- function(base_groups = NULL,
                          n_test      = 8,
                          t_run       = 5,
                          verbose     = TRUE) {

  if (is.null(base_groups)) { data(GroupInputs); base_groups <- GroupInputs }

  cfg <- get0("PAR_CONFIG", envir = .GlobalEnv)
  if (is.null(cfg)) stop("Call setup_parallel() before test_parallel().")
  if (cfg$n_cores < 2L) {
    message("[par_test] Only 1 core configured — skipping parallel test.")
    return(invisible(NULL))
  }

  message("\n", strrep("=", 62))
  message("ZooMSS SENSITIVITY — PARALLEL TEST  (n=", n_test,
          ", t=", t_run, " yr, method=", cfg$method, ", cores=", cfg$n_cores, ")")
  message(strrep("=", 62))

  set.seed(42)
  X <- lhs::randomLHS(n_test, N_PARAMS)
  colnames(X) <- PARAM_NAMES

  # --- Serial reference ---
  message("  Running serial reference...")
  t_serial <- proc.time()["elapsed"]
  serial_results <- lapply(seq_len(n_test), function(i) {
    x <- X[i, ]; names(x) <- PARAM_NAMES
    run_one(x, base_groups, chl = median(SA_CHL_GRADIENT), t_run = t_run)
  })
  t_serial <- proc.time()["elapsed"] - t_serial
  message(sprintf("  Serial: %.1fs total", t_serial))

  # --- Parallel run ---
  message("  Running in parallel...")
  t_par <- proc.time()["elapsed"]
  par_results <- par_lapply(
    seq_len(n_test),
    function(i) {
      x <- X[i, ]; names(x) <- PARAM_NAMES
      run_one(x, base_groups, chl = median(SA_CHL_GRADIENT), t_run = t_run)
    },
    extra_export = "base_groups"
  )
  t_par <- proc.time()["elapsed"] - t_par
  message(sprintf("  Parallel: %.1fs total", t_par))

  # --- Compare results ---
  message("\n  Comparing serial vs parallel results:")
  n_errors    <- 0L
  n_null_ser  <- sum(sapply(serial_results, is.null))
  n_null_par  <- sum(sapply(par_results,    is.null))

  if (n_null_ser > 0) message("  WARNING: ", n_null_ser, " serial run(s) returned NULL")
  if (n_null_par > 0) message("  WARNING: ", n_null_par, " parallel run(s) returned NULL")

  for (i in seq_len(n_test)) {
    s <- serial_results[[i]]
    p <- par_results[[i]]

    if (is.null(s) && is.null(p)) next
    if (is.null(s) != is.null(p)) {
      message(sprintf("  run %d: NULL mismatch (serial=%s, parallel=%s)",
                      i, is.null(s), is.null(p)))
      n_errors <- n_errors + 1L
      next
    }

    # Check key metrics agree within floating-point tolerance
    metrics <- c("fish_survive", "fish_stable", "fish_entropy")
    for (m in metrics) {
      diff_val <- abs(s[[m]] - p[[m]])
      if (diff_val > 1e-6) {
        message(sprintf("  run %d: '%s' mismatch — serial=%.6f, parallel=%.6f (diff=%.2e)",
                        i, m, s[[m]], p[[m]], diff_val))
        n_errors <- n_errors + 1L
      }
    }

    if (verbose) {
      ok_str <- if (n_errors == 0L) "OK" else "MISMATCH"
      message(sprintf("  run %d: survive=%d  stable=%d  entropy=%.4f  [%s]",
                      i, p$fish_survive, p$fish_stable, p$fish_entropy, ok_str))
    }
  }

  # --- Summary ---
  speedup <- t_serial / max(t_par, 0.001)
  ideal   <- min(cfg$n_cores, n_test)

  message(sprintf("\n  Speedup: %.1fx  (ideal with %d cores: %.1fx)",
                  speedup, cfg$n_cores, min(ideal, n_test / 1.0)))
  message(sprintf("  Efficiency: %.0f%%", 100 * speedup / ideal))

  if (speedup < 1.2 && cfg$n_cores > 2L) {
    message("  NOTE: Low speedup — likely dominated by object serialisation or ",
            "worker startup. Consider reducing base_groups size, or check that ",
            "large objects (kernels) are not being re-serialised each call.")
  }

  if (n_errors == 0L) {
    message("\n  PASSED: parallel results match serial. Ready for full batch runs.")
  } else {
    message("\n  FAILED: ", n_errors, " result mismatch(es). ",
            "Fix parallel export/package issues before running full batch.")
    message("  Common fixes:")
    message("    PSOCK: add missing object names to extra_export in par_lapply()")
    message("    PSOCK: add missing library() call to clusterEvalQ block in par_lapply()")
    message("    Both:  confirm capture.output() suppresses progress bar output")
  }

  invisible(list(
    serial   = serial_results,
    parallel = par_results,
    n_errors = n_errors,
    speedup  = speedup
  ))
}

# =============================================================================
# ZOOPLANKTON REFERENCE
# =============================================================================

compute_zoo_reference <- function(base_groups = NULL,
                                  chl_values  = SA_CHL_GRADIENT,
                                  sst         = SA_SST,
                                  t_run       = SA_T_RUN) {
  if (is.null(base_groups)) { data(GroupInputs); base_groups <- GroupInputs }
  zoo_idx <- which(base_groups$Type == "Zooplankton")
  message("[zoo_ref] Running ", length(chl_values), " reference chl points...")
  ref <- do.call(rbind, lapply(chl_values, function(cv) {
    bt  <- run_model_once(base_groups, chl = cv, sst = sst, t_run = t_run)
    bio <- bt[nrow(bt), zoo_idx]
    bio / (sum(bio) + 1e-30)
  }))
  assign("ZOO_REF_COMP", ref, envir = .GlobalEnv)
  message("[zoo_ref] Done. ZOO_REF_COMP stored in global environment.")
  invisible(ref)
}

# =============================================================================
# STAGE 1: MORRIS SCREENING
# =============================================================================
# Total runs = n_morris × (N_PARAMS + 1).  With n_morris=10, N_PARAMS=24: 250 runs.
# The binary metric (stable coexistence yes/no) is fast to compute.
# Chl is fixed at the median for speed; gradient sweep not needed here.

run_morris_stage <- function(base_groups = NULL,
                             n_morris    = 10,
                             outfile     = "morris_results.rds",
                             sst         = SA_SST,
                             t_run       = SA_T_RUN,
                             chl         = median(SA_CHL_GRADIENT)) {

  if (is.null(base_groups)) { data(GroupInputs); base_groups <- GroupInputs }
  cfg     <- get0("PAR_CONFIG", envir = .GlobalEnv)
  n_runs  <- n_morris * (N_PARAMS + 1L)
  n_cores <- if (is.null(cfg)) 1L else cfg$n_cores

  message("\n", strrep("=", 62))
  message("STAGE 1: Morris screening")
  message(sprintf("  n_morris=%d | N_PARAMS=%d | total_runs=%d | cores=%d",
                  n_morris, N_PARAMS, n_runs, n_cores))
  if (!is.null(outfile)) message("  Output: ", outfile)
  message(strrep("=", 62))

  # Morris calls this with matrix X [n_runs × N_PARAMS]; we return a numeric vector.
  morris_fn <- function(X) {
    unlist(par_lapply(
      seq_len(nrow(X)),
      function(i) {
        x <- X[i, ]; names(x) <- PARAM_NAMES
        r <- run_one(x, base_groups, chl = chl, sst = sst, t_run = t_run)
        if (is.null(r)) return(0L)
        as.integer(r$fish_survive == 1L && r$fish_stable == 1L)
      },
      extra_export = c("base_groups", "chl", "sst", "t_run")
    ))
  }

  t0 <- proc.time()["elapsed"]
  m  <- sensitivity::morris(
    model   = morris_fn,
    factors = PARAM_NAMES,
    r       = n_morris,
    design  = list(type = "oat", levels = 6, grid.jump = 3),
    binf    = rep(0, N_PARAMS),
    bsup    = rep(1, N_PARAMS),
    scale   = FALSE
  )
  elapsed <- proc.time()["elapsed"] - t0
  message(sprintf("\nMorris complete in %.1f minutes.", elapsed / 60))

  mu_star    <- apply(m$ee, 2, function(x) mean(abs(x)))
  sigma_ee   <- apply(m$ee, 2, sd)
  summary_df <- data.frame(
    param   = PARAM_NAMES,
    mu_star = mu_star,
    sigma   = sigma_ee,
    group   = sub("_[^_]+$", "", PARAM_NAMES),
    pname   = sub("^[^_]+_[^_]+_", "", PARAM_NAMES),
    stringsAsFactors = FALSE
  ) |> dplyr::arrange(dplyr::desc(mu_star))

  out <- list(morris = m, summary = summary_df,
              n_morris = n_morris, n_runs = n_runs,
              elapsed_s = elapsed, timestamp = Sys.time())

  if (!is.null(outfile)) { saveRDS(out, outfile); message("Saved: ", outfile) }
  message("\nTop 10 parameters by |mu*|:")
  print(head(summary_df[, c("param", "mu_star", "sigma")], 10))
  invisible(out)
}

# =============================================================================
# STAGE 2: LHS SAMPLING
# =============================================================================
# Checkpoints every `checkpoint_every` rows — partial results survive a crash.
# Resume from checkpoint with resume_lhs_stage().

run_lhs_stage <- function(base_groups       = NULL,
                          n_lhs             = 500,
                          outfile           = "lhs_results.rds",
                          checkpoint_every  = 50,
                          chl_sweep         = FALSE,   # adds zoo_mae_chl; ~7× slower
                          sst               = SA_SST,
                          t_run             = SA_T_RUN) {

  if (is.null(base_groups)) { data(GroupInputs); base_groups <- GroupInputs }
  cfg     <- get0("PAR_CONFIG", envir = .GlobalEnv)
  n_cores <- if (is.null(cfg)) 1L else cfg$n_cores

  message("\n", strrep("=", 62))
  message("STAGE 2: LHS sampling")
  message(sprintf("  n_lhs=%d | chl_sweep=%s | cores=%d | checkpoint every %d",
                  n_lhs, chl_sweep, n_cores, checkpoint_every))
  if (!is.null(outfile)) message("  Output: ", outfile)
  message(strrep("=", 62))

  set.seed(2025)
  X <- lhs::randomLHS(n_lhs, N_PARAMS)
  colnames(X) <- PARAM_NAMES

  all_Y    <- vector("list", n_lhs)
  n_chunks <- ceiling(n_lhs / checkpoint_every)
  t_start  <- proc.time()["elapsed"]

  for (chunk in seq_len(n_chunks)) {
    idx <- ((chunk - 1L) * checkpoint_every + 1L) :
           min(chunk * checkpoint_every, n_lhs)

    message(sprintf("  chunk %d/%d  (rows %d–%d) ...",
                    chunk, n_chunks, min(idx), max(idx)))

    chunk_res <- par_lapply(
      idx,
      function(i) {
        x <- X[i, ]; names(x) <- PARAM_NAMES
        run_one(x, base_groups,
                chl = median(SA_CHL_GRADIENT), sst = sst,
                t_run = t_run, chl_sweep = chl_sweep)
      },
      extra_export = c("base_groups", "sst", "t_run", "chl_sweep")
    )

    for (j in seq_along(idx)) all_Y[[idx[j]]] <- chunk_res[[j]]

    if (!is.null(outfile))
      saveRDS(list(X = X, Y = all_Y, completed = max(idx),
                   n_lhs = n_lhs, timestamp = Sys.time()),
              outfile)

    elapsed <- (proc.time()["elapsed"] - t_start) / 60
    rate    <- max(idx) / max(elapsed, 0.001)
    eta     <- (n_lhs - max(idx)) / max(rate, 0.001)
    message(sprintf("    %.1f min elapsed | %.0f runs/min | ETA %.1f min",
                    elapsed, rate, eta))
  }

  Y_df <- results_to_df(all_Y)
  out  <- list(X = X, Y = Y_df, n_lhs = n_lhs,
               elapsed_s = proc.time()["elapsed"] - t_start,
               timestamp = Sys.time())
  if (!is.null(outfile)) { saveRDS(out, outfile); message("Final saved: ", outfile) }

  message(sprintf("Survival: %.1f%%  |  Stable coexistence: %.1f%%",
                  mean(Y_df$fish_survive, na.rm = TRUE) * 100,
                  mean(Y_df$fish_stable,  na.rm = TRUE) * 100))
  invisible(out)
}

# Resume LHS from last checkpoint
resume_lhs_stage <- function(checkpoint_file,
                             base_groups      = NULL,
                             checkpoint_every = 50,
                             chl_sweep        = FALSE,
                             sst              = SA_SST,
                             t_run            = SA_T_RUN) {

  if (!file.exists(checkpoint_file)) stop("File not found: ", checkpoint_file)
  if (is.null(base_groups)) { data(GroupInputs); base_groups <- GroupInputs }

  partial <- readRDS(checkpoint_file)
  start   <- partial$completed + 1L
  n_lhs   <- partial$n_lhs

  if (start > n_lhs) { message("Run already complete."); return(invisible(partial)) }
  message(sprintf("[resume] Resuming from row %d / %d", start, n_lhs))

  X     <- partial$X
  all_Y <- partial$Y
  n_chunks <- ceiling((n_lhs - start + 1L) / checkpoint_every)

  for (chunk in seq_len(n_chunks)) {
    idx <- seq(from = (chunk - 1L) * checkpoint_every + start,
               to   = min(chunk * checkpoint_every + start - 1L, n_lhs))
    chunk_res <- par_lapply(
      idx,
      function(i) {
        x <- X[i, ]; names(x) <- PARAM_NAMES
        run_one(x, base_groups, chl = median(SA_CHL_GRADIENT),
                sst = sst, t_run = t_run, chl_sweep = chl_sweep)
      },
      extra_export = c("base_groups", "sst", "t_run", "chl_sweep")
    )
    for (j in seq_along(idx)) all_Y[[idx[j]]] <- chunk_res[[j]]
    partial$Y <- all_Y; partial$completed <- max(idx); partial$timestamp <- Sys.time()
    saveRDS(partial, checkpoint_file)
    message(sprintf("  rows %d–%d saved.", min(idx), max(idx)))
  }

  Y_df <- results_to_df(all_Y)
  out  <- list(X = X, Y = Y_df, n_lhs = n_lhs, timestamp = Sys.time())
  saveRDS(out, checkpoint_file)
  invisible(out)
}

# =============================================================================
# STAGE 3: SOBOL VARIANCE DECOMPOSITION
# =============================================================================

run_sobol_stage <- function(lhs_file    = NULL,
                            lhs_result  = NULL,
                            base_groups = NULL,
                            top_k       = 10,
                            metric      = "fish_stable",
                            n_sobol     = 500,
                            outfile     = "sobol_results.rds",
                            sst         = SA_SST,
                            t_run       = SA_T_RUN) {

  if (is.null(base_groups)) { data(GroupInputs); base_groups <- GroupInputs }
  if (is.null(lhs_result) && !is.null(lhs_file)) lhs_result <- readRDS(lhs_file)
  if (is.null(lhs_result)) stop("Provide lhs_file or lhs_result.")

  top_params <- get_top_params_lhs(lhs_result, top_k, metric)
  np         <- length(top_params)
  cfg        <- get0("PAR_CONFIG", envir = .GlobalEnv)
  n_cores    <- if (is.null(cfg)) 1L else cfg$n_cores

  message("\n", strrep("=", 62))
  message(sprintf("STAGE 3: Sobol  top_k=%d | n_sobol=%d | total~=%d | cores=%d",
                  np, n_sobol, n_sobol * (np + 2L), n_cores))
  message(strrep("=", 62))

  sobol_fn <- function(X_sub) {
    unlist(par_lapply(
      seq_len(nrow(X_sub)),
      function(i) {
        x_full <- setNames(rep(0.5, N_PARAMS), PARAM_NAMES)
        x_full[top_params] <- X_sub[i, ]
        r <- run_one(x_full, base_groups, sst = sst, t_run = t_run)
        if (is.null(r)) return(NA_real_)
        as.numeric(r[[metric]])
      },
      extra_export = c("base_groups", "top_params", "metric", "sst", "t_run")
    ))
  }

  set.seed(2025)
  mk <- function() as.data.frame(matrix(runif(n_sobol * np), nrow = n_sobol,
                                        dimnames = list(NULL, top_params)))
  t0  <- proc.time()["elapsed"]
  sob <- sensitivity::sobol2002(model = sobol_fn, X1 = mk(), X2 = mk(), nboot = 100)
  message(sprintf("\nSobol complete in %.1f minutes.",
                  (proc.time()["elapsed"] - t0) / 60))

  out <- list(sobol = sob, top_params = top_params, metric = metric,
              n_sobol = n_sobol, timestamp = Sys.time())
  if (!is.null(outfile)) { saveRDS(out, outfile); message("Saved: ", outfile) }
  print(plot(sob))
  invisible(out)
}

# =============================================================================
# STAGE 4: OPTIMISATION
# =============================================================================

run_optimise_stage <- function(lhs_file    = NULL,
                               lhs_result  = NULL,
                               base_groups = NULL,
                               n_restarts  = 5,
                               outfile     = "optimal_params.rds",
                               sst         = SA_SST,
                               t_run       = SA_T_RUN) {

  if (is.null(base_groups)) { data(GroupInputs); base_groups <- GroupInputs }
  if (is.null(lhs_result) && !is.null(lhs_file)) lhs_result <- readRDS(lhs_file)
  if (is.null(lhs_result)) stop("Provide lhs_file or lhs_result.")

  message("\n", strrep("=", 62))
  message("STAGE 4: Optimisation  (n_restarts=", n_restarts, ")")
  message(strrep("=", 62))

  Y <- lhs_result$Y; X <- lhs_result$X
  valid <- !is.na(Y[, "fish_stable"]) & Y[, "fish_stable"] == 1L
  if (!any(valid)) {
    valid <- !is.na(Y[, "fish_survive"]) & Y[, "fish_survive"] == 1L
    message("NOTE: no stable runs found; seeding from surviving runs only.")
  }
  if (!any(valid)) stop("No surviving runs in LHS — widen parameter bounds or check model.")

  starts <- X[head(which(valid)[order(Y[valid, "fish_entropy"],
                                      decreasing = TRUE)], n_restarts), ,
              drop = FALSE]
  zoo_ref <- get0("ZOO_REF_COMP", envir = .GlobalEnv)

  objective <- function(x_unit) {
    x_unit <- pmax(0, pmin(1, x_unit)); names(x_unit) <- PARAM_NAMES
    r <- run_one(x_unit, base_groups, sst = sst, t_run = t_run,
                 chl_sweep = !is.null(zoo_ref))
    if (is.null(r) || r$fish_survive == 0L) return(100)
    stab_pen <- if (r$fish_stable == 0L) 5 else 0
    zoo_pen  <- if (!is.na(r$zoo_mae_chl)) 10 * r$zoo_mae_chl else 0
    -(r$fish_entropy) + stab_pen + zoo_pen
  }

  best <- list(value = Inf)
  for (s in seq_len(nrow(starts))) {
    message(sprintf("  restart %d/%d ...", s, nrow(starts)))
    res <- optim(starts[s, ], objective, method = "L-BFGS-B",
                 lower = rep(0, N_PARAMS), upper = rep(1, N_PARAMS),
                 control = list(maxit = 200, factr = 1e9))
    message(sprintf("    objective = %.4f", res$value))
    if (res$value < best$value) best <- res
  }

  best_params <- unit_to_params(setNames(best$par, PARAM_NAMES))
  out <- list(optim_result = best,
              params_unit  = setNames(best$par, PARAM_NAMES),
              params_real  = best_params,
              objective_value = best$value, timestamp = Sys.time())
  if (!is.null(outfile)) { saveRDS(out, outfile); message("Saved: ", outfile) }
  message("\nOptimal fish parameters:")
  print(lapply(best_params, function(p) round(unlist(p), 4)))
  invisible(out)
}

# =============================================================================
# HELPERS
# =============================================================================

get_top_params_morris <- function(morris_result, top_k = 10) {
  s <- head(morris_result$summary, top_k)
  message("Top ", top_k, " parameters by |mu*|:")
  print(s[, c("param", "mu_star", "sigma")])
  invisible(as.character(s$param))
}

get_top_params_lhs <- function(lhs_result, top_k = 10, metric = "fish_stable") {
  Y  <- lhs_result$Y; X <- lhs_result$X
  ok <- complete.cases(cbind(X, Y[, metric]))
  cors <- apply(X[ok, ], 2, function(x)
    abs(cor(x, Y[ok, metric], method = "spearman")))
  top <- names(sort(cors, decreasing = TRUE))[seq_len(min(top_k, N_PARAMS))]
  message("Top ", top_k, " params (Spearman |r| with ", metric, "):")
  print(round(sort(cors, decreasing = TRUE)[seq_len(min(top_k, N_PARAMS))], 3))
  invisible(top)
}

save_optimal_params <- function(base_groups, opt_result,
                                path = "optimal_fish_params.csv") {
  gi       <- patch_group_inputs(base_groups, opt_result$params_real)
  fish_idx <- which(gi$Type == "Fish")
  out      <- gi[fish_idx, c("Species", "Wmat", "PPMR", "FeedWidth",
                              "f_M", "K_growth", "repro_eff", "ZSpre", "ZSexp")]
  out$R_frac <- 1 - out$f_M - out$K_growth
  write.csv(out, path, row.names = FALSE)
  message("Saved: ", path); print(out); invisible(out)
}

# =============================================================================
# PLOTS
# =============================================================================

plot_morris <- function(morris_result) {
  df        <- morris_result$summary
  df$label  <- sub("^[^_]+_[^_]+_", "", df$param)
  ggplot(df, aes(x = mu_star, y = sigma, label = label, colour = group)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20) +
    geom_abline(slope = c(0.5, 1, 2), linetype = "dashed", colour = "grey60") +
    labs(title    = "Morris Elementary Effects — ZooMSS fish parameters",
         subtitle = "Dashed: sigma = 0.5, 1, 2 × mu*",
         x = "Mean |EE|  (mu*)", y = "SD of EE  (sigma)", colour = "Fish group") +
    theme_bw(base_size = 12)
}

plot_lhs <- function(lhs_result, top_k = 6) {
  top <- get_top_params_lhs(lhs_result, top_k)
  df  <- cbind(as.data.frame(lhs_result$X[, top, drop = FALSE]),
               as.data.frame(lhs_result$Y[, c("fish_stable", "fish_entropy")]))
  plots <- lapply(top, function(p) {
    ggplot(df, aes_string(x = p, y = "fish_entropy",
                          colour = "factor(fish_stable)")) +
      geom_point(alpha = 0.4, size = 1) +
      scale_colour_manual(values = c("0" = "#d73027", "1" = "#1a9850"),
                          name = "Stable") +
      labs(x = sub("^[^_]+_[^_]+_", "", p), y = "Entropy") +
      theme_bw(base_size = 10)
  })
  patchwork::wrap_plots(plots, ncol = 3)
}

plot_zoo_gradient <- function(base_groups, opt_result = NULL,
                              chl_values = SA_CHL_GRADIENT,
                              sst = SA_SST, t_run = SA_T_RUN) {
  zoo_idx   <- which(base_groups$Type == "Zooplankton")
  zoo_names <- base_groups$Species[zoo_idx]
  zoo_cols  <- setNames(base_groups$PlotColour[zoo_idx], zoo_names)

  run_zoo <- function(gi, label) {
    dplyr::bind_rows(lapply(seq_along(chl_values), function(ci) {
      bt  <- run_model_once(gi, chl = chl_values[ci], sst = sst, t_run = t_run)
      bio <- bt[nrow(bt), zoo_idx]
      data.frame(chl   = log10(chl_values[ci]),
                 group = zoo_names,
                 prop  = bio / (sum(bio) + 1e-30),
                 run   = label)
    }))
  }

  df <- run_zoo(base_groups, "Default")
  if (!is.null(opt_result))
    df <- dplyr::bind_rows(df,
      run_zoo(patch_group_inputs(base_groups, opt_result$params_real), "Optimal"))

  ggplot(df, aes(x = chl, y = prop, fill = group)) +
    geom_area(position = "stack") +
    scale_fill_manual(values = zoo_cols) +
    facet_wrap(~run) +
    labs(x = "log10(Chl) [mg m⁻³]", y = "Biomass proportion",
         title = "Zooplankton composition along chl gradient") +
    theme_bw(base_size = 12)
}
