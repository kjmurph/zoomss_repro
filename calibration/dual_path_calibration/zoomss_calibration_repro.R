# =============================================================================
# Fish Reproduction Parameter Calibration Framework for ZooMSS
# =============================================================================
#
# Calibrates dual-pathway fish energy budget parameters to achieve:
#   1. Coexistence: All 3 fish groups persist across chl gradient
#   2. Stability: Oscillating steady-state (bounded CV)
#   3. Zoo conservation: Zooplankton composition matches legacy model
#   4. Biological realism: Parameters within literature-defensible ranges
#
# Architecture: LHS exploration -> filter -> L-BFGS-B refinement
# Follows pattern from zoomss_calibration.R (fishing mortality calibration)
#
# All parameters are group-specific (S = Small, M = Medium, L = Large):
#   PPMR (3), FeedWidth (3), K_growth (3), f_M (3), repro_eff (3),
#   Wmat (3), ZSpre (3), ZSexp (3)
#   Total free dimensions: 24
#
# Constraints:
#   Energy budget:  f_M_X + K_growth_X + R_frac_X = 1, R_frac_X >= 0.15
#                   (enforced independently per fish group)
#   Wmat ordering:  Wmat_S <= Wmat_M <= Wmat_L (biological realism)
#
# Rounding: PPMR rounded to integer, all other parameters to 2 d.p.
#
# Environmental forcing: Seasonal SST and chl via createEnviroData().
# =============================================================================


# --- Parameter Space Definition -----------------------------------------------

#' Define the calibration parameter space
#'
#' All 8 parameter types are group-specific (S, M, L).
#' Total: 24 free dimensions.
#'
#' @return data.frame with columns: name, lower, upper, group_idx
#' @export
repro_param_space <- function() {

  # Parameter types with bounds: c(lower, upper)
  # Wmat has group-specific bounds (different Wmax per group)
  param_defs <- list(
    PPMR      = c(50,    500),
    FeedWidth = c(0.8,   2.0),
    K_growth  = c(0.15,  0.45),
    f_M       = c(0.30,  0.70),
    repro_eff = c(0.001, 1.0),
    Wmat      = list(S = c(-2.0, 1.0), M = c(-2.0, 3.0), L = c(-2.0, 5.0)),
    ZSpre     = c(0.01,  1.0),
    ZSexp     = c(0.1,   1.0)
  )

  rows <- list()
  suffixes <- c("S", "M", "L")
  grp_idx  <- c(1L, 2L, 3L)

  for (pname in names(param_defs)) {
    bounds <- param_defs[[pname]]
    for (i in seq_along(suffixes)) {
      sfx <- suffixes[i]
      full_name <- paste0(pname, "_", sfx)
      if (is.list(bounds)) {
        b <- bounds[[sfx]]
      } else {
        b <- bounds
      }
      rows[[length(rows) + 1]] <- data.frame(
        name = full_name, lower = b[1], upper = b[2],
        group_idx = grp_idx[i], stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}


#' Round parameters for model application
#'
#' PPMR values are rounded to integers; all other parameters to 2 d.p.
#'
#' @param par Named numeric vector
#' @return Named numeric vector with rounded values
#' @export
round_params <- function(par) {
  ppmr_idx <- grep("^PPMR", names(par))
  if (length(ppmr_idx) > 0) par[ppmr_idx] <- round(par[ppmr_idx])
  other_idx <- setdiff(seq_along(par), ppmr_idx)
  if (length(other_idx) > 0) par[other_idx] <- round(par[other_idx], 2)
  par
}


#' Map a parameter vector to Groups data.frame modifications
#'
#' All parameters are group-specific. Rounded before application.
#'
#' @param par Named numeric vector of parameter values
#' @param Groups data.frame from getGroups()
#' @return Modified Groups data.frame
#' @export
apply_repro_params <- function(par, Groups) {

  fish_idx <- which(Groups$Type == "Fish")
  stopifnot(length(fish_idx) == 3)

  par <- round_params(par)

  grp_map <- c(S = 1L, M = 2L, L = 3L)
  param_types <- c("PPMR", "FeedWidth", "K_growth", "f_M",
                    "repro_eff", "Wmat", "ZSpre", "ZSexp")

  for (suffix in c("S", "M", "L")) {
    fi <- fish_idx[grp_map[suffix]]
    for (ptype in param_types) {
      pname <- paste0(ptype, "_", suffix)
      if (pname %in% names(par)) {
        Groups[[ptype]][fi] <- par[pname]
      }
    }
  }

  Groups$repro_on[fish_idx] <- 1L
  Groups
}


#' Enforce energy budget constraint per fish group
#'
#' Each group must independently satisfy: R_frac = 1 - f_M - K_growth >= 0.15
#'
#' @param par Named numeric vector
#' @param min_rfrac Minimum R_frac (default 0.15)
#' @param tol Floating-point tolerance (default 1e-10)
#' @return Logical: TRUE if all three groups satisfy constraint
#' @export
check_energy_constraint <- function(par, min_rfrac = 0.15, tol = 1e-10) {
  for (sfx in c("S", "M", "L")) {
    f_M <- par[paste0("f_M_", sfx)]
    K_growth <- par[paste0("K_growth_", sfx)]
    R_frac <- 1 - f_M - K_growth
    if (R_frac < (min_rfrac - tol)) return(FALSE)
  }
  TRUE
}


#' Enforce maturation size ordering: Wmat_S <= Wmat_M <= Wmat_L
#'
#' @param par Named numeric vector
#' @return Logical: TRUE if ordering constraint satisfied
#' @export
check_wmat_ordering <- function(par) {
  par["Wmat_S"] <= par["Wmat_M"] && par["Wmat_M"] <= par["Wmat_L"]
}


# --- Seasonal Environment Helper ---------------------------------------------

#' Create seasonal environmental forcing for a given chl level
#'
#' @param n_years Simulation length in years
#' @param dt Time step (default 0.1)
#' @param base_sst Base SST (default 15)
#' @param base_chl Base chlorophyll concentration (mg/m^3)
#' @param sst_amplitude SST seasonal amplitude (default 4)
#' @param chl_amplitude Chl seasonal amplitude (default 0.5)
#' @return Input parameters list for zoomss_model()
#' @export
create_seasonal_input <- function(n_years, dt = 0.1,
                                  base_sst = 15, base_chl = 1.0,
                                  sst_amplitude = 4, chl_amplitude = 0.5) {

  env_data <- createEnviroData(
    n_years       = n_years,
    dt            = dt,
    base_sst      = base_sst,
    base_chl      = base_chl,
    seasonal      = TRUE,
    sst_amplitude = sst_amplitude,
    chl_amplitude = chl_amplitude
  )

  createInputParams(
    time = env_data$time,
    sst  = env_data$sst,
    chl  = env_data$chl
  )
}


# --- Benchmark Generation -----------------------------------------------------

#' Generate legacy benchmark (repro_on = 0) across chl gradient
#'
#' @param chl_levels Numeric vector of chlorophyll concentrations (mg/m^3)
#' @param sst Numeric, base sea surface temperature (default 15)
#' @param n_years Numeric, simulation length in years (default 400)
#' @param dt Numeric, time step (default 0.1)
#' @param sst_amplitude SST seasonal amplitude (default 4)
#' @param chl_amplitude Chl seasonal amplitude (default 0.5)
#' @param assess_years Final N years for metric extraction (default 100)
#' @param cache_dir Character, directory for caching results
#' @param n_workers Integer, number of parallel workers (default 14)
#' @param force_rerun Logical, re-run even if cache exists
#' @return List with zoo_proportions, fish_biomass, and metadata
#' @export
generate_legacy_benchmark <- function(chl_levels,
                                      sst = 15,
                                      n_years = 400,
                                      dt = 0.1,
                                      sst_amplitude = 4,
                                      chl_amplitude = 0.5,
                                      assess_years = 100,
                                      cache_dir = NULL,
                                      n_workers = 14,
                                      force_rerun = FALSE) {

  if (is.null(cache_dir)) {
    cache_dir <- file.path(tempdir(), "zoomss_calib_benchmark")
  }
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  benchmark_file <- file.path(cache_dir, sprintf("benchmark_sst%.0f.rds", sst))
  if (file.exists(benchmark_file) && !force_rerun) {
    message("Loading cached benchmark: ", benchmark_file)
    return(readRDS(benchmark_file))
  }

  message("Generating legacy benchmark: ", length(chl_levels),
          " chl levels at SST = ", sst, "\u00b0C (seasonal, ",
          n_years, "yr, assess final ", assess_years, "yr)")

  Groups <- getGroups()
  fish_idx <- which(Groups$Type == "Fish")
  Groups$repro_on[fish_idx] <- 0L

  future::plan(future::multisession, workers = min(n_workers, length(chl_levels)))
  on.exit(future::plan(future::sequential), add = TRUE)

  results <- furrr::future_map(seq_along(chl_levels), function(i) {
    devtools::load_all(quiet = TRUE)
    chl <- chl_levels[i]
    cache_file <- file.path(cache_dir,
                            sprintf("legacy_sst%.0f_chl%.4f.rds", sst, log10(chl)))
    if (file.exists(cache_file) && !force_rerun) return(readRDS(cache_file))

    input_params <- create_seasonal_input(
      n_years = n_years, dt = dt,
      base_sst = sst, base_chl = chl,
      sst_amplitude = sst_amplitude, chl_amplitude = chl_amplitude
    )
    mdl <- zoomss_model(input_params = input_params, Groups = Groups, isave = 2)
    saveRDS(mdl, cache_file)
    mdl
  }, .options = furrr::furrr_options(seed = TRUE), .progress = TRUE)

  n_chl <- length(chl_levels)
  zoo_names <- Groups$Species[Groups$Type == "Zooplankton"]
  fish_names <- Groups$Species[Groups$Type == "Fish"]

  zoo_proportions <- matrix(NA, nrow = n_chl, ncol = length(zoo_names),
                            dimnames = list(NULL, zoo_names))
  fish_biomass <- matrix(NA, nrow = n_chl, ncol = length(fish_names),
                         dimnames = list(NULL, fish_names))

  for (i in seq_along(results)) {
    mdl <- results[[i]]
    avg <- averageTimeSeries(mdl, var = "biomass", n_years = assess_years)
    avg_total <- rowSums(avg)
    zoo_idx <- which(mdl$param$Groups$Type == "Zooplankton")
    zoo_bm <- avg_total[zoo_idx]
    zoo_total <- sum(zoo_bm)
    if (zoo_total > 0) zoo_proportions[i, ] <- zoo_bm / zoo_total
    fish_biomass[i, ] <- avg_total[mdl$param$fish_grps]
  }

  benchmark <- list(
    sst = sst, chl_levels = chl_levels, log10_chl = log10(chl_levels),
    zoo_proportions = zoo_proportions, fish_biomass = fish_biomass,
    zoo_names = zoo_names, fish_names = fish_names,
    n_years = n_years, dt = dt, assess_years = assess_years,
    sst_amplitude = sst_amplitude, chl_amplitude = chl_amplitude
  )
  saveRDS(benchmark, benchmark_file)
  message("Benchmark saved: ", benchmark_file)
  benchmark
}


# --- Objective Function -------------------------------------------------------

#' Evaluate a parameter set across multiple chlorophyll levels
#'
#' Runs the model at each chl level with seasonal forcing and computes a
#' composite score (lower is better). Assessment uses the final assess_years.
#' Enforces energy budget (per group) and Wmat ordering constraints.
#'
#' @param par Named numeric vector of parameter values
#' @param benchmark Output of generate_legacy_benchmark()
#' @param chl_indices Integer vector of indices to evaluate
#' @param n_years Simulation length (default 300)
#' @param dt Time step (default 0.1)
#' @param assess_years Final N years for metric extraction (default 100)
#' @param weights Named list of objective weights
#' @param return_details If TRUE, return detailed diagnostics
#' @return Numeric scalar or list if return_details = TRUE
#' @export
repro_objective <- function(par,
                            benchmark,
                            chl_indices = NULL,
                            n_years = 300,
                            dt = 0.1,
                            assess_years = 100,
                            weights = NULL,
                            return_details = FALSE) {

  if (is.null(weights)) {
    weights <- list(
      coexistence     = 5.0,
      stability       = 1.0,
      zoo_composition = 3.0,
      fish_ratio      = 0.0,
      spectrum_slope  = 1.0
    )
  }

  if (is.null(chl_indices)) chl_indices <- seq_along(benchmark$chl_levels)

  # Constraint checks
  if (!check_energy_constraint(par)) {
    if (return_details) return(list(score = 1e6, reason = "energy_constraint_violated"))
    return(1e6)
  }
  if (!check_wmat_ordering(par)) {
    if (return_details) return(list(score = 1e6, reason = "wmat_ordering_violated"))
    return(1e6)
  }

  Groups <- getGroups()
  Groups <- apply_repro_params(par, Groups)

  sst_amp <- if (!is.null(benchmark$sst_amplitude)) benchmark$sst_amplitude else 4
  chl_amp <- if (!is.null(benchmark$chl_amplitude)) benchmark$chl_amplitude else 0.5

  n_eval <- length(chl_indices)
  scores <- data.frame(
    chl_idx     = chl_indices,
    log10_chl   = benchmark$log10_chl[chl_indices],
    coexistence = rep(NA_real_, n_eval),
    stability   = rep(NA_real_, n_eval),
    zoo_comp    = rep(NA_real_, n_eval),
    fish_ratio  = rep(NA_real_, n_eval),
    spectrum    = rep(NA_real_, n_eval)
  )

  for (j in seq_along(chl_indices)) {
    ci <- chl_indices[j]
    chl <- benchmark$chl_levels[ci]

    tryCatch({
      input_params <- create_seasonal_input(
        n_years = n_years, dt = dt,
        base_sst = benchmark$sst, base_chl = chl,
        sst_amplitude = sst_amp, chl_amplitude = chl_amp
      )
      mdl <- zoomss_model(input_params = input_params, Groups = Groups, isave = 2)

      n_save <- length(mdl$time)
      dt_save <- mdl$param$isave * mdl$param$dt
      n_assess <- min(n_save, round(assess_years / dt_save))
      start_idx <- max(1, n_save - n_assess + 1)

      fish_grps <- mdl$param$fish_grps
      zoo_idx <- which(mdl$param$Groups$Type == "Zooplankton")
      n_fish <- length(fish_grps)

      # 1. COEXISTENCE
      fish_mean_bm <- sapply(seq_len(n_fish), function(f) {
        bm_slice <- mdl$biomass[start_idx:n_save, fish_grps[f], , drop = FALSE]
        mean(rowSums(bm_slice, dims = 2), na.rm = TRUE)
      })
      coexist_frac <- mean(fish_mean_bm > 1e-20)
      scores$coexistence[j] <- 1 - coexist_frac

      # 2. STABILITY (CV bounded)
      fish_cv <- sapply(seq_len(n_fish), function(f) {
        bm_slice <- mdl$biomass[start_idx:n_save, fish_grps[f], , drop = FALSE]
        bm <- rowSums(bm_slice, dims = 2)
        bm <- bm[bm > 0]
        if (length(bm) < 10) return(Inf)
        sd(bm) / mean(bm)
      })
      cv_penalty <- sapply(fish_cv, function(cv) {
        if (is.infinite(cv) || is.na(cv)) return(1)
        if (cv > 2.0) return(min(1, (cv - 2.0) / 5))
        if (cv < 0.001) return(0.1)
        return(0)
      })
      scores$stability[j] <- mean(cv_penalty)

      # 3. ZOO COMPOSITION (correlation with legacy)
      avg_bm <- averageTimeSeries(mdl, var = "biomass", n_years = assess_years)
      avg_bm_total <- rowSums(avg_bm)
      zoo_bm <- avg_bm_total[zoo_idx]
      zoo_total <- sum(zoo_bm)
      if (zoo_total > 0) {
        zoo_prop <- zoo_bm / zoo_total
        bench_prop <- benchmark$zoo_proportions[ci, ]
        if (all(!is.na(bench_prop)) && all(!is.na(zoo_prop))) {
          r <- cor(zoo_prop, bench_prop, method = "pearson")
          scores$zoo_comp[j] <- max(0, 1 - r)
        } else {
          scores$zoo_comp[j] <- 1
        }
      } else {
        scores$zoo_comp[j] <- 1
      }

      # 4. FISH BIOMASS RATIO
      bench_fish <- benchmark$fish_biomass[ci, ]
      model_fish <- avg_bm_total[fish_grps]
      if (all(bench_fish > 0) && all(model_fish > 0)) {
        log_ratios <- log10(model_fish / bench_fish)
        ratio_penalty <- sapply(log_ratios, function(lr) {
          if (abs(lr) <= 1) return(0)
          return(min(1, (abs(lr) - 1) / 2))
        })
        scores$fish_ratio[j] <- mean(ratio_penalty)
      } else {
        scores$fish_ratio[j] <- 0.5
      }

      # 5. SIZE SPECTRUM SLOPE
      all_abundance <- averageTimeSeries(mdl, var = "abundance",
                                         n_years = assess_years)
      w_vec <- mdl$param$w
      total_abund <- colSums(all_abundance)
      valid <- total_abund > 0
      if (sum(valid) > 5) {
        fit <- lm(log10(total_abund[valid]) ~ w_vec[valid])
        slope <- coef(fit)[2]
        if (slope > 0) {
          scores$spectrum[j] <- 1
        } else if (slope < -2) {
          scores$spectrum[j] <- min(1, (abs(slope) - 2) / 3)
        } else {
          scores$spectrum[j] <- 0
        }
      } else {
        scores$spectrum[j] <- 0.5
      }

    }, error = function(e) {
      scores$coexistence[j] <<- 1
      scores$stability[j]   <<- 1
      scores$zoo_comp[j]    <<- 1
      scores$fish_ratio[j]  <<- 1
      scores$spectrum[j]    <<- 1
    })
  }

  metric_scores <- c(
    coexistence = mean(scores$coexistence, na.rm = TRUE),
    stability   = mean(scores$stability, na.rm = TRUE),
    zoo_comp    = mean(scores$zoo_comp, na.rm = TRUE),
    fish_ratio  = mean(scores$fish_ratio, na.rm = TRUE),
    spectrum    = mean(scores$spectrum, na.rm = TRUE)
  )
  w <- c(weights$coexistence, weights$stability, weights$zoo_composition,
         weights$fish_ratio, weights$spectrum_slope)
  composite <- sum(metric_scores * w) / sum(w)

  if (return_details) {
    return(list(score = composite, metric_scores = metric_scores,
                per_chl_scores = scores, par = par))
  }
  composite
}


# --- LHS Exploration ----------------------------------------------------------

#' Generate Latin Hypercube Sample of parameter space
#'
#' Enforces per-group energy budget (R_frac >= 0.15) and Wmat ordering
#' (Wmat_S <= Wmat_M <= Wmat_L via sort).
#'
#' @param n_samples Number of samples (default 500)
#' @param param_space Output of repro_param_space()
#' @param seed Random seed
#' @return data.frame with one column per parameter
#' @export
generate_lhs_samples <- function(n_samples = 500,
                                 param_space = repro_param_space(),
                                 seed = 42) {
  set.seed(seed)
  n_params <- nrow(param_space)
  lhs_unit <- lhs::randomLHS(n_samples, n_params)
  colnames(lhs_unit) <- param_space$name

  lhs_scaled <- lhs_unit
  for (i in seq_len(n_params)) {
    lhs_scaled[, i] <- param_space$lower[i] +
      lhs_unit[, i] * (param_space$upper[i] - param_space$lower[i])
  }

  # Enforce energy constraint per group via rejection resampling
  for (row in seq_len(n_samples)) {
    for (sfx in c("S", "M", "L")) {
      f_M_col <- paste0("f_M_", sfx)
      K_col   <- paste0("K_growth_", sfx)
      f_M_idx <- which(param_space$name == f_M_col)
      K_idx   <- which(param_space$name == K_col)

      attempts <- 0
      while ((1 - lhs_scaled[row, f_M_idx] - lhs_scaled[row, K_idx]) < 0.15 &&
             attempts < 100) {
        lhs_scaled[row, f_M_idx] <- runif(1, param_space$lower[f_M_idx],
                                            param_space$upper[f_M_idx])
        lhs_scaled[row, K_idx]   <- runif(1, param_space$lower[K_idx],
                                            param_space$upper[K_idx])
        attempts <- attempts + 1
      }
      # Fallback: force compliance
      if ((1 - lhs_scaled[row, f_M_idx] - lhs_scaled[row, K_idx]) < 0.15) {
        lhs_scaled[row, K_idx] <- min(lhs_scaled[row, K_idx],
                                       0.85 - lhs_scaled[row, f_M_idx])
      }
    }
  }

  # Enforce Wmat ordering: sort so Wmat_S <= Wmat_M <= Wmat_L
  wmat_cols <- match(c("Wmat_S", "Wmat_M", "Wmat_L"), param_space$name)
  for (row in seq_len(n_samples)) {
    wmat_vals <- sort(lhs_scaled[row, wmat_cols])
    lhs_scaled[row, wmat_cols] <- wmat_vals
  }

  as.data.frame(lhs_scaled)
}


#' Run LHS exploration in parallel with restart capability
#'
#' @param lhs_samples data.frame from generate_lhs_samples()
#' @param benchmark Output of generate_legacy_benchmark()
#' @param chl_indices Indices of representative chl levels
#' @param n_years Simulation length for screening (default 300)
#' @param n_workers Number of parallel workers (default 14)
#' @param cache_dir Directory for caching results
#' @param batch_size Samples per batch (default 50)
#' @return data.frame with parameters and scores
#' @export
run_lhs_exploration <- function(lhs_samples,
                                benchmark,
                                chl_indices = NULL,
                                n_years = 300,
                                n_workers = 14,
                                cache_dir = NULL,
                                batch_size = 50) {

  if (is.null(cache_dir)) cache_dir <- file.path(tempdir(), "zoomss_calib_lhs")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  if (is.null(chl_indices)) {
    target_log10 <- c(-1.5, -1.0, -0.5, 0.0, 0.4)
    chl_indices <- sapply(target_log10, function(t) {
      which.min(abs(benchmark$log10_chl - t))
    })
    chl_indices <- unique(chl_indices)
  }

  n_total <- nrow(lhs_samples)
  results_file <- file.path(cache_dir, "lhs_results.rds")

  if (file.exists(results_file)) {
    existing <- readRDS(results_file)
    start_idx <- nrow(existing) + 1
    message("Resuming from sample ", start_idx, " of ", n_total)
  } else {
    existing <- NULL
    start_idx <- 1
  }

  if (start_idx > n_total) {
    message("All samples already evaluated.")
    return(existing)
  }

  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  all_results <- existing
  n_batches <- ceiling((n_total - start_idx + 1) / batch_size)

  for (batch in seq_len(n_batches)) {
    batch_start <- start_idx + (batch - 1) * batch_size
    batch_end <- min(batch_start + batch_size - 1, n_total)
    batch_idx <- batch_start:batch_end

    message(sprintf("Batch %d/%d: samples %d-%d",
                    batch, n_batches, batch_start, batch_end))

    batch_results <- furrr::future_map_dfr(batch_idx, function(i) {
      devtools::load_all(quiet = TRUE)
      par <- as.numeric(lhs_samples[i, ])
      names(par) <- names(lhs_samples)

      result <- repro_objective(
        par = par, benchmark = benchmark,
        chl_indices = chl_indices, n_years = n_years,
        return_details = TRUE
      )

      data.frame(
        sample_id   = i,
        score       = result$score,
        coexistence = result$metric_scores["coexistence"],
        stability   = result$metric_scores["stability"],
        zoo_comp    = result$metric_scores["zoo_comp"],
        fish_ratio  = result$metric_scores["fish_ratio"],
        spectrum    = result$metric_scores["spectrum"],
        t(par), stringsAsFactors = FALSE
      )
    }, .options = furrr::furrr_options(seed = TRUE), .progress = TRUE)

    all_results <- rbind(all_results, batch_results)
    saveRDS(all_results, results_file)
    message(sprintf("  Checkpoint: %d/%d complete", batch_end, n_total))
  }
  all_results
}


# --- Filtering and Refinement -------------------------------------------------

#' Filter LHS results by hard constraints
#'
#' @param lhs_results data.frame from run_lhs_exploration()
#' @param max_coexistence Max coexistence penalty (0 = all coexist)
#' @param max_zoo_comp Max zoo composition penalty
#' @param max_score Max composite score
#' @param top_n Number of top candidates to return
#' @return Filtered and sorted data.frame
#' @export
filter_lhs_candidates <- function(lhs_results,
                                  max_coexistence = 0.01,
                                  max_zoo_comp = 0.3,
                                  max_score = NULL,
                                  top_n = 20) {

  candidates <- lhs_results[
    lhs_results$coexistence <= max_coexistence &
    lhs_results$zoo_comp <= max_zoo_comp, ]

  if (!is.null(max_score)) candidates <- candidates[candidates$score <= max_score, ]
  candidates <- candidates[order(candidates$score), ]
  if (nrow(candidates) > top_n) candidates <- candidates[1:top_n, ]

  message(sprintf("Filtered: %d candidates from %d samples",
                  nrow(candidates), nrow(lhs_results)))
  candidates
}


#' Refine a candidate with L-BFGS-B optimisation
#'
#' @param par_init Named numeric vector (starting point)
#' @param benchmark Output of generate_legacy_benchmark()
#' @param chl_indices Indices for evaluation (NULL = all)
#' @param n_years Simulation length (default 400)
#' @param param_space Output of repro_param_space()
#' @param maxit Maximum L-BFGS-B iterations (default 50)
#' @return List with optimised parameters and diagnostics
#' @export
refine_candidate <- function(par_init,
                             benchmark,
                             chl_indices = NULL,
                             n_years = 400,
                             param_space = repro_param_space(),
                             maxit = 50) {

  if (is.null(chl_indices)) chl_indices <- seq_along(benchmark$chl_levels)

  lower <- setNames(param_space$lower, param_space$name)
  upper <- setNames(param_space$upper, param_space$name)

  obj_fn <- function(par) {
    names(par) <- names(par_init)
    if (!check_energy_constraint(par)) return(1e6)
    if (!check_wmat_ordering(par)) return(1e6)
    repro_objective(par = par, benchmark = benchmark,
                    chl_indices = chl_indices, n_years = n_years)
  }

  result <- optim(
    par = par_init, fn = obj_fn, method = "L-BFGS-B",
    lower = lower[names(par_init)], upper = upper[names(par_init)],
    control = list(maxit = maxit, trace = 1)
  )

  names(result$par) <- names(par_init)
  result$par <- round_params(result$par)

  details <- repro_objective(
    par = result$par, benchmark = benchmark,
    chl_indices = chl_indices, n_years = n_years,
    return_details = TRUE
  )

  list(par = result$par, score = result$value,
       convergence = result$convergence, details = details,
       optim_result = result)
}


# --- Yield Curve Validation (Post-hoc) ----------------------------------------

#' Generate yield curves for a calibrated parameter set
#'
#' @param par Named numeric vector of calibrated parameters
#' @param fmort_levels Numeric vector of fishing mortality rates
#' @param chl Chlorophyll (mg/m^3, default 1.0)
#' @param sst Base SST (default 15)
#' @param n_years Simulation length (default 400)
#' @param dt Time step (default 0.1)
#' @param sst_amplitude SST seasonal amplitude (default 4)
#' @param chl_amplitude Chl seasonal amplitude (default 0.5)
#' @param assess_years Final N years for assessment (default 100)
#' @return data.frame with Fmort, Fish_Group, Biomass, Yield
#' @export
yield_curve_validation <- function(par,
                                   fmort_levels = seq(0, 2, by = 0.1),
                                   chl = 1.0, sst = 15,
                                   n_years = 400, dt = 0.1,
                                   sst_amplitude = 4,
                                   chl_amplitude = 0.5,
                                   assess_years = 100) {

  Groups <- getGroups()
  Groups <- apply_repro_params(par, Groups)
  fish_idx <- which(Groups$Type == "Fish")
  fish_names <- Groups$Species[fish_idx]
  results <- data.frame()

  for (fm in fmort_levels) {
    input_params <- create_seasonal_input(
      n_years = n_years, dt = dt,
      base_sst = sst, base_chl = chl,
      sst_amplitude = sst_amplitude, chl_amplitude = chl_amplitude
    )
    mdl_Groups <- Groups
    if ("Fmort" %in% names(mdl_Groups)) mdl_Groups$Fmort[fish_idx] <- fm

    tryCatch({
      mdl <- zoomss_model(input_params = input_params, Groups = mdl_Groups, isave = 2)
      avg_bm <- averageTimeSeries(mdl, var = "biomass", n_years = assess_years)
      avg_bm_total <- rowSums(avg_bm)
      for (f in seq_along(fish_names)) {
        results <- rbind(results, data.frame(
          Fmort = fm, Fish_Group = fish_names[f],
          Biomass = avg_bm_total[fish_idx[f]],
          Yield = fm * avg_bm_total[fish_idx[f]],
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      for (f in seq_along(fish_names)) {
        results <<- rbind(results, data.frame(
          Fmort = fm, Fish_Group = fish_names[f],
          Biomass = NA, Yield = NA, stringsAsFactors = FALSE
        ))
      }
    })
  }
  results
}


# --- Pipeline Wrapper ---------------------------------------------------------

#' Run the complete calibration pipeline
#'
#' @param n_samples LHS samples (default 500)
#' @param n_workers Parallel workers (default 14)
#' @param cache_dir Base cache directory
#' @param sst Base temperature (default 15)
#' @param sst_amplitude SST seasonal amplitude (default 4)
#' @param chl_amplitude Chl seasonal amplitude (default 0.5)
#' @param benchmark_years Benchmark simulation length (default 400)
#' @param screening_years LHS screening run length (default 300)
#' @param refinement_years Refinement run length (default 400)
#' @param assess_years Final years for metric extraction (default 100)
#' @param top_n_refine Candidates to refine (default 5)
#' @param seed Random seed
#' @return List with all calibration results
#' @export
run_repro_calibration <- function(n_samples = 500,
                                  n_workers = 14,
                                  cache_dir = "calibration_repro_cache",
                                  sst = 15,
                                  sst_amplitude = 4,
                                  chl_amplitude = 0.5,
                                  benchmark_years = 400,
                                  screening_years = 300,
                                  refinement_years = 400,
                                  assess_years = 100,
                                  top_n_refine = 5,
                                  seed = 42) {

  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  message("=== Phase 1: Generating Legacy Benchmark ===")
  log10_chl_seq <- seq(-1.7, 0.5, by = 0.1)
  chl_levels <- 10^log10_chl_seq

  benchmark <- generate_legacy_benchmark(
    chl_levels = chl_levels, sst = sst,
    n_years = benchmark_years, assess_years = assess_years,
    sst_amplitude = sst_amplitude, chl_amplitude = chl_amplitude,
    n_workers = n_workers,
    cache_dir = file.path(cache_dir, "benchmark")
  )

  message("\n=== Phase 2: LHS Parameter Exploration ===")
  lhs_samples <- generate_lhs_samples(n_samples = n_samples, seed = seed)
  lhs_results <- run_lhs_exploration(
    lhs_samples = lhs_samples, benchmark = benchmark,
    n_years = screening_years, n_workers = n_workers,
    cache_dir = file.path(cache_dir, "lhs")
  )

  message("\n=== Phase 2b: Filtering Candidates ===")
  candidates <- filter_lhs_candidates(lhs_results, top_n = top_n_refine * 2)

  message("\n=== Phase 3: Refining Top Candidates ===")
  n_refine <- min(top_n_refine, nrow(candidates))
  refined <- list()
  param_names <- names(lhs_samples)

  for (k in seq_len(n_refine)) {
    message(sprintf("\nRefining candidate %d/%d (LHS score: %.4f)",
                    k, n_refine, candidates$score[k]))
    par_init <- as.numeric(candidates[k, param_names])
    names(par_init) <- param_names
    refined[[k]] <- refine_candidate(
      par_init = par_init, benchmark = benchmark,
      n_years = refinement_years
    )
  }

  calibration <- list(
    benchmark = benchmark,
    lhs_samples = lhs_samples,
    lhs_results = lhs_results,
    candidates = candidates,
    refined = refined,
    best = refined[[which.min(sapply(refined, function(x) x$score))]],
    param_space = repro_param_space(),
    settings = list(n_samples = n_samples, sst = sst,
                    sst_amplitude = sst_amplitude,
                    chl_amplitude = chl_amplitude,
                    benchmark_years = benchmark_years,
                    screening_years = screening_years,
                    refinement_years = refinement_years,
                    assess_years = assess_years, seed = seed)
  )

  saveRDS(calibration, file.path(cache_dir, "calibration_results.rds"))
  message("\n=== Calibration Complete ===")
  message(sprintf("Best score: %.4f", calibration$best$score))
  calibration
}
