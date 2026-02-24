# =============================================================================
# ZooMSS Fishing Calibration Framework
# =============================================================================
# Calibrate catchability (q) and selectivity threshold (w_min) parameters
# by fitting predicted catch to observed catch data.
#
# Workflow (following DBPM pattern):
#   1. Define parameter bounds per fish group
#   2. Latin Hypercube Sampling (LHS) across parameter space
#   3. Evaluate each sample: correlation + RMSE against observed catch
#   4. Filter: correlation >= threshold, then select minimum RMSE
#   5. Refine best LHS parameters with L-BFGS-B optimisation
#
# Designed for per-LME calibration with LME-averaged inputs,
# with calibrated values then transferable to grid-cell runs.
# =============================================================================


#' Default Parameter Bounds for Fishing Calibration
#'
#' @title Get default parameter bounds for fishing calibration
#' @description Returns a list of lower and upper bounds for catchability (q)
#'   and minimum selectivity weight (w_min, in log10 grams) for each fish group.
#'   Bounds are based on the ZooMSS fishing calibration framework (see conceptual figure).
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{q_lower}: Named numeric vector of lower bounds for q (length 3)
#'     \item \code{q_upper}: Named numeric vector of upper bounds for q (length 3)
#'     \item \code{wmin_lower}: Named numeric vector of lower bounds for log10(w_min) (length 3)
#'     \item \code{wmin_upper}: Named numeric vector of upper bounds for log10(w_min) (length 3)
#'   }
#' @export
#'
#' @examples
#' bounds <- getDefaultBounds()
#' bounds$q_lower
#' bounds$q_upper
#'
getDefaultBounds <- function() {
  list(
    # Catchability bounds
    q_lower = c(Fish_Small = 0.001, Fish_Med = 0.005, Fish_Large = 0.01),
    q_upper = c(Fish_Small = 0.05,  Fish_Med = 0.08,  Fish_Large = 0.10),
    # Minimum selectivity weight bounds (log10 grams)
    wmin_lower = c(Fish_Small = 0.0,  Fish_Med = 1.0, Fish_Large = 2.3),
    wmin_upper = c(Fish_Small = 1.0,  Fish_Med = 2.0, Fish_Large = 2.7)
  )
}


#' Fishing Calibration Objective Function
#'
#' @title Evaluate predicted vs observed catch for a given parameter set
#' @description Runs ZooMSS with the specified catchability and selectivity parameters,
#'   extracts predicted catch over the calibration period, and computes the sum of
#'   squared errors (SSE) against observed catch. This is the function minimised
#'   during calibration.
#'
#' @param par Numeric vector of length 6: c(q_small, q_med, q_large,
#'   wmin_small, wmin_med, wmin_large). The w_min values are in log10 grams.
#' @param input_params Data frame of environmental + effort time series
#'   (created by \code{createInputParams()} + \code{addFishingEffort()}).
#' @param Groups Base Groups data frame (from \code{getGroups()}).
#' @param observed_catch Data frame with columns: \code{time}, \code{catch_small},
#'   \code{catch_med}, \code{catch_large}. Observed catch in g wet weight,
#'   at the same temporal resolution as model output.
#' @param calibration_period Numeric vector of length 2: c(start_year, end_year).
#'   Only time steps within this range are used for SSE calculation.
#' @param isave Integer. Save frequency for model output (default: 2).
#' @param metric Character. Error metric to return: "SSE" (default), "RMSE",
#'   or "negcor" (negative correlation, for maximising correlation via minimisation).
#' @param verbose Logical. If TRUE, print progress during evaluation (default: FALSE).
#'
#' @return Numeric scalar: the error metric value. Returns \code{Inf} if the
#'   model fails or produces invalid output.
#' @export
#'
#' @examples
#' \dontrun{
#' par <- c(0.01, 0.02, 0.03, 0.5, 1.5, 2.5)
#' sse <- fishingObjective(par, input_params, Groups, observed_catch,
#'                         calibration_period = c(1960, 2010))
#' }
#'
fishingObjective <- function(par, input_params, Groups, observed_catch,
                             calibration_period = c(1960, 2010),
                             isave = 2, metric = "SSE", verbose = FALSE) {

  # Unpack parameters
  q_small  <- par[1]
  q_med    <- par[2]
  q_large  <- par[3]
  wmin_small <- par[4]
  wmin_med   <- par[5]
  wmin_large <- par[6]

  if (verbose) {
    cat(sprintf("  q = [%.4f, %.4f, %.4f], w_min = [%.2f, %.2f, %.2f]\n",
                q_small, q_med, q_large, wmin_small, wmin_med, wmin_large))
  }

  # Update Groups with current parameter values
  Groups_cal <- tryCatch({
    setFishingParams(Groups,
                     q_small = q_small, q_med = q_med, q_large = q_large,
                     w_min_small = wmin_small, w_min_med = wmin_med,
                     w_min_large = wmin_large)
  }, error = function(e) {
    if (verbose) cat("  setFishingParams failed:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(Groups_cal)) return(Inf)

  # Run the model (suppress output during calibration)
  mdl <- tryCatch({
    suppressMessages(
      zoomss_model(input_params = input_params, Groups = Groups_cal, isave = isave)
    )
  }, error = function(e) {
    if (verbose) cat("  Model run failed:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(mdl) || is.null(mdl$catch)) return(Inf)

  # Extract predicted catch within calibration period
  time_idx <- which(mdl$time >= calibration_period[1] & mdl$time <= calibration_period[2])
  if (length(time_idx) == 0) {
    if (verbose) cat("  No time steps in calibration period\n")
    return(Inf)
  }

  pred_catch <- mdl$catch[time_idx, , drop = FALSE]  # nsave_cal x 3
  pred_time  <- mdl$time[time_idx]

  # Align observed catch to model time steps (nearest-neighbour matching)
  obs_aligned <- alignObservedCatch(observed_catch, pred_time)
  if (is.null(obs_aligned)) {
    if (verbose) cat("  Failed to align observed catch\n")
    return(Inf)
  }

  # Check for valid values
  if (any(!is.finite(pred_catch)) || any(!is.finite(obs_aligned))) {
    if (verbose) cat("  Non-finite values in catch\n")
    return(Inf)
  }

  # Calculate error metric
  residuals <- pred_catch - obs_aligned

  result <- switch(metric,
    "SSE"    = sum(residuals^2),
    "RMSE"   = sqrt(mean(residuals^2)),
    "negcor" = {
      # Negative mean correlation across fish groups (for minimisation)
      cors <- sapply(1:ncol(pred_catch), function(f) {
        if (sd(pred_catch[, f]) == 0 || sd(obs_aligned[, f]) == 0) return(0)
        cor(pred_catch[, f], obs_aligned[, f])
      })
      -mean(cors, na.rm = TRUE)
    },
    stop("Unknown metric: ", metric, ". Use 'SSE', 'RMSE', or 'negcor'.")
  )

  if (verbose) cat(sprintf("  %s = %.6g\n", metric, result))
  return(result)
}


#' Align Observed Catch to Model Time Steps
#'
#' @title Match observed catch data to model output time points
#' @description Uses nearest-neighbour matching to align observed catch time series
#'   to the model's output time grid. This handles differences in temporal resolution
#'   between observed data (typically annual) and model output (sub-annual).
#'
#' @param observed_catch Data frame with columns: \code{time}, \code{catch_small},
#'   \code{catch_med}, \code{catch_large}.
#' @param model_time Numeric vector of model output time points.
#'
#' @return Matrix (length(model_time) x 3) of aligned observed catch values,
#'   or NULL if alignment fails.
#' @noRd
#'
alignObservedCatch <- function(observed_catch, model_time) {

  required_cols <- c("time", "catch_small", "catch_med", "catch_large")
  if (!all(required_cols %in% names(observed_catch))) {
    warning("observed_catch must have columns: ", paste(required_cols, collapse = ", "))
    return(NULL)
  }

  obs_time <- observed_catch$time
  obs_mat  <- as.matrix(observed_catch[, c("catch_small", "catch_med", "catch_large")])

  # Nearest-neighbour matching: for each model time, find closest observed time
  aligned <- matrix(NA, nrow = length(model_time), ncol = 3)
  for (i in seq_along(model_time)) {
    nearest_idx <- which.min(abs(obs_time - model_time[i]))
    aligned[i, ] <- obs_mat[nearest_idx, ]
  }

  return(aligned)
}


#' Latin Hypercube Search for Fishing Parameters
#'
#' @title Explore parameter space using Latin Hypercube Sampling
#' @description Generates parameter combinations using LHS, evaluates each by
#'   running ZooMSS and computing correlation and RMSE against observed catch.
#'   This is the first step of the DBPM-style calibration workflow: broad
#'   exploration before targeted optimisation.
#'
#' @param input_params Data frame of environmental + effort time series.
#' @param Groups Base Groups data frame.
#' @param observed_catch Data frame with columns: time, catch_small, catch_med, catch_large.
#' @param calibration_period Numeric vector c(start_year, end_year).
#' @param n_samples Integer. Number of LHS samples to evaluate (default: 500).
#' @param bounds List of parameter bounds (from \code{getDefaultBounds()}).
#' @param isave Integer. Model save frequency (default: 2).
#' @param cores Integer. Number of cores for parallel evaluation (default: 1, sequential).
#' @param verbose Logical. Print progress (default: TRUE).
#'
#' @return Data frame with columns: q_small, q_med, q_large, wmin_small, wmin_med,
#'   wmin_large, SSE, RMSE, cor_small, cor_med, cor_large, cor_mean. Each row
#'   corresponds to one LHS sample.
#' @export
#'
#' @examples
#' \dontrun{
#' lhs_results <- lhsSearch(input_params, Groups, observed_catch,
#'                           calibration_period = c(1960, 2010),
#'                           n_samples = 200)
#'
#' # Filter: correlation >= 0.5, then select minimum RMSE
#' good <- lhs_results[lhs_results$cor_mean >= 0.5, ]
#' best <- good[which.min(good$RMSE), ]
#' }
#'
lhsSearch <- function(input_params, Groups, observed_catch,
                      calibration_period = c(1960, 2010),
                      n_samples = 500, bounds = NULL,
                      isave = 2, cores = 1, verbose = TRUE) {

  if (is.null(bounds)) bounds <- getDefaultBounds()

  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("Package 'lhs' is required for Latin Hypercube Sampling.\n",
         "Install with: install.packages('lhs')")
  }

  # Generate LHS design (6 parameters, n_samples points)
  lhs_design <- lhs::randomLHS(n_samples, 6)

  # Scale to parameter bounds
  # Columns: q_small, q_med, q_large, wmin_small, wmin_med, wmin_large
  lower <- c(bounds$q_lower, bounds$wmin_lower)
  upper <- c(bounds$q_upper, bounds$wmin_upper)

  par_matrix <- sweep(sweep(lhs_design, 2, upper - lower, "*"), 2, lower, "+")
  colnames(par_matrix) <- c("q_small", "q_med", "q_large",
                             "wmin_small", "wmin_med", "wmin_large")

  if (verbose) {
    cat("LHS Search: evaluating", n_samples, "parameter combinations\n")
    cat("Parameter bounds:\n")
    for (i in 1:6) {
      cat(sprintf("  %s: [%.4f, %.4f]\n", colnames(par_matrix)[i], lower[i], upper[i]))
    }
  }

  # Evaluation function for a single parameter set
  eval_one <- function(i) {
    par_i <- par_matrix[i, ]

    # Run model and get predicted catch
    Groups_i <- tryCatch(
      suppressMessages(
        setFishingParams(Groups,
                         q_small = par_i[1], q_med = par_i[2], q_large = par_i[3],
                         w_min_small = par_i[4], w_min_med = par_i[5],
                         w_min_large = par_i[6])
      ),
      error = function(e) NULL
    )
    if (is.null(Groups_i)) return(rep(NA, 6))  # SSE, RMSE, cor_s, cor_m, cor_l, cor_mean

    mdl <- tryCatch(
      suppressMessages(
        zoomss_model(input_params = input_params, Groups = Groups_i, isave = isave)
      ),
      error = function(e) NULL
    )
    if (is.null(mdl) || is.null(mdl$catch)) return(rep(NA, 6))

    # Extract calibration period
    time_idx <- which(mdl$time >= calibration_period[1] & mdl$time <= calibration_period[2])
    if (length(time_idx) < 2) return(rep(NA, 6))

    pred_catch <- mdl$catch[time_idx, , drop = FALSE]
    pred_time  <- mdl$time[time_idx]

    obs_aligned <- alignObservedCatch(observed_catch, pred_time)
    if (is.null(obs_aligned)) return(rep(NA, 6))

    # Check validity
    if (any(!is.finite(pred_catch)) || any(!is.finite(obs_aligned))) return(rep(NA, 6))

    # Calculate metrics
    residuals <- pred_catch - obs_aligned
    sse  <- sum(residuals^2)
    rmse <- sqrt(mean(residuals^2))

    # Per-group correlations
    cors <- sapply(1:3, function(f) {
      if (sd(pred_catch[, f]) == 0 || sd(obs_aligned[, f]) == 0) return(NA)
      cor(pred_catch[, f], obs_aligned[, f])
    })

    c(sse, rmse, cors, mean(cors, na.rm = TRUE))
  }

  # Run evaluations (parallel or sequential)
  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    if (verbose) cat("Running in parallel on", cores, "cores\n")
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export required objects to workers
    parallel::clusterExport(cl, c("input_params", "Groups", "observed_catch",
                                   "calibration_period", "isave", "par_matrix"),
                            envir = environment())
    parallel::clusterEvalQ(cl, library(zoomss))

    results_list <- parallel::parLapply(cl, 1:n_samples, eval_one)
  } else {
    if (verbose) pb <- progress::progress_bar$new(
      format = "LHS [:bar] :current/:total (:percent) eta: :eta",
      total = n_samples, width = 60, show_after = 0
    )
    results_list <- lapply(1:n_samples, function(i) {
      if (verbose) pb$tick()
      eval_one(i)
    })
  }

  # Assemble results
  metrics_mat <- do.call(rbind, results_list)
  colnames(metrics_mat) <- c("SSE", "RMSE", "cor_small", "cor_med", "cor_large", "cor_mean")

  results_df <- as.data.frame(cbind(par_matrix, metrics_mat))

  # Remove failed runs
  n_failed <- sum(is.na(results_df$SSE))
  if (verbose && n_failed > 0) {
    cat(sprintf("  %d of %d runs failed (%.0f%%)\n",
                n_failed, n_samples, 100 * n_failed / n_samples))
  }

  if (verbose) {
    valid <- results_df[!is.na(results_df$SSE), ]
    if (nrow(valid) > 0) {
      cat(sprintf("  Best RMSE: %.4g | Best mean correlation: %.3f\n",
                  min(valid$RMSE), max(valid$cor_mean, na.rm = TRUE)))
    }
  }

  return(results_df)
}


#' Calibrate Fishing Parameters
#'
#' @title Full calibration workflow: LHS exploration + L-BFGS-B refinement
#' @description Runs the complete DBPM-style calibration workflow:
#'   1. Latin Hypercube Sampling to explore parameter space
#'   2. Filter by correlation threshold and minimum RMSE
#'   3. Refine best parameters using L-BFGS-B bounded optimisation
#'
#'   Designed for per-LME calibration. Calibrated parameters can then be used
#'   as initial values for grid-cell-level runs.
#'
#' @param input_params Data frame of environmental + effort time series.
#' @param Groups Base Groups data frame (from \code{getGroups()}).
#' @param observed_catch Data frame with columns: time, catch_small, catch_med, catch_large.
#' @param calibration_period Numeric vector c(start_year, end_year) (default: c(1960, 2010)).
#' @param n_lhs Integer. Number of LHS samples (default: 500).
#' @param cor_threshold Numeric. Minimum mean correlation to pass LHS filter (default: 0.5).
#' @param bounds List of parameter bounds (default: from \code{getDefaultBounds()}).
#' @param refine Logical. If TRUE (default), refine best LHS result with L-BFGS-B.
#' @param isave Integer. Model save frequency (default: 2).
#' @param cores Integer. Number of cores for parallel LHS evaluation (default: 1).
#' @param verbose Logical. Print progress (default: TRUE).
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{par}: Named numeric vector of calibrated parameters
#'       (q_small, q_med, q_large, wmin_small, wmin_med, wmin_large)
#'     \item \code{q}: Named numeric vector of calibrated catchability values
#'     \item \code{wmin}: Named numeric vector of calibrated selectivity thresholds (log10 g)
#'     \item \code{SSE}: Final SSE value
#'     \item \code{RMSE}: Final RMSE value
#'     \item \code{correlations}: Per-group correlations from the final parameters
#'     \item \code{lhs_results}: Full LHS search results data frame
#'     \item \code{lhs_best}: Best LHS parameters (before refinement)
#'     \item \code{optim_result}: Output from \code{optim()} (if refine = TRUE)
#'     \item \code{convergence}: 0 = converged (from optim), -1 = LHS only
#'     \item \code{bounds}: Parameter bounds used
#'     \item \code{calibration_period}: Time period used for calibration
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare observed catch data
#' obs_catch <- data.frame(
#'   time = 1960:2010,
#'   catch_small = runif(51, 100, 500),
#'   catch_med   = runif(51, 200, 800),
#'   catch_large = runif(51, 50, 300)
#' )
#'
#' # Run calibration
#' cal <- calibrateFishing(input_params, Groups, obs_catch,
#'                          n_lhs = 200, cores = 4)
#'
#' # Apply calibrated parameters
#' Groups_cal <- setFishingParams(Groups,
#'                                q_small = cal$q["Fish_Small"],
#'                                q_med   = cal$q["Fish_Med"],
#'                                q_large = cal$q["Fish_Large"],
#'                                w_min_small = cal$wmin["Fish_Small"],
#'                                w_min_med   = cal$wmin["Fish_Med"],
#'                                w_min_large = cal$wmin["Fish_Large"])
#' }
#'
calibrateFishing <- function(input_params, Groups, observed_catch,
                             calibration_period = c(1960, 2010),
                             n_lhs = 500, cor_threshold = 0.5,
                             bounds = NULL, refine = TRUE,
                             isave = 2, cores = 1, verbose = TRUE) {

  if (is.null(bounds)) bounds <- getDefaultBounds()

  # =========================================================================
  # STEP 1: Latin Hypercube Sampling
  # =========================================================================
  if (verbose) cat("\n=== Step 1: LHS Exploration (", n_lhs, "samples) ===\n")

  lhs_results <- lhsSearch(
    input_params = input_params,
    Groups = Groups,
    observed_catch = observed_catch,
    calibration_period = calibration_period,
    n_samples = n_lhs,
    bounds = bounds,
    isave = isave,
    cores = cores,
    verbose = verbose
  )

  # =========================================================================
  # STEP 2: Filter by correlation threshold, select minimum RMSE
  # =========================================================================
  if (verbose) cat("\n=== Step 2: Filter & Select Best LHS Parameters ===\n")

  valid <- lhs_results[!is.na(lhs_results$SSE), ]

  if (nrow(valid) == 0) {
    stop("All LHS samples failed. Check input data and parameter bounds.")
  }

  # Apply correlation filter
  passing <- valid[!is.na(valid$cor_mean) & valid$cor_mean >= cor_threshold, ]

  if (nrow(passing) == 0) {
    if (verbose) {
      cat(sprintf("  No samples passed cor_threshold = %.2f\n", cor_threshold))
      cat(sprintf("  Best mean correlation achieved: %.3f\n", max(valid$cor_mean, na.rm = TRUE)))
      cat("  Falling back to minimum RMSE without correlation filter\n")
    }
    passing <- valid
  } else {
    if (verbose) {
      cat(sprintf("  %d of %d valid samples passed correlation threshold (>= %.2f)\n",
                  nrow(passing), nrow(valid), cor_threshold))
    }
  }

  # Select best by minimum RMSE
  best_idx <- which.min(passing$RMSE)
  lhs_best <- passing[best_idx, ]

  par_lhs <- as.numeric(lhs_best[1:6])
  names(par_lhs) <- c("q_small", "q_med", "q_large", "wmin_small", "wmin_med", "wmin_large")

  if (verbose) {
    cat("  Best LHS parameters:\n")
    cat(sprintf("    q:     [%.4f, %.4f, %.4f]\n", par_lhs[1], par_lhs[2], par_lhs[3]))
    cat(sprintf("    w_min: [%.2f, %.2f, %.2f] (log10 g)\n", par_lhs[4], par_lhs[5], par_lhs[6]))
    cat(sprintf("    RMSE: %.4g | Mean cor: %.3f\n", lhs_best$RMSE, lhs_best$cor_mean))
  }

  # =========================================================================
  # STEP 3: L-BFGS-B Refinement
  # =========================================================================
  optim_result <- NULL
  par_final <- par_lhs

  if (refine) {
    if (verbose) cat("\n=== Step 3: L-BFGS-B Refinement ===\n")

    lower <- c(bounds$q_lower, bounds$wmin_lower)
    upper <- c(bounds$q_upper, bounds$wmin_upper)

    optim_result <- tryCatch({
      optim(par = par_lhs,
            fn = fishingObjective,
            method = "L-BFGS-B",
            lower = lower,
            upper = upper,
            input_params = input_params,
            Groups = Groups,
            observed_catch = observed_catch,
            calibration_period = calibration_period,
            isave = isave,
            metric = "SSE",
            verbose = FALSE,
            control = list(maxit = 100, factr = 1e7))
    }, error = function(e) {
      if (verbose) cat("  L-BFGS-B failed:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(optim_result)) {
      par_final <- optim_result$par
      names(par_final) <- names(par_lhs)

      if (verbose) {
        cat("  Optimisation", ifelse(optim_result$convergence == 0, "converged", "did not converge"), "\n")
        cat("  Refined parameters:\n")
        cat(sprintf("    q:     [%.4f, %.4f, %.4f]\n", par_final[1], par_final[2], par_final[3]))
        cat(sprintf("    w_min: [%.2f, %.2f, %.2f] (log10 g)\n", par_final[4], par_final[5], par_final[6]))
        cat(sprintf("    SSE: %.4g (was %.4g from LHS)\n", optim_result$value, lhs_best$SSE))
      }
    } else {
      if (verbose) cat("  Refinement failed — using LHS best parameters\n")
    }
  }

  # =========================================================================
  # STEP 4: Compute final diagnostics
  # =========================================================================
  if (verbose) cat("\n=== Final Evaluation ===\n")

  # Run final model to get correlations
  final_eval <- evaluateCalibration(par_final, input_params, Groups,
                                     observed_catch, calibration_period, isave)

  # Package results
  q_final <- par_final[1:3]
  names(q_final) <- c("Fish_Small", "Fish_Med", "Fish_Large")
  wmin_final <- par_final[4:6]
  names(wmin_final) <- c("Fish_Small", "Fish_Med", "Fish_Large")

  if (verbose) {
    cat("  Final RMSE:", round(final_eval$RMSE, 4), "\n")
    cat("  Correlations:", paste(names(final_eval$correlations), "=",
                                  round(final_eval$correlations, 3), collapse = ", "), "\n")
    cat("\nCalibration complete.\n")
  }

  list(
    par = par_final,
    q = q_final,
    wmin = wmin_final,
    SSE = final_eval$SSE,
    RMSE = final_eval$RMSE,
    correlations = final_eval$correlations,
    lhs_results = lhs_results,
    lhs_best = lhs_best,
    optim_result = optim_result,
    convergence = if (!is.null(optim_result)) optim_result$convergence else -1L,
    bounds = bounds,
    calibration_period = calibration_period
  )
}


#' Evaluate Calibration Quality
#'
#' @title Compute diagnostic metrics for a calibrated parameter set
#' @description Runs ZooMSS with the given parameters and computes SSE, RMSE,
#'   and per-group correlations against observed catch.
#'
#' @param par Numeric vector of length 6 (q_small, q_med, q_large, wmin_small, wmin_med, wmin_large).
#' @param input_params Data frame of environmental + effort time series.
#' @param Groups Base Groups data frame.
#' @param observed_catch Data frame with time, catch_small, catch_med, catch_large columns.
#' @param calibration_period Numeric vector c(start_year, end_year).
#' @param isave Integer. Model save frequency.
#'
#' @return A list with: SSE, RMSE, correlations (named vector), predicted_catch, observed_catch_aligned, time.
#' @export
#'
evaluateCalibration <- function(par, input_params, Groups, observed_catch,
                                calibration_period = c(1960, 2010), isave = 2) {

  Groups_cal <- suppressMessages(
    setFishingParams(Groups,
                     q_small = par[1], q_med = par[2], q_large = par[3],
                     w_min_small = par[4], w_min_med = par[5], w_min_large = par[6])
  )

  mdl <- suppressMessages(
    zoomss_model(input_params = input_params, Groups = Groups_cal, isave = isave)
  )

  time_idx <- which(mdl$time >= calibration_period[1] & mdl$time <= calibration_period[2])
  pred_catch <- mdl$catch[time_idx, , drop = FALSE]
  pred_time  <- mdl$time[time_idx]
  obs_aligned <- alignObservedCatch(observed_catch, pred_time)

  residuals <- pred_catch - obs_aligned
  sse  <- sum(residuals^2)
  rmse <- sqrt(mean(residuals^2))

  fish_names <- c("Fish_Small", "Fish_Med", "Fish_Large")
  cors <- sapply(1:3, function(f) {
    if (sd(pred_catch[, f]) == 0 || sd(obs_aligned[, f]) == 0) return(NA)
    cor(pred_catch[, f], obs_aligned[, f])
  })
  names(cors) <- fish_names

  list(
    SSE = sse,
    RMSE = rmse,
    correlations = cors,
    predicted_catch = pred_catch,
    observed_catch_aligned = obs_aligned,
    time = pred_time
  )
}


#' Plot Calibration Results
#'
#' @title Visualize predicted vs observed catch from calibration
#' @description Creates a multi-panel plot comparing predicted and observed catch
#'   time series for each fish group, with correlation and RMSE annotations.
#'
#' @param cal_eval Output from \code{evaluateCalibration()}.
#' @param colours Optional named vector of colours for fish groups.
#'
#' @return A ggplot2 object (requires ggplot2 and patchwork).
#' @export
#'
plotCalibration <- function(cal_eval, colours = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotCalibration.")
  }

  fish_names <- c("Fish_Small", "Fish_Med", "Fish_Large")
  n_time <- length(cal_eval$time)

  # Build data frame
  df <- data.frame(
    time = rep(cal_eval$time, times = 6),
    catch = c(as.vector(cal_eval$predicted_catch),
              as.vector(cal_eval$observed_catch_aligned)),
    group = rep(rep(fish_names, each = n_time), 2),
    source = rep(c("Predicted", "Observed"), each = n_time * 3)
  )

  if (is.null(colours)) {
    colours <- c("Fish_Small" = "#E69F00", "Fish_Med" = "#56B4E9", "Fish_Large" = "#009E73")
  }

  # Annotation labels
  cor_labels <- paste0("r = ", round(cal_eval$correlations, 3))
  names(cor_labels) <- fish_names

  plots <- lapply(fish_names, function(fn) {
    df_sub <- df[df$group == fn, ]
    ggplot2::ggplot(df_sub, ggplot2::aes(x = .data$time, y = .data$catch,
                                          linetype = .data$source)) +
      ggplot2::geom_line(colour = colours[fn], linewidth = 0.8) +
      ggplot2::scale_linetype_manual(values = c("Predicted" = "solid", "Observed" = "dashed")) +
      ggplot2::theme_bw() +
      ggplot2::annotate("text", x = min(df_sub$time) + 2, y = max(df_sub$catch) * 0.95,
                         label = cor_labels[fn], hjust = 0, size = 3.5) +
      ggplot2::labs(x = "Time (years)", y = "Catch (g WW)",
                     title = fn, linetype = "Source")
  })

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plots, ncol = 1, guides = "collect") +
      patchwork::plot_annotation(
        title = "Fishing Calibration: Predicted vs Observed Catch",
        subtitle = paste0("Overall RMSE: ", round(cal_eval$RMSE, 4))
      )
  } else {
    plots[[1]]
  }
}


#' Plot LHS Search Results
#'
#' @title Visualize the LHS parameter search landscape
#' @description Creates diagnostic plots of the LHS search showing how error
#'   metrics vary across the parameter space.
#'
#' @param lhs_results Data frame from \code{lhsSearch()} or \code{calibrateFishing()$lhs_results}.
#' @param cor_threshold Numeric. Correlation threshold to highlight (default: 0.5).
#'
#' @return A ggplot2 object.
#' @export
#'
plotLHSResults <- function(lhs_results, cor_threshold = 0.5) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotLHSResults.")
  }

  valid <- lhs_results[!is.na(lhs_results$SSE), ]
  valid$passes_filter <- !is.na(valid$cor_mean) & valid$cor_mean >= cor_threshold

  # q vs RMSE for each fish group
  par_names <- c("q_small", "q_med", "q_large")
  fish_labels <- c("Fish_Small", "Fish_Med", "Fish_Large")

  plots <- lapply(seq_along(par_names), function(i) {
    ggplot2::ggplot(valid, ggplot2::aes(x = .data[[par_names[i]]], y = .data$RMSE,
                                         colour = .data$cor_mean)) +
      ggplot2::geom_point(alpha = 0.5, size = 1.5) +
      ggplot2::scale_colour_viridis_c(name = "Mean\ncorrelation", limits = c(-1, 1)) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = paste0("q (", fish_labels[i], ")"), y = "RMSE",
                     title = fish_labels[i])
  })

  # Correlation histogram
  p_hist <- ggplot2::ggplot(valid, ggplot2::aes(x = .data$cor_mean)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = cor_threshold, linetype = "dashed", colour = "red") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Mean Correlation", y = "Count",
                   title = paste0("LHS Correlation Distribution (threshold = ", cor_threshold, ")"))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(c(plots, list(p_hist)), ncol = 2) +
      patchwork::plot_annotation(title = "LHS Parameter Search Results")
  } else {
    plots[[1]]
  }
}
