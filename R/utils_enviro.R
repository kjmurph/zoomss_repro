#' Create ZooMSS Input Parameters Object
#'
#' @title Create input parameters data frame for ZooMSS model runs
#' @description Creates a properly formatted input parameters data frame for ZooMSS model
#'   simulations, combining temporal parameters with environmental time series data.
#' @details This function combines environmental time series (SST and chlorophyll) with
#'   time data to create the input_params object required by zoomss_model().
#'   The function performs validation checks using assertthat to ensure:
#'   - All input vectors are numeric and of equal length
#'   - SST values are within reasonable ocean range (-2 to 35 deg C)
#'   - Chlorophyll values are positive and within typical range (0 to 50 mg/m^3)
#'   - Time values are increasing and reasonable
#'
#'
#' @param time Numeric vector of time values in years (must be increasing and uniform, can start at any value)
#' @param sst Numeric vector of sea surface temperature values in deg C
#' @param chl Numeric vector of chlorophyll concentration values in mg/m^3
#' @param cellID Optional numeric vector of cell identifiers for spatial data (default: NULL)
#'
#' @return Data frame with columns: time, time_step, sst, chl, and cellID (if provided)
#' @export
#'
#' @examples
#' \dontrun{
#' # Create simple environmental time series
#' time_vec <- seq(0, 10, 0.01)  # 10 years with 0.01 year time steps
#' sst_vec <- 15 + 3*sin(2*pi*time_vec/1)  # annual cycle
#' chl_vec <- 0.5 + 0.2*cos(2*pi*time_vec/1)  # annual cycle
#'
#' # Create input parameters object
#' input_params <- createInputParams(time_vec, sst_vec, chl_vec)
#'
#' # Use with ZooMSS model
#' results <- zoomss_model(input_params, Groups, isave = 50)
#' }
#'
createInputParams <- function(time,
                              sst,
                              chl,
                              cellID = NULL) {

  # Load assertthat package for validation
  if (!requireNamespace("assertthat", quietly = TRUE)) {
    stop("assertthat package required for input validation")
  }

  # Validate input data types and structure
  assertthat::assert_that(is.numeric(time), msg = "time must be numeric")
  assertthat::assert_that(is.numeric(sst), msg = "sst must be numeric")
  assertthat::assert_that(is.numeric(chl), msg = "chl must be numeric")

  # Validate equal lengths if length of sst and chl > 1
  if (length(sst) > 1 && length(chl) > 1){
    assertthat::assert_that(length(time) == length(sst),
                            msg = "time and sst must have the same length")
    assertthat::assert_that(length(time) == length(chl),
                            msg = "time and chl must have the same length")
  }

  # Validate cellID if provided
  if (!is.null(cellID)) {
    assertthat::assert_that(is.numeric(cellID), msg = "cellID must be numeric")
    assertthat::assert_that(length(cellID) == length(time),
                            msg = "cellID must have the same length as time")
  }

  # Validate time vector properties
  assertthat::assert_that(length(time) > 1, msg = "time must have at least 2 values")
  assertthat::assert_that(all(!is.na(time)), msg = "time cannot contain NA values")
  assertthat::assert_that(all(diff(time) > 0), msg = "time values must be increasing")

  # Calculate dt and tmax from time vector
  dt_values <- diff(time)
  dt <- dt_values[1]  # Use first time step as dt

  # Check if time steps are uniform - ERROR if not consistent
  max_dt_diff <- max(abs(dt_values - dt))
  if (max_dt_diff > dt * 0.001) {  # Allow only 0.1% variation (much stricter)
    stop("Time steps are not uniform. Maximum deviation: ", round(max_dt_diff, 6),
         " (", round(100 * max_dt_diff / dt, 2), "% of dt). ",
         "ZooMSS requires uniform time steps for accurate results.")
  }

  tmax <- max(time)  # Maximum time value (not duration)

  # Validate temporal parameters
  assertthat::assert_that(dt > 0, msg = "calculated dt must be positive")
  # Note: tmax can be any value (positive, negative, or zero) as it's the final time point

  # Validate environmental data ranges
  assertthat::assert_that(all(!is.na(sst)), msg = "sst cannot contain NA values")
  assertthat::assert_that(all(!is.na(chl)), msg = "chl cannot contain NA values")
  assertthat::assert_that(all(sst >= -2 & sst <= 35),
                          msg = "sst values must be within ocean range (-2 to 35 deg C)")
  assertthat::assert_that(all(chl >= 0 & chl <= 50),
                          msg = "chl values must be within range (0 to 50 mg/m^3)")

  # Create formatted data frame
  if (is.null(cellID)) {
    formatted_data <- data.frame(
      time = time,
      time_step = seq_along(time),
      sst = sst,
      chl = chl
    )
  } else {
    formatted_data <- data.frame(
      time = time,
      time_step = seq_along(time),
      sst = sst,
      chl = chl,
      cellID = cellID
    )
  }

  # Provide summary information
  n_time_points <- nrow(formatted_data)
  n_time_steps <- n_time_points - 1

  cat("ZooMSS input parameters created:\n")
  cat("- Time points:", n_time_points, "(time values provided)\n")
  cat("- Time steps:", n_time_steps, "(intervals to simulate)\n")
  cat("- Time range:", round(min(formatted_data$time), 3), "to",
      round(max(formatted_data$time), 3), "years\n")
  cat("- dt =", round(dt, 4), "years\n")
  cat("- SST range:", round(min(formatted_data$sst), 1), "to",
      round(max(formatted_data$sst), 1), "deg C\n")
  cat("- Chlorophyll range:", round(min(formatted_data$chl), 2), "to",
      round(max(formatted_data$chl), 2), "mg/m^3\n")

  # Helpful reminder about time vector interpretation
  if (length(time) > 1 && all(diff(time) == 1) && min(time) %% 1 == 0 && max(time) %% 1 == 0) {
    cat("- Note: Time vector", min(time), ":", max(time), "creates", n_time_steps,
        "time steps (intervals) from", length(time), "time points.\n")
  }

  return(formatted_data)
}


#' Create Environmental Time Series
#'
#' @title Generate synthetic or custom environmental time series for ZooMSS
#' @description
#' Creates environmental time series (SST and chlorophyll) for ZooMSS testing using either:
#'  (a) internally generated synthetic series (static or seasonal), or
#'  (b) user-provided custom series, used raw or with seasonal and/or stochastic effects applied on top.
#'
#' @details
#' SST stochasticity is additive Gaussian (optional AR(1)), scaled so sst_noise_sd is the **marginal** SD.
#' Chl stochasticity is multiplicative lognormal (optional AR(1)) to keep values positive and mean-preserving.
#' A strictly positive chlorophyll floor is enforced; by default it equals base_chl (chl_floor_mode = "base").
#'
#' @param n_years,dt,time,base_sst,base_chl,seasonal,sst_amplitude,chl_amplitude,stochastic,
#'   sst_noise_sd,chl_noise_cv,ar1_phi,seed,enforce_bounds,use_custom,custom_sst,custom_chl,custom_use_raw,
#'   chl_floor_mode,chl_floor_value See discussion above.
#' @return data.frame(time, sst, chl)
#' @export
createEnviroData <- function(n_years,
                             dt,
                             time = NULL,
                             base_sst = 15,
                             base_chl = 0.5,
                             seasonal = TRUE,
                             sst_amplitude = 3,
                             chl_amplitude = 0.2,
                             stochastic = FALSE,
                             sst_noise_sd = 0,
                             chl_noise_cv = 0,
                             ar1_phi = 0,
                             seed = NULL,
                             enforce_bounds = FALSE,
                             use_custom = FALSE,
                             custom_sst = NULL,
                             custom_chl = NULL,
                             custom_use_raw = TRUE,
                             chl_floor_mode = c("base","epsilon","fixed"),
                             chl_floor_value = 1e-8) {

  # --- Basic checks ---
  if (!is.logical(seasonal) || length(seasonal) != 1L) stop("seasonal must be a single logical.")
  if (!is.logical(stochastic) || length(stochastic) != 1L) stop("stochastic must be a single logical.")
  if (!is.logical(use_custom) || length(use_custom) != 1L) stop("use_custom must be a single logical.")
  if (!is.logical(custom_use_raw) || length(custom_use_raw) != 1L) stop("custom_use_raw must be a single logical.")
  if (!is.logical(enforce_bounds) || length(enforce_bounds) != 1L) stop("enforce_bounds must be a single logical.")
  if (!is.numeric(ar1_phi) || length(ar1_phi) != 1L || ar1_phi <= -0.99 || ar1_phi >= 0.99) {
    stop("ar1_phi must be a single numeric value in (-0.99, 0.99).")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) stop("seed must be a single integer.")
    set.seed(as.integer(seed))
  }

  # --- Time handling ---
  if (!is.null(time)) {
    if (!is.numeric(time) || anyNA(time) || any(diff(time) <= 0)) {
      stop("time must be numeric, non-NA, and strictly increasing.")
    }
    time_years <- time
  } else {
    time_years <- seq(0, n_years, by = dt)
  }
  n <- length(time_years)

  # --- Chlorophyll floor ---
  chl_floor_mode <- match.arg(chl_floor_mode)
  chl_floor <- switch(
    chl_floor_mode,
    base    = max(base_chl, 1e-8),
    epsilon = max(1e-8, .Machine$double.eps),
    fixed   = {
      if (!is.numeric(chl_floor_value) || length(chl_floor_value) != 1L) {
        stop("chl_floor_value must be a single numeric when chl_floor_mode = 'fixed'.")
      }
      max(chl_floor_value, 1e-12)
    }
  )

  # --- Determine backbone and flags (avoid double seasonal) ---
  do_seasonal   <- seasonal
  do_stochastic <- stochastic

  if (isTRUE(use_custom)) {
    if (!is.null(custom_sst) && length(custom_sst) != n) stop("custom_sst length must match length(time).")
    if (!is.null(custom_chl) && length(custom_chl) != n) stop("custom_chl length must match length(time).")

    sst_values <- if (!is.null(custom_sst)) custom_sst else rep(base_sst, n)
    chl_values <- if (!is.null(custom_chl)) custom_chl else rep(base_chl, n)

    if (isTRUE(custom_use_raw)) {
      # Return custom series (plus floors/bounds), no seasonal/noise added
      do_seasonal   <- FALSE
      do_stochastic <- FALSE
    }
  } else {
    # Synthetic backbone: apply seasonal here (if requested) and DO NOT add it again later
    if (isTRUE(seasonal)) {
      sst_values <- base_sst + sst_amplitude * sin(2 * pi * time_years)
      chl_values <- base_chl + chl_amplitude * sin(2 * pi * time_years + pi)
    } else {
      sst_values <- rep(base_sst, n)
      chl_values <- rep(base_chl, n)
    }
    # Prevent double-adding seasonal in the generic step below
    do_seasonal <- FALSE
  }

  # If we are modifying custom baselines, add seasonal now
  if (isTRUE(do_seasonal)) {
    sst_values <- sst_values + sst_amplitude * sin(2 * pi * time_years)
    chl_values <- chl_values + chl_amplitude * sin(2 * pi * time_years + pi)
  }

  # Ensure strict positivity BEFORE multiplicative operations
  chl_values <- pmax(chl_values, chl_floor)

  # --- Optional stochasticity ---
  if (isTRUE(do_stochastic)) {
    # SST: additive Gaussian noise (optional AR(1)) with correct marginal SD
    if (is.numeric(sst_noise_sd) && length(sst_noise_sd) == 1L && sst_noise_sd > 0) {
      if (ar1_phi == 0) {
        eps_sst <- stats::rnorm(n, mean = 0, sd = sst_noise_sd)
      } else {
        sst_raw <- as.numeric(stats::arima.sim(model = list(ar = ar1_phi), n = n, sd = 1))
        eps_sst <- (sst_raw - mean(sst_raw)) / stats::sd(sst_raw) * sst_noise_sd
      }
      sst_values <- sst_values + eps_sst
    }

    # Chl: multiplicative lognormal noise (optional AR(1)), mean-preserving
    if (is.numeric(chl_noise_cv) && length(chl_noise_cv) == 1L && chl_noise_cv > 0) {
      sigma <- sqrt(log(1 + chl_noise_cv^2))
      mu <- -0.5 * sigma^2
      if (ar1_phi == 0) {
        Z <- stats::rnorm(n, mean = mu, sd = sigma)
      } else {
        z_raw <- as.numeric(stats::arima.sim(model = list(ar = ar1_phi), n = n, sd = 1))
        z_raw <- (z_raw - mean(z_raw)) / stats::sd(z_raw)
        Z <- mu + sigma * z_raw
      }
      mult <- exp(Z)
      chl_values <- chl_values * mult
    }
  }

  # Ensure strict positivity AFTER noise as well
  chl_values <- pmax(chl_values, chl_floor)

  # --- Optional bounds (with optional warnings) ---
  if (isTRUE(enforce_bounds)) {
    sst_before <- sst_values
    chl_before <- chl_values

    sst_values <- pmin(pmax(sst_values, -2), 35)
    chl_values <- pmin(pmax(chl_values, chl_floor), 50)

    # (Optional) Uncomment if you want warnings with clamp counts:
    # n_clamp_sst <- sum(sst_before != sst_values)
    # n_clamp_chl <- sum(chl_before != chl_values)
    # if (n_clamp_sst > 0) warning(n_clamp_sst, " SST values clamped to [-2, 35] °C")
    # if (n_clamp_chl > 0) warning(n_clamp_chl, " Chl values clamped to [", signif(chl_floor, 3), ", 50] mg/m^3")
  }

  return(data.frame(
    time = time_years,
    sst  = sst_values,
    chl  = chl_values
  ))
}
