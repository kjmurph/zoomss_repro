#' Create Environmental Time Series
#'
#' @title Generate synthetic or custom environmental time series for ZooMSS
#' @description
#' Creates environmental time series (SST and chlorophyll) for ZooMSS testing using either:
#'  (a) internally generated synthetic series (static or seasonal), or
#'  (b) user-provided custom series, used raw or with seasonal and/or stochastic
#'      effects applied on top.
#'
#' @details
#' Deterministic seasonal cycles use sinusoids with SST and Chl out of phase.
#' Optional stochasticity:
#'   * SST: additive Gaussian noise (optionally AR(1)-correlated).
#'   * Chl: multiplicative lognormal noise (optionally AR(1)-correlated), preserving positivity
#'     and approximately preserving the deterministic mean via mean-1 multipliers.
#'
#' If `use_custom = TRUE`, you may supply `custom_sst` and/or `custom_chl`. If `custom_use_raw = TRUE`,
#' the custom series are returned as-is (optionally clamped). If `custom_use_raw = FALSE`, seasonal
#' signals and/or stochasticity are applied on top of the custom series.
#'
#' @param n_years Number of years to generate (ignored if `time` is provided).
#' @param dt Time step size in years (ignored if `time` is provided).
#' @param time Optional numeric vector of time (years). If provided, length must
#'   match custom vectors if `use_custom = TRUE`. Must be strictly increasing and (ideally) uniform.
#' @param base_sst Base sea surface temperature in deg C (default: 15).
#' @param base_chl Base chlorophyll concentration in mg/m^3 (default: 0.5).
#' @param seasonal Logical, whether to add seasonal variation (default: TRUE).
#' @param sst_amplitude Amplitude of SST seasonal variations in deg C (default: 3).
#' @param chl_amplitude Amplitude of chlorophyll seasonal variations in mg/m^3 (default: 0.2).
#' @param stochastic Logical, whether to add stochastic variation (default: FALSE).
#' @param sst_noise_sd Standard deviation of additive SST noise in deg C (default: 0).
#' @param chl_noise_cv Coefficient of variation for multiplicative lognormal noise on chlorophyll (default: 0).
#' @param ar1_phi AR(1) autocorrelation parameter for noise (applied to both SST and Chl noise if used).
#'   Must be in (-0.99, 0.99). Set 0 for white noise (default: 0).
#' @param seed Optional integer seed for reproducibility (default: NULL).
#' @param enforce_bounds Logical, clamp SST to [-2, 35] and Chl to [0, 50] (default: TRUE).
#' @param use_custom Logical, if TRUE use `custom_sst` / `custom_chl` (default: FALSE).
#' @param custom_sst Optional numeric vector of custom SST values. If `use_custom=TRUE` and provided,
#'   must match `length(time)` (if `time` provided) or implied length from `n_years/dt`.
#' @param custom_chl Optional numeric vector of custom Chl values (same length rules as `custom_sst`).
#' @param custom_use_raw Logical, if TRUE use custom vectors as-is (no seasonal/noise). If FALSE,
#'   apply seasonal/noise on top of the custom series (default: TRUE).
#'
#' @return Data frame with columns: time, sst, chl
#' @export
#'
#' @examples
#' # A) Deterministic seasonal data (original behavior)
#' env1 <- createEnviroData(n_years = 2, dt = 0.01, seasonal = TRUE)
#'
#' # B) Use user-supplied custom series (raw)
#' t  <- seq(0, 1, by = 0.01)
#' cs <- 18 + sin(2*pi*t)        # arbitrary custom SST
#' cc <- 0.8 + 0.1*cos(2*pi*t)   # arbitrary custom Chl
#' env2 <- createEnviroData(time = t, use_custom = TRUE,
#'                          custom_sst = cs, custom_chl = cc,
#'                          custom_use_raw = TRUE)
#'
#' # C) Custom baseline + seasonal added on top
#' env3 <- createEnviroData(time = t, use_custom = TRUE,
#'                          custom_sst = cs, custom_chl = cc,
#'                          custom_use_raw = FALSE,
#'                          seasonal = TRUE, sst_amplitude = 0.5, chl_amplitude = 0.05)
#'
#' # D) Custom baseline + stochastic noise only
#' env4 <- createEnviroData(time = t, use_custom = TRUE,
#'                          custom_sst = cs, custom_chl = cc,
#'                          custom_use_raw = FALSE,
#'                          seasonal = FALSE, stochastic = TRUE, seed = 123,
#'                          sst_noise_sd = 0.3, chl_noise_cv = 0.25, ar1_phi = 0.6)
#'
#' # E) Custom baseline + seasonal + stochastic
#' env5 <- createEnviroData(time = t, use_custom = TRUE,
#'                          custom_sst = cs, custom_chl = cc,
#'                          custom_use_raw = FALSE,
#'                          seasonal = TRUE, stochastic = TRUE, seed = 42,
#'                          sst_amplitude = 0.5, chl_amplitude = 0.05,
#'                          sst_noise_sd = 0.2, chl_noise_cv = 0.2)
#'
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
                             enforce_bounds = TRUE,
                             use_custom = FALSE,
                             custom_sst = NULL,
                             custom_chl = NULL,
                             custom_use_raw = TRUE) {

  # Basic checks for noise options
  if (!is.logical(stochastic) || length(stochastic) != 1L) {
    stop("stochastic must be a single logical value.")
  }
  if (!is.numeric(ar1_phi) || length(ar1_phi) != 1L || ar1_phi <= -0.99 || ar1_phi >= 0.99) {
    stop("ar1_phi must be a single numeric value in (-0.99, 0.99).")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) stop("seed must be a single integer.")
    set.seed(as.integer(seed))
  }

  # Time handling
  if (!is.null(time)) {
    if (!is.numeric(time) || anyNA(time) || any(diff(time) <= 0)) {
      stop("time must be numeric, non-NA, and strictly increasing.")
    }
    time_years <- time
  } else {
    time_years <- seq(0, n_years, by = dt)
  }
  n <- length(time_years)

  # Build deterministic backbone
  if (isTRUE(use_custom)) {
    # Validate custom lengths if provided
    if (!is.null(custom_sst) && length(custom_sst) != n) {
      stop("custom_sst length must match length(time).")
    }
    if (!is.null(custom_chl) && length(custom_chl) != n) {
      stop("custom_chl length must match length(time).")
    }

    # Start from either raw custom or base values (if one of custom series is missing)
    sst_values <- if (!is.null(custom_sst)) custom_sst else rep(base_sst, n)
    chl_values <- if (!is.null(custom_chl)) custom_chl else rep(base_chl, n)

    if (!isTRUE(custom_use_raw)) {
      # Add seasonal structure if requested
      if (isTRUE(seasonal)) {
        sst_values <- sst_values + sst_amplitude * sin(2 * pi * time_years)
        chl_values <- chl_values + chl_amplitude * sin(2 * pi * time_years + pi)
      }
    } else {
      # If using raw custom, ignore seasonal/stochastic switches
      seasonal   <- FALSE
      stochastic <- FALSE
    }

  } else {
    # Synthetic backbone (original behavior)
    if (isTRUE(seasonal)) {
      sst_values <- base_sst + sst_amplitude * sin(2 * pi * time_years)
      chl_values <- base_chl + chl_amplitude * sin(2 * pi * time_years + pi)
    } else {
      sst_values <- rep(base_sst, n)
      chl_values <- rep(base_chl, n)
    }
  }

  # Optional stochasticity (applied to current sst_values/chl_values)
  if (isTRUE(stochastic)) {
    # SST: additive Gaussian (optional AR(1))
    if (is.numeric(sst_noise_sd) && length(sst_noise_sd) == 1L && sst_noise_sd > 0) {
      if (ar1_phi == 0) {
        eps_sst <- stats::rnorm(n, mean = 0, sd = sst_noise_sd)
      } else {
        # Simulate AR(1) with unit innovation, then center/scale to target marginal SD
        sst_raw <- as.numeric(stats::arima.sim(model = list(ar = ar1_phi), n = n, sd = 1))
        eps_sst <- (sst_raw - mean(sst_raw)) / stats::sd(sst_raw) * sst_noise_sd
      }
      sst_values <- sst_values + eps_sst
    }

    # Chl: multiplicative lognormal (optional AR(1)), mean-preserving
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

  # Optional bounds to keep within downstream validation ranges
  if (isTRUE(enforce_bounds)) {
    sst_before <- sst_values
    chl_before <- chl_values

    sst_values <- pmin(pmax(sst_values, -2), 35)
    chl_values <- pmin(pmax(chl_values, 0), 50)

    n_clamp_sst <- sum(sst_before != sst_values)
    n_clamp_chl <- sum(chl_before != chl_values)
    if (n_clamp_sst > 0) {
      warning(n_clamp_sst, " SST values clamped to [-2, 35] °C")
    }
    if (n_clamp_chl > 0) {
      warning(n_clamp_chl, " Chl values clamped to [0, 50] mg/m^3")
    }
  }

  return(data.frame(
    time = time_years,
    sst  = sst_values,
    chl  = chl_values
  ))
}
