#' Calculate Size Selectivity for Fishing Gear
#'
#' @description
#' Calculates size-dependent fishing selectivity, representing the probability
#' that an organism of a given size is vulnerable to fishing gear. This function
#' supports multiple selectivity curve types commonly used in fisheries science.
#'
#' @param w Numeric vector of body masses. By default interpreted as log10(mass in grams)
#'   to match ZooMSS conventions, but can be linear scale if `log_scale = FALSE`.
#' @param type Character string specifying selectivity type. Options:
#'   \itemize{
#'     \item `"knife_edge"`: Step function at threshold size (default)
#'     \item `"logistic"`: Smooth logistic curve
#'     \item `"custom"`: User-defined function (must be provided in params$fun)
#'   }
#' @param params Named list of selectivity parameters. Required parameters depend on type:
#'   \itemize{
#'     \item For `"knife_edge"`: `w_min` (minimum vulnerable size)
#'     \item For `"logistic"`: `L50` (size at 50\% selectivity), `L95` (size at 95\% selectivity)
#'     \item For `"custom"`: `fun` (function taking w as argument and returning selectivity)
#'   }
#' @param log_scale Logical indicating whether `w` and size parameters in `params`
#'   are in log10 scale. Default `TRUE` to match ZooMSS size class convention
#'   (`w_log10`). Set to `FALSE` if working with linear body masses.
#'
#' @return Numeric vector of selectivity values between 0 and 1, same length as `w`.
#'
#' @details
#' **Knife-edge selectivity:**
#' All organisms above `w_min` are fully vulnerable (selectivity = 1),
#' all below are invulnerable (selectivity = 0). This represents trawl nets or
#' other gear with a sharp size cutoff.
#'
#' **Logistic selectivity:**
#' Smooth transition from 0 to 1 selectivity, parameterized by `L50` (size at
#' 50\% selectivity) and `L95` (size at 95\% selectivity). The curve follows:
#' \deqn{S(w) = \frac{1}{1 + \exp(-\log(19) \cdot (w - L_{50}) / (L_{95} - L_{50}))}}
#'
#' This represents gillnets or hooks where vulnerability increases gradually with size.
#'
#' **Custom selectivity:**
#' Allows user-defined selectivity functions for specialized gear types.
#' The custom function should accept a numeric vector and return selectivity values.
#'
#' @note
#' When using `log_scale = TRUE` (default), all size parameters (`w_min`, `L50`, `L95`)
#' must be in log10 grams to match the ZooMSS size class grid (`w_log10`).
#'
#' @examples
#' # Knife-edge selectivity (log10 scale, matching ZooMSS convention)
#' w_log10 <- seq(-2, 3, by = 0.1)  # Size classes from 0.01g to 1000g
#' s_knife <- calc_selectivity(w_log10, type = "knife_edge",
#'                              params = list(w_min = 1))  # Vulnerable above 10g
#' plot(w_log10, s_knife, type = "l", xlab = "log10(mass in g)", ylab = "Selectivity")
#'
#' # Logistic selectivity
#' s_logistic <- calc_selectivity(w_log10, type = "logistic",
#'                                 params = list(L50 = 1.5, L95 = 2.0))
#' plot(w_log10, s_logistic, type = "l", xlab = "log10(mass in g)", ylab = "Selectivity")
#'
#' # Custom selectivity (dome-shaped for size-selective gear)
#' dome_selectivity <- function(w) {
#'   exp(-((w - 1.5)^2) / 0.5)  # Peak selectivity at log10(mass) = 1.5
#' }
#' s_custom <- calc_selectivity(w_log10, type = "custom",
#'                               params = list(fun = dome_selectivity))
#' plot(w_log10, s_custom, type = "l", xlab = "log10(mass in g)", ylab = "Selectivity")
#'
#' @export
calc_selectivity <- function(w, type = "knife_edge", params = list(), log_scale = TRUE) {
  
  # Input validation
  if (!is.numeric(w) || length(w) == 0) {
    stop("'w' must be a non-empty numeric vector")
  }
  
  if (!is.character(type) || length(type) != 1) {
    stop("'type' must be a single character string")
  }
  
  if (!is.list(params)) {
    stop("'params' must be a named list")
  }
  
  # Calculate selectivity based on type
  selectivity <- switch(
    type,
    knife_edge = {
      if (!"w_min" %in% names(params)) {
        stop("'params' must contain 'w_min' for knife-edge selectivity")
      }
      selectivity_knife_edge(w, params$w_min)
    },
    logistic = {
      if (!all(c("L50", "L95") %in% names(params))) {
        stop("'params' must contain 'L50' and 'L95' for logistic selectivity")
      }
      if (params$L95 <= params$L50) {
        stop("'L95' must be greater than 'L50'")
      }
      selectivity_logistic(w, params$L50, params$L95)
    },
    custom = {
      if (!"fun" %in% names(params)) {
        stop("'params' must contain 'fun' for custom selectivity")
      }
      if (!is.function(params$fun)) {
        stop("'params$fun' must be a function")
      }
      params$fun(w)
    },
    stop("Unknown selectivity type '", type, "'. Must be 'knife_edge', 'logistic', or 'custom'")
  )
  
  # Validate output
  if (!is.numeric(selectivity) || length(selectivity) != length(w)) {
    stop("Selectivity function must return a numeric vector of the same length as 'w'")
  }
  
  if (any(selectivity < 0 | selectivity > 1, na.rm = TRUE)) {
    warning("Selectivity values outside [0,1] detected. Clamping to valid range.")
    selectivity <- pmax(0, pmin(1, selectivity))
  }
  
  return(selectivity)
}


#' Knife-Edge Selectivity Function
#'
#' Internal helper function implementing knife-edge (step function) selectivity.
#' Not exported - used internally by calc_selectivity.
#'
#' @param w Numeric vector of body masses (or log10 body masses)
#' @param w_min Numeric scalar, minimum vulnerable size (same scale as w)
#'
#' @return Numeric vector of selectivity values (0 or 1)
#'
#' @details
#' Implements a step function where selectivity is 0 for w less than w_min
#' (invulnerable) and 1 for w greater than or equal to w_min (fully vulnerable).
#'
#' @keywords internal
selectivity_knife_edge <- function(w, w_min) {
  as.numeric(w >= w_min)
}


#' Logistic Selectivity Function
#'
#' Internal helper function implementing logistic (sigmoid) selectivity curve.
#' Not exported - used internally by calc_selectivity.
#'
#' @param w Numeric vector of body masses (or log10 body masses)
#' @param L50 Numeric scalar, size at 50 percent selectivity (same scale as w)
#' @param L95 Numeric scalar, size at 95 percent selectivity (same scale as w)
#'
#' @return Numeric vector of selectivity values between 0 and 1
#'
#' @details
#' Implements a logistic curve parameterized by two quantiles.
#' The coefficient log(19) ensures that selectivity equals 0.95 at L95
#' when it equals 0.5 at L50, following standard fisheries parameterization.
#'
#' This produces a smooth S-shaped curve representing gradual increase in
#' vulnerability with size, typical of gillnets and hook-based gear.
#'
#' @keywords internal
selectivity_logistic <- function(w, L50, L95) {
  # Logistic curve ensuring S(L50) = 0.5 and S(L95) = 0.95
  # Coefficient log(19) derived from: 0.95 / (1-0.95) = 19
  steepness <- log(19) / (L95 - L50)
  1 / (1 + exp(-steepness * (w - L50)))
}
