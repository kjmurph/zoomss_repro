#' Set Fishing Parameters for Effort-Driven Fishing
#'
#' @title Configure catchability and selectivity parameters for effort-driven fishing
#' @description Adds catchability coefficient (q) to the Groups data frame and optionally
#'   updates selectivity bounds (Fmort_W0, Fmort_Wmax). This function prepares the Groups
#'   data for use with effort time series in zoomss_model().
#'
#' @details When effort time series are provided in input_params, fishing mortality is
#'   calculated dynamically at each time step as:
#'
#'   F(w, t) = Effort(t) * q * Selectivity(w)
#'
#'   where Selectivity is a knife-edge function: 1 for w in [Fmort_W0, Fmort_Wmax], 0 otherwise.
#'
#'   The catchability coefficient q controls the efficiency of fishing effort at removing
#'   biomass. It is the primary parameter to be estimated during calibration against
#'   observed catch data.
#'
#' @param Groups Data frame of functional groups (from getGroups())
#' @param q_small Numeric. Catchability for Fish_Small group (default: 0)
#' @param q_med Numeric. Catchability for Fish_Med group (default: 0)
#' @param q_large Numeric. Catchability for Fish_Large group (default: 0)
#' @param w_min_small Numeric. Log10 minimum weight (g) for Fish_Small selectivity.
#'   If NULL (default), uses existing Fmort_W0 value.
#' @param w_min_med Numeric. Log10 minimum weight (g) for Fish_Med selectivity.
#'   If NULL (default), uses existing Fmort_W0 value.
#' @param w_min_large Numeric. Log10 minimum weight (g) for Fish_Large selectivity.
#'   If NULL (default), uses existing Fmort_W0 value.
#'
#' @return Modified Groups data frame with q column added/updated and optionally
#'   updated Fmort_W0 values
#' @export
#'
#' @examples
#' \dontrun{
#' Groups <- getGroups()
#' Groups <- setFishingParams(Groups, q_small = 0.01, q_med = 0.02, q_large = 0.03)
#'
#' # With custom selectivity thresholds
#' Groups <- setFishingParams(Groups,
#'                            q_small = 0.01, q_med = 0.02, q_large = 0.03,
#'                            w_min_small = 0.5, w_min_med = 1.5, w_min_large = 3.0)
#' }
#'
setFishingParams <- function(Groups,
                             q_small = 0, q_med = 0, q_large = 0,
                             w_min_small = NULL, w_min_med = NULL, w_min_large = NULL) {

  # Validate Groups has fish
  fish_mask <- Groups$Type == "Fish"
  if (sum(fish_mask) == 0) {
    stop("No fish groups found in Groups data frame")
  }

  # Initialize q column (0 for all groups — zooplankton are never fished in this framework)
  if (!"q" %in% names(Groups)) {
    Groups$q <- 0
  }

  # Map q values to fish groups by name
  fish_names <- Groups$Species[fish_mask]
  q_values <- c(q_small, q_med, q_large)

  if (length(fish_names) != 3) {
    warning("Expected 3 fish groups (Fish_Small, Fish_Med, Fish_Large), found ",
            length(fish_names), ". Assigning q values in order.")
  }

  n_fish <- min(length(fish_names), length(q_values))
  fish_idx <- which(fish_mask)
  for (i in 1:n_fish) {
    Groups$q[fish_idx[i]] <- q_values[i]
  }

  # Optionally update selectivity bounds (Fmort_W0)
  w_mins <- list(w_min_small, w_min_med, w_min_large)
  for (i in 1:n_fish) {
    if (!is.null(w_mins[[i]])) {
      Groups$Fmort_W0[fish_idx[i]] <- w_mins[[i]]
    }
  }

  # Validate selectivity bounds
  for (i in 1:n_fish) {
    fg <- fish_idx[i]
    if (Groups$Fmort_W0[fg] < Groups$W0[fg] || Groups$Fmort_W0[fg] > Groups$Wmax[fg]) {
      warning("Fmort_W0 for ", Groups$Species[fg], " (", Groups$Fmort_W0[fg],
              ") is outside group size range [", Groups$W0[fg], ", ", Groups$Wmax[fg], "]")
    }
    if (Groups$Fmort_Wmax[fg] < Groups$Fmort_W0[fg]) {
      warning("Fmort_Wmax < Fmort_W0 for ", Groups$Species[fg])
    }
  }

  cat("Fishing parameters set:\n")
  for (i in 1:n_fish) {
    fg <- fish_idx[i]
    cat("  ", Groups$Species[fg], ": q =", Groups$q[fg],
        ", selectivity range [", Groups$Fmort_W0[fg], ",", Groups$Fmort_Wmax[fg], "] log10 g\n")
  }

  return(Groups)
}


#' Add Fishing Effort Time Series to Input Parameters
#'
#' @title Append effort columns to input_params for effort-driven fishing
#' @description Adds effort_small, effort_med, and effort_large columns to the
#'   input_params data frame. When these columns are present, zoomss_model()
#'   automatically enables effort-driven fishing mortality.
#'
#' @details Effort values are dimensionless multipliers that combine with catchability (q)
#'   and selectivity to determine fishing mortality: F(w,t) = Effort(t) * q * Selectivity(w).
#'
#'   Effort time series must have the same number of rows as input_params (one value per
#'   time step). Values should be non-negative. Zero effort means no fishing for that
#'   group at that time step.
#'
#' @param input_params Data frame created by createInputParams() with time, sst, chl columns
#' @param effort_small Numeric vector. Effort time series for Fish_Small group.
#'   If scalar, replicated to all time steps.
#' @param effort_med Numeric vector. Effort time series for Fish_Med group.
#'   If scalar, replicated to all time steps.
#' @param effort_large Numeric vector. Effort time series for Fish_Large group.
#'   If scalar, replicated to all time steps.
#'
#' @return Modified input_params data frame with effort columns added
#' @export
#'
#' @examples
#' \dontrun{
#' # Constant effort
#' input_params <- addFishingEffort(input_params, effort_small = 1, effort_med = 1, effort_large = 1)
#'
#' # Time-varying effort (e.g., from FishMIP)
#' input_params <- addFishingEffort(input_params,
#'                                  effort_small = fishmip_effort$small,
#'                                  effort_med = fishmip_effort$med,
#'                                  effort_large = fishmip_effort$large)
#' }
#'
addFishingEffort <- function(input_params, effort_small = 0, effort_med = 0, effort_large = 0) {

  n <- nrow(input_params)

  # Helper: expand scalar to vector
  expand_effort <- function(eff, name) {
    if (length(eff) == 1) {
      rep(eff, n)
    } else if (length(eff) == n) {
      eff
    } else {
      stop(name, " must be either length 1 (constant) or length ", n,
           " (matching input_params rows). Got length ", length(eff))
    }
  }

  input_params$effort_small <- expand_effort(effort_small, "effort_small")
  input_params$effort_med   <- expand_effort(effort_med, "effort_med")
  input_params$effort_large <- expand_effort(effort_large, "effort_large")

  # Validate non-negative
  if (any(input_params$effort_small < 0) || any(input_params$effort_med < 0) ||
      any(input_params$effort_large < 0)) {
    stop("Effort values must be non-negative")
  }

  cat("Fishing effort added to input_params:\n")
  cat("  effort_small: range [", range(input_params$effort_small), "]\n")
  cat("  effort_med:   range [", range(input_params$effort_med), "]\n")
  cat("  effort_large: range [", range(input_params$effort_large), "]\n")

  return(input_params)
}
