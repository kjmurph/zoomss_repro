#' Calculate Fishing Mortality from Effort and Selectivity
#'
#' @description
#' Calculates size- and group-specific fishing mortality rates from fishing effort,
#' catchability coefficients, and selectivity curves. This function implements the
#' standard fisheries equation: F(w,g) = E(t) × q(g) × S(w,g)
#'
#' @param effort Numeric scalar, fishing effort at current time step. Can represent
#'   various effort metrics (vessel-days, trawl-hours, hook-hours, etc.) depending
#'   on the application. Units should be consistent with catchability coefficient.
#' @param q Numeric vector of length ngrps, catchability coefficient for each
#'   functional group (units: 1/yr per unit effort). Represents the efficiency
#'   with which effort translates to fishing mortality.
#' @param selectivity Matrix of dimension (ngrps x ngrid), size-specific selectivity
#'   for each functional group. Values should be between 0 and 1, where 0 = invulnerable
#'   and 1 = fully vulnerable. Typically generated using calc_selectivity.
#'
#' @return Matrix of dimension (ngrps × ngrid) containing fishing mortality rates
#'   (yr⁻¹) for each functional group and size class.
#'
#' @details
#' **Fishing mortality equation:**
#' \deqn{F(w,g,t) = E(t) \times q(g) \times S(w,g)}
#'
#' where:
#' \itemize{
#'   \item \eqn{E(t)} = fishing effort at time t
#'   \item \eqn{q(g)} = catchability of functional group g
#'   \item \eqn{S(w,g)} = selectivity at body mass w for group g
#' }
#'
#' This formulation allows fishing mortality to vary:
#' \itemize{
#'   \item **Over time:** through changes in effort E(t)
#'   \item **Among groups:** through different catchability coefficients q(g)
#'   \item **With body size:** through selectivity curves S(w,g)
#' }
#'
#' **Relationship to total mortality:**
#' Fishing mortality is combined with natural mortality components:
#' \deqn{Z(w,g,t) = M_{senescence}(w,g) + M_{predation}(w,g,t) + F(w,g,t)}
#'
#' @note
#' \itemize{
#'   \item When effort = 0, fishing mortality is zero (unfished scenario)
#'   \item Catchability q should be calibrated to match observed catches
#'   \item Selectivity should be pre-calculated using \code{\link{calc_selectivity}}
#'     for each functional group
#' }
#'
#' @seealso
#' \code{\link{calc_selectivity}} for generating selectivity matrices,
#' \code{\link{calc_catch}} for calculating resulting catch
#'
#' @examples
#' # Set up parameters for 3 functional groups, 50 size classes
#' ngrps <- 3
#' ngrid <- 50
#' w_log10 <- seq(-2, 3, length.out = ngrid)
#'
#' # Define catchability (different vulnerability by group)
#' q <- c(0.01, 0.02, 0.015)  # Small, Medium, Large fish
#'
#' # Create selectivity matrix (knife-edge at different sizes)
#' selectivity <- matrix(0, nrow = ngrps, ncol = ngrid)
#' selectivity[1, ] <- calc_selectivity(w_log10, "knife_edge", list(w_min = 0.5))
#' selectivity[2, ] <- calc_selectivity(w_log10, "knife_edge", list(w_min = 1.0))
#' selectivity[3, ] <- calc_selectivity(w_log10, "knife_edge", list(w_min = 1.5))
#'
#' # Calculate fishing mortality for current effort level
#' effort_current <- 100  # e.g., 100 vessel-days
#' F_mortality <- calc_fishing_mortality(effort_current, q, selectivity)
#'
#' # Examine resulting F for each group
#' matplot(w_log10, t(F_mortality), type = "l", lty = 1,
#'         xlab = "log10(mass in g)", ylab = "Fishing mortality (1/yr)",
#'         main = "Size-specific fishing mortality by group")
#' legend("topleft", legend = c("Small", "Medium", "Large"), col = 1:3, lty = 1)
#'
#' @export
calc_fishing_mortality <- function(effort, q, selectivity) {
  
  # Input validation
  if (!is.numeric(effort) || length(effort) != 1) {
    stop("'effort' must be a single numeric value")
  }
  
  if (effort < 0) {
    stop("'effort' must be non-negative")
  }
  
  if (!is.numeric(q) || length(q) == 0) {
    stop("'q' must be a numeric vector of catchability coefficients")
  }
  
  if (any(q < 0)) {
    stop("'q' must contain non-negative values")
  }
  
  if (!is.matrix(selectivity)) {
    stop("'selectivity' must be a matrix")
  }
  
  if (nrow(selectivity) != length(q)) {
    stop("'selectivity' must have nrow = length(q) (one row per group)")
  }
  
  if (any(selectivity < 0 | selectivity > 1, na.rm = TRUE)) {
    warning("'selectivity' contains values outside [0, 1]. Check selectivity calculation.")
  }
  
  # Calculate F(w,g) = E(t) * q(g) * S(w,g)
  # Use sweep to efficiently multiply each row by effort * q[g]
  F_matrix <- sweep(selectivity, 1, effort * q, "*")
  
  return(F_matrix)
}


#' Calculate Catch from Fishing Mortality and Abundance
#'
#' @description
#' Calculates catch biomass (yield) from fishing mortality, abundance, and body size
#' using the Baranov catch equation. This is the biomass removed from the population
#' by fishing.
#'
#' @param F_w Matrix of dimension (ngrps × ngrid), fishing mortality rates (yr⁻¹)
#'   for each functional group and size class. Typically from \code{\link{calc_fishing_mortality}}.
#' @param N_w Matrix of dimension (ngrps × ngrid), abundance density (numbers per m²)
#'   for each functional group and size class. This is the population density
#'   at the current time step.
#' @param w Numeric vector of length ngrid, body mass of each size class (grams).
#'   Should match the size classes used in ZooMSS (\code{param$w}).
#' @param dw Numeric scalar, width of size classes in log10 space. Default is 0.1
#'   to match ZooMSS convention (\code{param$dx}).
#' @param dt Numeric scalar, time step in years. Default is 1.0 for annual catch.
#'   Use smaller values (e.g., 1/12 for monthly time steps) if needed.
#' @param by_size Logical, if TRUE returns catch by size class (matrix), if FALSE
#'   returns total catch per group (vector). Default FALSE.
#'
#' @return
#' If \code{by_size = FALSE} (default): Numeric vector of length ngrps containing
#' total catch per functional group (g m⁻² yr⁻¹).
#'
#' If \code{by_size = TRUE}: Matrix of dimension (ngrps × ngrid) containing
#' size-resolved catch for each group (g m⁻² yr⁻¹).
#'
#' @details
#' **Catch equation:**
#' \deqn{C(w,g) = F(w,g) \times N(w,g) \times w \times \Delta w \times \Delta t}
#'
#' where:
#' \itemize{
#'   \item \eqn{F(w,g)} = fishing mortality rate (yr⁻¹)
#'   \item \eqn{N(w,g)} = abundance density (numbers m⁻²)
#'   \item \eqn{w} = individual body mass (g)
#'   \item \eqn{\Delta w} = size class width (accounting for log10 spacing)
#'   \item \eqn{\Delta t} = time step duration (yr)
#' }
#'
#' **Size class integration:**
#' Since ZooMSS uses log10-spaced size classes, the width in linear space is:
#' \deqn{\Delta w_{linear} = w \times \log(10) \times dw}
#'
#' **Total catch per group:**
#' \deqn{C_{total}(g) = \sum_{w} C(w,g)}
#'
#' This represents the instantaneous catch rate. For finite time periods with
#' changing populations, the Baranov catch equation would require integration,
#' but this approximation is valid for small time steps or low F values.
#'
#' @note
#' \itemize{
#'   \item Catch units are biomass per area per time (g m⁻² yr⁻¹)
#'   \item To convert to total regional catch, multiply by region area
#'   \item This gives instantaneous catch; for seasonal or pulse fishing,
#'     adjust \code{dt} accordingly
#'   \item Returns zero when F = 0 (unfished) or N = 0 (extinct)
#' }
#'
#' @seealso
#' \code{\link{calc_fishing_mortality}} for calculating F(w,g)
#'
#' @examples
#' # Set up example population
#' ngrps <- 3
#' ngrid <- 50
#' w <- 10^seq(-2, 3, length.out = ngrid)  # Body masses from 0.01g to 1000g
#' w_log10 <- log10(w)
#'
#' # Create abundance distribution (decreasing with size)
#' N <- matrix(0, nrow = ngrps, ncol = ngrid)
#' for (g in 1:ngrps) {
#'   N[g, ] <- 1e6 * w^(-2)  # Power-law size spectrum
#' }
#'
#' # Set up fishing
#' q <- c(0.01, 0.02, 0.015)
#' selectivity <- matrix(0, nrow = ngrps, ncol = ngrid)
#' for (g in 1:ngrps) {
#'   selectivity[g, ] <- calc_selectivity(w_log10, "knife_edge",
#'                                        list(w_min = g * 0.5))
#' }
#'
#' # Calculate fishing mortality and catch
#' effort <- 100
#' F_w <- calc_fishing_mortality(effort, q, selectivity)
#' catch_total <- calc_catch(F_w, N, w, dw = 0.1)
#'
#' print(catch_total)  # Total catch per group
#'
#' # Get size-resolved catch
#' catch_by_size <- calc_catch(F_w, N, w, dw = 0.1, by_size = TRUE)
#' matplot(w_log10, t(catch_by_size), type = "l", lty = 1, log = "y",
#'         xlab = "log10(mass in g)", ylab = "Catch (g/m²/yr)",
#'         main = "Size-resolved catch by group")
#'
#' @export
calc_catch <- function(F_w, N_w, w, dw = 0.1, dt = 1.0, by_size = FALSE) {
  
  # Input validation
  if (!is.matrix(F_w)) {
    stop("'F_w' must be a matrix")
  }
  
  if (!is.matrix(N_w)) {
    stop("'N_w' must be a matrix")
  }
  
  if (!identical(dim(F_w), dim(N_w))) {
    stop("'F_w' and 'N_w' must have identical dimensions")
  }
  
  if (!is.numeric(w) || length(w) != ncol(F_w)) {
    stop("'w' must be a numeric vector with length = ncol(F_w)")
  }
  
  if (any(w <= 0)) {
    stop("'w' must contain positive values (body masses)")
  }
  
  if (!is.numeric(dw) || length(dw) != 1 || dw <= 0) {
    stop("'dw' must be a single positive number")
  }
  
  if (!is.numeric(dt) || length(dt) != 1 || dt <= 0) {
    stop("'dt' must be a single positive number")
  }
  
  if (any(F_w < 0, na.rm = TRUE)) {
    warning("Negative fishing mortality detected. Setting to zero.")
    F_w[F_w < 0] <- 0
  }
  
  if (any(N_w < 0, na.rm = TRUE)) {
    warning("Negative abundance detected. Setting to zero.")
    N_w[N_w < 0] <- 0
  }
  
  # Calculate catch per size class
  # C(w,g) = F(w,g) * N(w,g) * w * Δw * Δt
  # Note: Δw in log space, so actual width = w * log(10) * dw
  
  # Replicate w across rows to match matrix dimensions
  w_matrix <- matrix(w, nrow = nrow(F_w), ncol = ncol(F_w), byrow = TRUE)
  
  # Catch biomass per size class (g m⁻² yr⁻¹, or per time step if dt ≠ 1)
  catch_matrix <- F_w * N_w * w_matrix * log(10) * dw * dt
  
  if (by_size) {
    return(catch_matrix)
  } else {
    # Sum across size classes to get total catch per group
    catch_total <- rowSums(catch_matrix, na.rm = TRUE)
    return(catch_total)
  }
}
