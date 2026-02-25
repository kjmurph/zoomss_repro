#' Set Up ZooMSS Model Parameters
#'
#' @title Initialize and validate ZooMSS model parameters
#' @description Sets up the complete parameter list for ZooMSS model runs, including
#'   functional group parameters, model dimensions, and environmental forcing data.
#' @details This function creates a comprehensive parameter object that contains:
#'
#'   **Static Parameters (fixed across time steps):**
#'   - Model dimensions (number of groups, size classes, time steps)
#'   - Biological parameters (growth efficiency, mortality rates)
#'   - Size class definitions and ranges for each functional group
#'   - Phytoplankton size spectrum parameters
#'
#'   **Dynamic Parameters (calculated from environmental data):**
#'   - Phytoplankton abundance time series based on chlorophyll
#'   - Temperature effects on metabolism for zooplankton and fish
#'   - Environmental forcing validation and interpolation
#'
#'   The function validates that environmental time series data covers the full
#'   simulation period and pre-calculates time-varying parameters to optimize
#'   model performance during the main simulation loop.
#'
#' @param Groups Data frame containing functional group definitions with columns:
#'   Species, Type, W0 (log min size), Wmax (log max size), and various biological parameters
#' @param input_params Data frame with model parameters including:
#'   time (time vector in years), sst (sea surface temperature), and chl (chlorophyll).
#'   The time vector can start at any value and the model automatically calculates dt (time step) and tmax (maximum time).
#' @param isave Save frequency in time steps (default: 50)
#'
#' @return List containing comprehensive model parameters:
#'   \itemize{
#'     \item Groups: Functional group definitions
#'     \item ngrps: Number of functional groups
#'     \item ngrid: Number of size classes
#'     \item w: Size class weights (g)
#'     \item tmax, dt, isave: Temporal parameters
#'     \item zoo_grps, fish_grps: Indices for different organism types
#'     \item phyto_int, phyto_slope: Time series of phytoplankton parameters
#'     \item temp_eff_zoo, temp_eff_fish: Time series of temperature effects
#'     \item Additional biological and physical parameters
#'   }
#'
#' @examples
#' \dontrun{
#' # Load functional groups
#' data(Groups)
#'
#' # Create environmental time series
#' env_data <- createEnviroData(10, 0.01)
#' input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
#'
#' # Generate parameter list
#' params <- zoomss_params(Groups, input_params, isave = 50)
#' }
#'
#' @noRd
#'
zoomss_params <- function(Groups, input_params, isave){

  # Calculate dt and tmax from time column in input_params
  time_values <- input_params$time
  dt_calc <- time_values[2] - time_values[1]  # Use first time step difference
  tmax_calc <- max(time_values)  # Maximum time value (not duration)

  # Check if time steps are uniform - ERROR if not consistent
  if (length(time_values) > 2) {
    dt_diffs <- diff(time_values)
    max_diff <- max(abs(dt_diffs - dt_calc))
    if (max_diff > dt_calc * 0.001) {  # Allow only 0.1% variation
      stop("Time steps are not uniform in input_params. Maximum deviation: ",
           round(max_diff, 6), " (", round(100 * max_diff / dt_calc, 2), "% of dt). ",
           "ZooMSS requires uniform time steps for accurate results.")
    }
  }

  param <- list(
    Groups = Groups, # Read in functional group specific parameters from file
    ngrps = dim(Groups)[1], # no. of Groups
    dx = 0.1, # log10 weight step
    day = 12, # day length (hours of each day in sun)
    tmax = tmax_calc, # max years - calculated from time
    dt = dt_calc, # timestep - calculated from time
    w0 = 10^(min(Groups$W0)),		# minimum size class
    wMax = 10^(max(Groups$Wmax)),# maximum size class
    # ZSpre = 1, # senescence mortality prefactor
    # ZSexp = 0.3, # senescence mortality exponent
    w0_phyto = 10^(-14.5), # minimum phytoplankton size class (1um)
    # wMax_phyto will be calculated from time series
    zoo_grps = which(Groups$Type == "Zooplankton"), # Which rows are zooplankton
    fish_grps = which(Groups$Type == "Fish"), # Which rows are fish
    num_zoo = sum(Groups$Type == "Zooplankton"), # How many zooplankton
    num_fish = sum(Groups$Type == "Fish"), # How many fish
    cc_phyto = 0.1, # Carbon content of phytoplankton size classes
    Carbon_max = max(Groups$Carbon), # Maximum Carbon content across all groups (for defecation scaling)
    isave = isave # how often to save results every 'isave' time steps
  )

  # Process provided time series
  n_time_steps <- nrow(input_params)

  ## Add additional parameters which are based on the parameter set
  param2 <- list(
    nsave  = max(1, floor(n_time_steps / param$isave)), # Number of time slots to save (regular intervals only, minimum 1)
    ntime = n_time_steps, # Total number of time steps to simulate
    itimemax = n_time_steps # Number of iterations for time loop (one for each time point)
  )
  # Pre-calculate phytoplankton parameters for each time step

  # Initialize arrays to store time-varying parameters (optimize memory allocation)
  param2$phyto_int <- numeric(n_time_steps)
  param2$phyto_slope <- numeric(n_time_steps)
  param2$phyto_max <- numeric(n_time_steps)
  param2$wMax_phyto <- numeric(n_time_steps)

  # Pre-calculate temperature effects for each time step (more efficient than in-loop)
  param2$temp_eff_zoo <- matrix(NA, nrow = n_time_steps, ncol = param$num_zoo)
  param2$temp_eff_fish <- matrix(NA, nrow = n_time_steps, ncol = param$num_fish)

  # Check if phytoplankton parameters are already in input_params (expanded format)
  if (nrow(input_params) == n_time_steps && all(c("phyto_int", "phyto_slope", "phyto_max") %in% names(input_params))) {
    cat("Using pre-calculated phytoplankton parameters from input_params\n")

    # Use pre-calculated values (most efficient path)
    param2$phyto_int <- input_params$phyto_int
    param2$phyto_slope <- input_params$phyto_slope
    param2$phyto_max <- input_params$phyto_max
    param2$wMax_phyto <- 10^input_params$phyto_max

    # Calculate temperature effects for each time step - consistent temperature formula throughout model
    temp_factor <- 2^((input_params$sst - 30)/10)
    param2$temp_eff_zoo[,] <- temp_factor # preallocated above for all groups
    param2$temp_eff_fish[,] <- temp_factor # preallocated above for all groups

  } else {
    cat("Calculating phytoplankton parameters from environmental time series\n")

    # Calculate phytoplankton parameters for each time step

    # Calculate phytoplankton parameters using the existing function
    phyto_params <- calculatePhytoParam(input_params)

    param2$phyto_int <- phyto_params$phyto_int
    param2$phyto_slope <- phyto_params$phyto_slope
    param2$phyto_max <- phyto_params$phyto_max
    param2$wMax_phyto <- 10^phyto_params$phyto_max

    # Temperature effects (same calculation throughout model)
    temp_factor <- 2^((input_params$sst - 30)/10)
    param2$temp_eff_zoo[,] <- temp_factor # preallocated above for all groups
    param2$temp_eff_fish[,] <- temp_factor # preallocated above for all groups
  }

  # Store the maximum wMax_phyto from time series for grid creation
  param2$wMax_phyto <- max(param2$wMax_phyto)

  # Add final parameters that depend on the complete parameter set (calculate only once)
  param2$w_log10 <- round(seq(from = min(Groups$W0), to = max(Groups$Wmax), param$dx), digits = 2) # Set up log10 weight grid
  param2$w <- 10^(seq(from = min(Groups$W0), to = max(Groups$Wmax), param$dx)) # Set up weight grid
  param2$w_phyto <- 10^(seq(from = log10(param$w0_phyto), to = log10(param2$wMax_phyto), param$dx)) # Set up phytoplankton size classes
  param2$ngrid <- length(param2$w) # total number of size classes for zoo and fish
  param2$ngridPP <- length(param2$w_phyto) # total number of size classes for phyto

  # =============================================================================
  # FISH ENERGY BUDGET DERIVED PARAMETERS
  # =============================================================================
  # Only fish groups use the explicit energy budget.
  # Zooplankton use the original GGE pathway (GrossGEscale * Carbon).

  # Derived reproduction fraction for fish (R_frac = 1 - f_M - K_growth)
  # Set to 0 for zooplankton (they don't use this pathway)
  param2$R_frac <- rep(0, param$ngrps)
  names(param2$R_frac) <- Groups$Species

  for (f in seq_along(param$fish_grps)) {
    fg <- param$fish_grps[f]
    param2$R_frac[fg] <- 1 - Groups$f_M[fg] - Groups$K_growth[fg]
  }

  # Validate fish energy budget closure
  for (f in seq_along(param$fish_grps)) {
    fg <- param$fish_grps[f]
    budget_sum <- Groups$f_M[fg] + Groups$K_growth[fg] + param2$R_frac[fg]
    if (abs(budget_sum - 1.0) > 1e-10) {
      stop(paste("Fish energy budget does not sum to 1.0 for:", Groups$Species[fg]))
    }
    if (param2$R_frac[fg] < 0) {
      stop(paste("R_frac is negative (f_M + K_growth > 1) for:", Groups$Species[fg]))
    }
  }

  cat("Dual-pathway model:\n")
  cat("  Zooplankton: GGE pathway (GrossGEscale * Carbon)\n")
  cat("  Fish: Explicit energy budget\n")
  if (param$num_fish > 0) {
    cat("    f_M range:      ", range(Groups$f_M[param$fish_grps]), "\n")
    cat("    K_growth range:  ", range(Groups$K_growth[param$fish_grps]), "\n")
    cat("    R_frac range:    ", range(param2$R_frac[param$fish_grps]), "\n")
  }

  # Phytoplankton defecation for fish groups — continuous scaling based on cc_phyto
  # def = def_high + (def_low - def_high) * (1 - Carbon / Carbon_max)
  # Use fish group values (they should be identical across fish)
  if (param$num_fish > 0) {
    fg1 <- param$fish_grps[1]
    param2$def_phyto <- Groups$def_high[fg1] +
      (Groups$def_low[fg1] - Groups$def_high[fg1]) * (1 - param$cc_phyto / param$Carbon_max)
  } else {
    param2$def_phyto <- 0.3  # fallback
  }

  # Maturation size indices for each group (index of Wmat in w_log10 grid)
  param2$mat_size_idx <- sapply(1:param$ngrps, function(g) {
    which.min(abs(param2$w_log10 - Groups$Wmat[g]))
  })

  # Minimum size class index for each group (for recruitment boundary)
  param2$min_size_idx <- sapply(1:param$ngrps, function(g) {
    which(param2$w_log10 == Groups$W0[g])
  })

  # Maximum size class index for each group
  param2$max_size_idx <- sapply(1:param$ngrps, function(g) {
    which(param2$w_log10 == Groups$Wmax[g])
  })

  # Create maturity ogive matrix (ngrps x ngrid)
  # mat_ogive = 1 / (1 + exp(-slope * (log10_w - Wmat)))
  # Gives smooth transition from 0 (fully immature) to 1 (fully mature) around Wmat
  param2$mat_ogive <- matrix(0, nrow = param$ngrps, ncol = param2$ngrid)
  for (g in 1:param$ngrps) {
    param2$mat_ogive[g, ] <- 1 / (1 + exp(-Groups$mat_ogive_slope[g] *
                                            (param2$w_log10 - Groups$Wmat[g])))
    # Zero out outside group's size range
    param2$mat_ogive[g, param2$w_log10 < Groups$W0[g]] <- 0
    param2$mat_ogive[g, param2$w_log10 > Groups$Wmax[g]] <- 0
  }

  # =============================================================================
  # EFFORT-DRIVEN FISHING PARAMETERS
  # =============================================================================
  # Detect effort time series columns dynamically based on fish groups present.
  # Expected column naming convention: effort_{Species} where Species is the
  # fish group name converted to lowercase with spaces replaced by underscores.
  # e.g., Fish_Small -> effort_fish_small, Fish_Med -> effort_fish_med
  #
  # When present, fishing mortality is calculated dynamically each time step as:
  #   F(w,t) = Effort(t) * q * Selectivity(w)
  # When absent, the model falls back to the existing static Fmort pathway.

  if (param$num_fish > 0) {
    # Build expected column names from fish group Species names
    fish_species <- Groups$Species[param$fish_grps]
    effort_cols <- paste0("effort_", tolower(gsub(" ", "_", fish_species)))
    names(effort_cols) <- fish_species

    effort_available <- all(effort_cols %in% names(input_params))
  } else {
    effort_cols <- character(0)
    effort_available <- FALSE
  }

  if (effort_available) {
    cat("Effort time series detected — enabling effort-driven fishing\n")
    cat("  Fish groups:", paste(fish_species, collapse = ", "), "\n")
    cat("  Effort columns:", paste(effort_cols, collapse = ", "), "\n")

    # Validate q column exists in Groups
    if (!"q" %in% names(Groups)) {
      stop("Groups must contain a 'q' (catchability) column when effort data is provided.\n",
           "  Add q values to Groups or use setFishingParams() to set up fishing parameters.")
    }

    # Store effort matrix (ntime x num_fish)
    param2$effort <- as.matrix(input_params[, effort_cols, drop = FALSE])
    colnames(param2$effort) <- effort_cols

    # Store catchability per fish group (vector, length = num_fish)
    param2$q <- Groups$q[param$fish_grps]
    names(param2$q) <- fish_species

    # Pre-calculate selectivity vectors per fish group (knife-edge at Fmort_W0)
    # selectivity[f, ] = 1 where w_log10 >= Fmort_W0 AND w_log10 <= Fmort_Wmax, else 0
    # Stored as full ngrps x ngrid matrix (zeros for non-fish groups)
    param2$selectivity <- matrix(0, nrow = param$ngrps, ncol = param2$ngrid)
    for (f in 1:param$num_fish) {
      fg <- param$fish_grps[f]
      sel_idx <- which(param2$w_log10 >= Groups$Fmort_W0[fg] &
                         param2$w_log10 <= Groups$Fmort_Wmax[fg])
      param2$selectivity[fg, sel_idx] <- 1
    }

    param2$effort_fishing <- TRUE
    cat("  Catchability (q):", paste(names(param2$q), "=", param2$q, collapse = ", "), "\n")

  } else {
    param2$effort_fishing <- FALSE
    # Static Fmort pathway will be used (existing behaviour)
    if (param$num_fish > 0 && any(Groups$Fmort[param$fish_grps] > 0)) {
      cat("Using static fishing mortality (no effort time series)\n")
    }
  }

  # Final parameter combination
  # Exclude time series vectors from input_params since they're now stored as _ts arrays in param2
  input_params_filtered <- input_params[!names(input_params) %in% c("tmax", "dt", "isave", "time_step", "phyto_int", "phyto_slope", "phyto_max",
                                                                    effort_cols)]

  param_final <- c(input_params_filtered, param, param2)
  return(param_final)
}
