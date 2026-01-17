#' Run Complete ZooMSS Model Simulation
#'
#' @title Main ZooMSS model function for complete simulations
#' @description This is the main wrapper function that orchestrates a complete ZooMSS
#'   model simulation from parameter setup through model execution to output processing.
#' @details This function coordinates the entire ZooMSS modeling workflow:
#'   1. Validates that environmental time series data is provided
#'   2. Sets up model parameters using the Groups data and input parameters
#'   3. Initializes the model structure and feeding kernels
#'   4. Runs the model forward in time with dynamic environmental forcing
#'   5. Processes outputs by averaging the final 50% of the simulation
#'   6. Returns organized results including abundances, diets, growth, and mortality
#'
#'   This is the primary entry point for
#'   running ZooMSS simulations with environmental forcing.
#'
#' @param input_params Data frame containing model parameters and environmental time series.
#'   Must include columns: time (time vector in years), sst (sea surface temperature),
#'   and chl (chlorophyll). Can optionally include cellID for spatial data. The time step (dt)
#'   and maximum time (tmax) are automatically calculated from the time vector. Can be created using createInputParams().
#' @param Groups Data frame defining functional groups with their biological parameters.
#'   Must include columns defining species characteristics, size ranges, and feeding parameters.
#'   If NULL, uses default ZooMSS functional groups. Can be obtained/customized using
#'   getGroups().
#' @param isave Save frequency in time steps (default: 10)
#' @param fishing_params Optional list specifying dynamic fishing mortality. If NULL (default),
#'   uses static fishing mortality from Groups$Fmort. If provided, must contain:
#'   \itemize{
#'     \item effort: numeric vector of fishing effort time series (length = nrow(input_params))
#'     \item catchability: named numeric vector of catchability coefficients (q) for each group
#'     \item selectivity: list with named elements for each group, each containing:
#'       \itemize{
#'         \item type: "knife_edge", "logistic", or "custom"
#'         \item params: list of parameters (e.g., list(w_min = 1.0) for knife_edge)
#'         \item log_scale: TRUE (default) if params are in log10 grams
#'       }
#'   }
#' @param save_catch_by_size Logical, if TRUE saves size-resolved catch (nsave x ngrps x ngrid),
#'   if FALSE (default) saves only total catch by group (nsave x ngrps). Size-resolved catch
#'   provides detailed information about catch size structure but increases memory usage.
#'
#' @return Complete ZooMSS model results object containing:
#'   \itemize{
#'     \item param: Model parameters and environmental forcing data
#'     \item time: Time values corresponding to saved results (accounting for isave)
#'     \item abundance: Abundance time series (time x groups x size classes)
#'     \item growth: Growth rate time series
#'     \item mortality: Mortality rate time series
#'     \item diet: Diet composition time series
#'     \item Additional model structure and kernel data
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default groups
#' env_data <- createEnviroData(10, 0.01)
#' input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
#' results <- zoomss_model(input_params, isave = 50)
#'
#' # Using custom groups
#' Groups <- getGroups()  # Get default groups
#' Groups$W0[1] <- -12.5          # Modify a parameter
#' results <- zoomss_model(input_params, Groups, isave = 100)
#'
#' # Loading groups from file
#' custom_groups <- getGroups(source = "file", file = "my_groups.csv")
#' results <- zoomss_model(input_params, custom_groups)
#' }
#'
zoomss_model <- function(input_params, Groups = NULL, isave = 1, fishing_params = NULL, save_catch_by_size = FALSE){

  # Handle default Groups parameter
  if (is.null(Groups)) {
    Groups <- getGroups(source = "default")
  } else {
    # Validate user-provided Groups
    validateGroups(Groups)
  }

  input_params <- untibble(input_params)

  # Validate that input_params has the required environmental data
  if (nrow(input_params) > 1 && all(c("time", "sst", "chl") %in% names(input_params))) {
    # Environmental time series found - proceed with model
  } else {
    stop("No environmental time series provided and input_params doesn't contain expanded time series data")
  }

  ################### RUN THE MODEL ###################
  param <- zoomss_params(Groups, input_params, isave, fishing_params, save_catch_by_size) # Set up parameter list
  model <- zoomss_setup(param) # Set up model equation stuff
  model_output <- zoomss_run(model) # Run the model

  my_rename <- function(.x, ..., .strict = TRUE) {
    pos <- tidyselect::eval_rename(quote(c(...)), .x, strict = .strict)
    names(.x)[pos] <- names(pos)
    .x
  }

  model_output <- model_output %>%
    my_rename(abundance = "N", growth = "gg", mortality = "Z", 
              fishing_mortality = "F_save") # Make variables more easily readable
  
  # Rename catch outputs based on save_catch_by_size setting
  if (param$save_catch_by_size) {
    model_output <- model_output %>%
      my_rename(catch_by_size = "catch_by_size_save", catch = "catch_save")
  } else {
    model_output <- model_output %>%
      my_rename(catch = "catch_save")
  }

  model_output$biomass = getBiomass(model_output, units = "ww")
  model_output$biomassC <- getBiomass(model_output, units = "carbon")

  return(model_output)
}
