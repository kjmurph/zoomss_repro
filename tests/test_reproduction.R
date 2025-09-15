# test_reproduction.R
# Test script to verify fish reproduction implementation

library(tidyverse)

# Source new functions
source("R/zoomss_reproduction.R")
source("R/update_groups_for_reproduction.R")

#' Run reproduction test
run_reproduction_test <- function(test_years = 10, 
                                 with_reproduction = TRUE,
                                 plot_results = TRUE) {
  
  message("Starting reproduction test...")
  
  # Load and prepare Groups data
  Groups <- read_csv("data-raw/GroupInputs.csv", show_col_types = FALSE)
  
  if(with_reproduction) {
    Groups <- add_reproduction_params(Groups)
    message("Reproduction enabled")
  } else {
    Groups$Repro <- 0  # Disable reproduction
    message("Reproduction disabled (control)")
  }
  
  # Create simple test environment
  n_time <- test_years * 12  # Monthly timesteps
  test_env <- data.frame(
    time = seq(0, test_years, length.out = n_time),
    sst = 15 + 5 * sin(2 * pi * seq(0, test_years, length.out = n_time)),  # Seasonal temperature
    chla = exp(-3 + 0.5 * sin(2 * pi * seq(0, test_years, length.out = n_time)))  # Seasonal chlorophyll
  )
  
  # Set up parameters
  params <- zoomss_params(Groups, test_env)
  
  # Initialize model
  model <- zoomss_setup(params)
  
  # Run simulation
  results <- zoomss_run(model)
  
  # Extract diagnostics
  fish_grps <- which(Groups$Type == "Fish")
  
  if(with_reproduction && length(fish_grps) > 0) {
    diagnostics <- check_reproduction_diagnostics(results, fish_grps)
    
    message("\nReproduction Diagnostics:")
    message("Mean SSB by group: ", paste(round(diagnostics$SSB_mean, 2), collapse = ", "))
    message("Mean recruitment: ", paste(round(diagnostics$recruitment_mean, 6), collapse = ", "))
    message("Reproductive energy fraction: ", 
            paste(round(diagnostics$mean_repro_investment, 3), collapse = ", "))
  }
  
  return(list(
    results = results,
    params = params,
    diagnostics = if(with_reproduction) diagnostics else NULL
  ))
}

# Run quick test if script is sourced
if(interactive()) {
  message("Running quick reproduction test...")
  quick_test <- run_reproduction_test(test_years = 5)
}

