# update_groups_for_reproduction.R
# Script to add reproduction parameters to GroupInputs

library(readr)
library(dplyr)

#' Add reproduction parameters to Groups dataframe
#' Can be called in zoomss_params() or used to update CSV
add_reproduction_params <- function(Groups) {
  
  # Check if reproduction parameters already exist
  if(!"ReproInvest" %in% names(Groups)) {
    Groups$ReproInvest <- NA_real_
  }
  
  if(!"ReproExp" %in% names(Groups)) {
    Groups$ReproExp <- NA_real_
  }
  
  # Identify fish groups
  fish_grps <- which(Groups$Type == "Fish")
  
  # Set reproductive parameters for fish
  if(length(fish_grps) > 0) {
    # Reproductive investment at maturity size
    # Higher for smaller fish (r-selected) vs larger fish (K-selected)
    # These values are based on Gunderson 1997 and similar to DBPM
    if(length(fish_grps) >= 3) {
      Groups$ReproInvest[fish_grps] <- c(0.10, 0.08, 0.06)  # Small, Medium, Large
    } else {
      Groups$ReproInvest[fish_grps] <- rep(0.08, length(fish_grps))  # Default value
    }
    
    # Size scaling exponent (negative = decreases with size above maturity)
    Groups$ReproExp[fish_grps] <- rep(-0.25, length(fish_grps))
    
    # Enable reproduction by setting Repro flag to 1
    Groups$Repro[fish_grps] <- 1
    
    message("Added fish reproduction parameters to Groups dataframe")
    message("Fish groups with reproduction: ", paste(Groups$Species[fish_grps], collapse = ", "))
  }
  
  return(Groups)
}

#' Update the CSV file with reproduction parameters
update_groups_csv <- function(csv_path = "data-raw/GroupInputs.csv", 
                             output_path = NULL) {
  
  # Read existing GroupInputs
  Groups <- read_csv(csv_path, show_col_types = FALSE)
  
  # Add reproduction parameters
  Groups <- add_reproduction_params(Groups)
  
  # Set output path
  if(is.null(output_path)) {
    output_path <- csv_path
  }
  
  # Write updated CSV
  write_csv(Groups, output_path)
  
  message("Updated CSV file: ", output_path)
  
  return(Groups)
}

#' Validate reproduction parameters
validate_reproduction_params <- function(Groups) {
  
  fish_grps <- which(Groups$Type == "Fish")
  
  if(length(fish_grps) == 0) {
    warning("No fish groups found in Groups dataframe")
    return(FALSE)
  }
  
  # Check all fish have reproduction parameters
  if(any(is.na(Groups$ReproInvest[fish_grps]))) {
    warning("Some fish groups missing ReproInvest values")
  }
  
  if(any(is.na(Groups$ReproExp[fish_grps]))) {
    warning("Some fish groups missing ReproExp values")
  }
  
  # Check parameter ranges are sensible
  if(any(Groups$ReproInvest[fish_grps] < 0 | Groups$ReproInvest[fish_grps] > 0.5)) {
    warning("ReproInvest should be between 0 and 0.5")
  }
  
  if(any(Groups$ReproExp[fish_grps] > 0)) {
    warning("ReproExp should be negative (reproduction decreases with size)")
  }
  
  # Check maturity sizes are set
  if(any(is.na(Groups$Wmat[fish_grps]))) {
    stop("Fish groups must have Wmat (maturity size) defined")
  }
  
  # Check egg sizes (W0) are set
  if(any(is.na(Groups$W0[fish_grps]))) {
    stop("Fish groups must have W0 (egg size) defined")
  }
  
  message("Reproduction parameters validated successfully")
  
  return(TRUE)
}

