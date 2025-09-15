# Fix ZooMSS Reproduction Parameters
# Adjust parameters to achieve more realistic reproduction levels

library(devtools)
load_all()
source("R/zoomss_reproduction.R")
source("R/update_groups_for_reproduction.R")

cat("=== FIXING REPRODUCTION PARAMETERS ===\n")

# Function to update reproduction parameters with more realistic values
update_reproduction_params_realistic <- function(Groups) {
  
  fish_grps <- which(Groups$Type == "Fish")
  
  if(length(fish_grps) > 0) {
    cat("Updating reproduction parameters for more realistic reproduction...\n")
    
    # Increase reproductive investment (original was very conservative)
    # Make smaller fish invest more heavily in reproduction (r-strategy)
    Groups$ReproInvest[fish_grps] <- c(0.25, 0.20, 0.15)  # Increased from 0.10, 0.08, 0.06
    
    # Reduce the negative size scaling to maintain reasonable investment across sizes
    Groups$ReproExp[fish_grps] <- rep(-0.15, length(fish_grps))  # Less negative than -0.25
    
    # Consider making maturity sizes smaller for better reproduction
    # Move Large fish maturity down slightly
    Groups$Wmat[fish_grps[3]] <- 3  # From 4 to 3 (10^3 = 1000g instead of 10000g)
    
    cat("Updated parameters:\n")
    for(i in seq_along(fish_grps)) {
      grp <- fish_grps[i]
      cat("  ", Groups$Species[grp], ": ReproInvest =", Groups$ReproInvest[grp], 
          ", ReproExp =", Groups$ReproExp[grp], ", Wmat =", Groups$Wmat[grp], "\n")
    }
  }
  
  return(Groups)
}

# Test the updated parameters
Groups <- getGroups()
Groups <- add_reproduction_params(Groups)
Groups <- update_reproduction_params_realistic(Groups)

# Save updated Groups to CSV
write.csv(Groups, "data-raw/GroupInputs_realistic_repro.csv", row.names = FALSE)
cat("Saved updated parameters to: data-raw/GroupInputs_realistic_repro.csv\n")

# Test reproduction with new parameters
cat("\n=== TESTING UPDATED PARAMETERS ===\n")

# Create test environment
test_years <- 3
test_time <- seq(0, test_years, length.out = test_years * 12)
test_env <- data.frame(
  time = test_time,
  sst = 15 + 3 * sin(2 * pi * test_time),
  chl = exp(-3 + 0.3 * sin(2 * pi * test_time))
)

input_params <- createInputParams(test_env$time, test_env$sst, test_env$chl)
params <- zoomss_params(Groups, input_params, isave = 4)

# Run model with updated parameters
cat("Running model with updated reproduction parameters...\n")
model <- zoomss_setup(params)
results <- zoomss_run(model)

if(!is.null(results$SSB)) {
  cat("Updated SSB results:\n")
  fish_grps <- which(Groups$Type == "Fish")
  
  for(i in seq_along(fish_grps)) {
    grp_name <- Groups$Species[fish_grps[i]]
    SSB_timeseries <- results$SSB[i, ]
    cat("  ", grp_name, ": Mean SSB =", round(mean(SSB_timeseries, na.rm = TRUE), 4), 
        ", Final SSB =", round(SSB_timeseries[length(SSB_timeseries)], 4), "\n")
  }
  
  # Check if SSB values are more reasonable
  total_SSB <- sum(results$SSB[, ncol(results$SSB)])
  cat("Total final SSB:", round(total_SSB, 4), "\n")
  
  if(total_SSB > 0.01) {
    cat("SUCCESS: SSB levels are now realistic!\n")
  } else if(total_SSB > 1e-4) {
    cat("IMPROVED: SSB levels are better but could be higher\n")
  } else {
    cat("STILL LOW: SSB levels need further adjustment\n")
  }
} else {
  cat("ERROR: SSB not found in results\n")
}

# Test recruitment function with new SSB levels
if(!is.null(results$SSB)) {
  cat("\n=== TESTING RECRUITMENT FUNCTION ===\n")
  
  fish_grps <- which(Groups$Type == "Fish")
  recruit_params <- initialize_recruitment_params(Groups, fish_grps, model$N[1,,], params$w)
  
  final_SSB <- results$SSB[, ncol(results$SSB)]
  
  for(i in seq_along(fish_grps)) {
    grp_name <- Groups$Species[fish_grps[i]]
    
    # Test Beverton-Holt recruitment
    recruitment_BH <- calculate_beverton_holt_recruitment(final_SSB[i], recruit_params[[i]])
    
    # Test Ricker recruitment  
    recruitment_R <- calculate_ricker_recruitment(final_SSB[i], recruit_params[[i]])
    
    cat("  ", grp_name, ":\n")
    cat("    SSB:", round(final_SSB[i], 4), "\n")
    cat("    Beverton-Holt recruitment:", round(recruitment_BH, 4), "\n")
    cat("    Ricker recruitment:", round(recruitment_R, 4), "\n")
    cat("    Virgin recruitment (R0):", round(recruit_params[[i]]$R0, 4), "\n")
  }
}

cat("\n=== PARAMETER ADJUSTMENT COMPLETE ===\n")
cat("Use the updated CSV file for more realistic reproduction\n")