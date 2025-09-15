# Comprehensive ZooMSS Reproduction Analysis
# Long-term simulation with proper visualization and comparison

library(devtools)
load_all()
library(magrittr)  # For %>% operator

cat("=== COMPREHENSIVE ZOOMSS REPRODUCTION ANALYSIS ===\n")
cat("Testing with long-term simulation (200 years) for stability\n\n")

# Create long-term environmental forcing (200 years, seasonal cycles)
sim_years <- 200
time_steps_per_year <- 12  # Monthly resolution
total_steps <- sim_years * time_steps_per_year

test_time <- seq(0, sim_years, length.out = total_steps)
test_env <- data.frame(
  time = test_time,
  sst = 15 + 3 * sin(2 * pi * test_time),  # Seasonal temperature (12-18°C)
  chl = exp(-3 + 0.5 * sin(2 * pi * test_time))  # Seasonal chlorophyll
)

cat("Environmental forcing created:\n")
cat("- Duration:", sim_years, "years\n")
cat("- Time steps:", total_steps, "\n")
cat("- SST range:", round(range(test_env$sst), 1), "°C\n")
cat("- Chl range:", round(range(test_env$chl), 4), "mg/m³\n\n")

# Load Groups and add reproduction parameters
Groups <- getGroups()
source("R/update_groups_for_reproduction.R")
Groups <- add_reproduction_params(Groups)

cat("Groups loaded with reproduction parameters\n")
fish_names <- Groups$Species[Groups$Type == "Fish"]
zoo_names <- Groups$Species[Groups$Type == "Zooplankton"]
cat("Fish groups:", paste(fish_names, collapse = ", "), "\n")
cat("Zooplankton groups:", paste(zoo_names, collapse = ", "), "\n\n")

# Create input parameters
input_params <- createInputParams(test_env$time, test_env$sst, test_env$chl)

# Setup model parameters with appropriate saving frequency
# Save every 6 months for detailed analysis (isave = 6)
cat("Setting up model parameters...\n")
params <- zoomss_params(Groups, input_params, isave = 6)

cat("Model parameters:\n")
cat("- Fish groups:", params$num_fish, "\n") 
cat("- Zooplankton groups:", params$num_zoo, "\n")
cat("- Size classes:", params$ngrid, "\n")
cat("- Save frequency: every", params$isave, "time steps\n")
cat("- Total saves:", params$nsave, "\n\n")

# Setup and run model
cat("Initializing model structure...\n")
model <- zoomss_setup(params)

cat("Starting long-term simulation...\n")
cat("This may take several minutes for", sim_years, "years...\n")
start_time <- Sys.time()

# Run the full simulation
results <- zoomss_run(model)

# Apply the same processing as zoomss_model to get biomass
# Rename arrays for compatibility with getBiomass function
my_rename <- function(.x, ..., .strict = TRUE) {
  pos <- tidyselect::eval_rename(quote(c(...)), .x, strict = .strict)
  names(.x)[pos] <- names(pos)
  .x
}

results <- results %>%
  my_rename(abundance = "N", growth = "gg", mortality = "Z")

# Calculate biomass arrays
results$biomass <- getBiomass(results, units = "ww")
results$biomassC <- getBiomass(results, units = "carbon")

end_time <- Sys.time()
cat("Simulation completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n\n")

# === ANALYSIS SECTION ===
cat("=== ANALYSIS RESULTS ===\n\n")

# Check reproduction status
if(!is.null(results$reproduction_enabled) && results$reproduction_enabled) {
  cat("✅ Fish reproduction was ENABLED\n")
  
  # Analyze SSB dynamics
  if(!is.null(results$SSB)) {
    cat("\n--- Spawning Stock Biomass (SSB) Analysis ---\n")
    SSB_final <- results$SSB[, ncol(results$SSB)]
    SSB_mean <- rowMeans(results$SSB[, (ncol(results$SSB)-10):ncol(results$SSB)], na.rm = TRUE)
    
    for(i in seq_along(fish_names)) {
      cat(sprintf("%s: Final SSB = %.2e, Mean (last 10 saves) = %.2e\n", 
                  fish_names[i], SSB_final[i], SSB_mean[i]))
    }
    
    # Check for SSB trends
    cat("\nSSB Time Series Summary:\n")
    cat("- SSB dimensions:", paste(dim(results$SSB), collapse = " x "), "\n")
    cat("- Non-zero SSB values:", sum(results$SSB > 0, na.rm = TRUE), "out of", 
        length(results$SSB), "total\n")
    
    # Calculate SSB coefficient of variation (stability)
    for(i in seq_along(fish_names)) {
      if(SSB_mean[i] > 0) {
        cv <- sd(results$SSB[i, ], na.rm = TRUE) / SSB_mean[i]
        cat(sprintf("%s CV: %.2f (%.1f%% variation)\n", fish_names[i], cv, cv*100))
      }
    }
  }
} else {
  cat("❌ Fish reproduction was DISABLED\n")
}

# Check biomass units and scales
cat("\n--- Biomass and Abundance Analysis ---\n")

# Check if biomass arrays exist in results
if("biomass" %in% names(results)) {
  final_biomass <- results$biomass
  cat("Using pre-calculated biomass from results\n")
} else if("N" %in% names(results)) {
  final_biomass <- getBiomass(results, units = "ww")
  cat("Calculating biomass from abundance\n")
} else {
  cat("⚠️  No biomass or abundance data found in results\n")
  cat("Available result components:", paste(names(results), collapse = ", "), "\n")
  final_biomass <- NULL
}

if(!is.null(final_biomass)) {
  cat("Final biomass summary (wet weight, g/m³):\n")
  total_biomass <- rowSums(final_biomass[dim(final_biomass)[1],,], na.rm = TRUE)
  for(i in seq_along(Groups$Species)) {
    cat(sprintf("%s: %.2e g/m³\n", Groups$Species[i], total_biomass[i]))
  }
  
  # Also check carbon biomass if available
  if("biomassC" %in% names(results)) {
    final_biomassC <- results$biomassC
    total_biomassC <- rowSums(final_biomassC[dim(final_biomassC)[1],,], na.rm = TRUE)
    cat("\nFinal biomass summary (carbon, g C/m³):\n") 
    for(i in seq_along(Groups$Species)) {
      cat(sprintf("%s: %.2e g C/m³\n", Groups$Species[i], total_biomassC[i]))
    }
    
    cat("\nTotal ecosystem biomass:\n")
    cat("- Wet weight:", sprintf("%.2e g/m³", sum(total_biomass, na.rm = TRUE)), "\n")
    cat("- Carbon:", sprintf("%.2e g C/m³", sum(total_biomassC, na.rm = TRUE)), "\n")
  } else {
    cat("\nTotal ecosystem biomass (wet weight):", sprintf("%.2e g/m³", sum(total_biomass, na.rm = TRUE)), "\n")
  }
  
  # Model stability analysis
  cat("\n--- Model Stability Analysis ---\n")
  n_time <- dim(final_biomass)[1]
  early_period <- 1:min(50, floor(n_time/4))  # First quarter or 50 saves
  late_period <- max(1, n_time-49):n_time     # Last 50 saves
  
  early_biomass <- apply(final_biomass[early_period,,], c(2), mean, na.rm = TRUE)
  late_biomass <- apply(final_biomass[late_period,,], c(2), mean, na.rm = TRUE)
  
  cat("Biomass change from early to late simulation:\n")
  for(i in seq_along(Groups$Species)) {
    if(early_biomass[i] > 0) {
      change <- (late_biomass[i] - early_biomass[i]) / early_biomass[i] * 100
      cat(sprintf("%s: %.1f%% change\n", Groups$Species[i], change))
    }
  }
} else {
  cat("Skipping biomass analysis due to missing data\n")
}

# Save results for plotting
cat("\n--- Saving Results for Visualization ---\n")
save(results, Groups, params, file = "long_term_reproduction_results.RData")
cat("Results saved to 'long_term_reproduction_results.RData'\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Run the visualization script next to create detailed plots\n")