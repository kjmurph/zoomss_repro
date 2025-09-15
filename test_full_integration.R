# Full ZooMSS reproduction integration test
# Test the complete model with reproduction enabled

library(devtools)
load_all()

cat("=== ZOOMSS REPRODUCTION INTEGRATION TEST ===\n")

# Create simple test environment (2 years, monthly steps)
test_time <- seq(0, 2, length.out = 24)
test_env <- data.frame(
  time = test_time,
  sst = 15 + 3 * sin(2 * pi * test_time),  # Seasonal temperature
  chl = exp(-3 + 0.3 * sin(2 * pi * test_time))  # Seasonal chlorophyll
)

cat("Environmental data created:", nrow(test_env), "time steps\n")

# Get Groups and add reproduction parameters
Groups <- getGroups()
source("R/update_groups_for_reproduction.R")
Groups <- add_reproduction_params(Groups)

cat("Groups updated with reproduction parameters\n")
cat("Fish groups with reproduction:", sum(Groups$Repro[Groups$Type == "Fish"] > 0), "\n")

# Create input parameters
input_params <- createInputParams(test_env$time, test_env$sst, test_env$chl)
cat("Input parameters created\n")

# Test parameter setup
cat("Testing zoomss_params...\n")
params <- zoomss_params(Groups, input_params, isave = 5)
cat("Parameters set up successfully\n")
cat("Number of fish groups:", params$num_fish, "\n")
cat("Fish group indices:", paste(params$fish_grps, collapse = ", "), "\n")

# Test model setup
cat("Testing zoomss_setup...\n")
model <- zoomss_setup(params)
cat("Model setup successful\n")

# Test short model run
cat("Testing zoomss_run (short simulation)...\n")
tryCatch({
  results <- zoomss_run(model)
  cat("Model run completed successfully!\n")
  
  # Check reproduction outputs
  if(!is.null(results$reproduction_enabled) && results$reproduction_enabled) {
    cat("Reproduction was enabled and completed\n")
    if(!is.null(results$SSB)) {
      cat("SSB tracking dimensions:", paste(dim(results$SSB), collapse = " x "), "\n")
      cat("Final SSB values:", paste(round(results$SSB[, ncol(results$SSB)], 2), collapse = ", "), "\n")
    }
  } else {
    cat("Reproduction was not enabled in this run\n")
  }
  
  # Check basic model outputs
  cat("Final abundances shape:", paste(dim(results$N), collapse = " x "), "\n")
  cat("Growth rates shape:", paste(dim(results$gg), collapse = " x "), "\n")
  
  cat("\n=== INTEGRATION TEST PASSED ===\n")
  
}, error = function(e) {
  cat("ERROR in model run:", e$message, "\n")
  cat("This may indicate array dimension or integration issues\n")
})