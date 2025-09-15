# Quick diagnostic to check results structure
library(devtools)
load_all()

# Run a short test to see result structure
Groups <- getGroups()
source("R/update_groups_for_reproduction.R")
Groups <- add_reproduction_params(Groups)

# Short test environment
test_env <- data.frame(
  time = seq(0, 2, length.out = 24),
  sst = rep(15, 24),
  chl = rep(0.05, 24)
)

input_params <- createInputParams(test_env$time, test_env$sst, test_env$chl)
params <- zoomss_params(Groups, input_params, isave = 5)
model <- zoomss_setup(params)
results <- zoomss_run(model)

cat("=== RESULTS OBJECT STRUCTURE ===\n")
cat("Names in results object:\n")
print(names(results))

cat("\nDimensions of key arrays:\n")
for(name in names(results)) {
  if(is.array(results[[name]]) || is.matrix(results[[name]])) {
    cat(sprintf("%s: %s\n", name, paste(dim(results[[name]]), collapse = " x ")))
  }
}

cat("\n=== TESTING BIOMASS CALCULATION ===\n")
# Try different biomass calculations
tryCatch({
  biomass_ww <- getBiomass(results, units = "ww")
  cat("✅ getBiomass with 'ww' works\n")
  cat("Biomass dimensions:", paste(dim(biomass_ww), collapse = " x "), "\n")
}, error = function(e) {
  cat("❌ getBiomass with 'ww' failed:", e$message, "\n")
})

tryCatch({
  biomass_c <- getBiomass(results, units = "carbon") 
  cat("✅ getBiomass with 'carbon' works\n")
}, error = function(e) {
  cat("❌ getBiomass with 'carbon' failed:", e$message, "\n")
})

# Check if biomass is pre-calculated
if("biomass" %in% names(results)) {
  cat("✅ Pre-calculated 'biomass' found in results\n")
}
if("biomassC" %in% names(results)) {
  cat("✅ Pre-calculated 'biomassC' found in results\n")
}