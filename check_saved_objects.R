# Check saved objects
load("long_term_reproduction_results.RData")
cat("Available objects:\n")
print(ls())
cat("\nObject details:\n")
for(obj in ls()) {
  cat(obj, ":", class(get(obj)), "\n")
  if(is.array(get(obj)) || is.matrix(get(obj))) {
    cat("  Dimensions:", dim(get(obj)), "\n")
  }
}

# Check results structure
cat("\nResults structure:\n")
print(names(results))
cat("\nGroups structure:\n")
print(names(Groups))
print(head(Groups))

# Check if we have SSB data
if("ssb" %in% names(results)) {
  cat("\nSSB dimensions:", dim(results$ssb), "\n")
}
if("abundance" %in% names(results)) {
  cat("\nAbundance dimensions:", dim(results$abundance), "\n")
}
if("biomass" %in% names(results)) {
  cat("\nBiomass dimensions:", dim(results$biomass), "\n")
}