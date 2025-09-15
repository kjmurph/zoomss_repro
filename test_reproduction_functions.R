# Quick test script for fish reproduction
# Test reproduction functions with minimal setup

# Load the package and functions
library(devtools)
load_all()
source("R/zoomss_reproduction.R")
source("R/update_groups_for_reproduction.R")

# Get Groups with reproduction parameters
Groups <- getGroups()
Groups <- add_reproduction_params(Groups)

# Validate reproduction setup
cat("=== REPRODUCTION SETUP TEST ===\n")
cat("Fish groups found:", sum(Groups$Type == "Fish"), "\n")
cat("Reproduction enabled:", sum(Groups$Repro[Groups$Type == "Fish"] > 0), "\n")

# Test individual functions
cat("\n=== FUNCTION TESTS ===\n")

# Test maturity ogive
w <- 10^seq(-3, 7, 0.1)  # Weight vector
fish_grps <- which(Groups$Type == "Fish")
test_grp <- fish_grps[1]

cat("Testing maturity ogive...\n")
mat_ogive <- calculate_maturity_ogive(w, Groups$Wmat[test_grp])
cat("Maturity ogive range:", round(range(mat_ogive), 3), "\n")

cat("Testing reproductive investment...\n")
repro_inv <- calculate_reproductive_investment(w, Groups$ReproInvest[test_grp], 
                                             Groups$ReproExp[test_grp], Groups$Wmat[test_grp])
cat("Reproductive investment range:", round(range(repro_inv), 3), "\n")

cat("Testing recruitment parameters...\n")
# Create dummy abundance matrix
N_test <- matrix(runif(length(Groups$Species) * length(w), 0.1, 1), 
                nrow = length(Groups$Species), ncol = length(w))

recruit_params <- initialize_recruitment_params(Groups, fish_grps, N_test, w)
cat("Recruitment parameters initialized for", length(recruit_params), "fish groups\n")

cat("Testing SSB calculation...\n")
repro_energy_test <- matrix(0, nrow = length(Groups$Species), ncol = length(w))
repro_energy_test[fish_grps, ] <- 0.01  # Small reproductive energy

SSB_test <- calculate_spawning_stock_biomass(N_test, repro_energy_test, w, fish_grps)
cat("SSB calculated:", round(SSB_test, 6), "\n")

cat("Testing Beverton-Holt recruitment...\n")
recruitment_test <- calculate_beverton_holt_recruitment(SSB_test[1], recruit_params[[1]])
cat("Recruitment value:", round(recruitment_test, 6), "\n")

cat("\n=== ALL FUNCTION TESTS PASSED ===\n")