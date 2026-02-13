# Step 1.1: Update GroupInputs.csv and rebuild package data
gi <- read.csv("data-raw/GroupInputs.csv", stringsAsFactors = FALSE)

# Add new energy budget columns
gi$AssimCategory <- c("Protist", "Protist", "Gelatinous", "Crustacean",
                      "Crustacean", "Crustacean", "MuscularInvert",
                      "Gelatinous", "Gelatinous", "Fish", "Fish", "Fish")
gi$Kappa <- c(0.7, 0.7, 0.6, 0.7, 0.7, 0.7, 0.7, 0.5, 0.6, NA, NA, NA)
gi$MetabConst <- 0
gi$MetabExp <- 0.75
gi$StarvSens <- 0.3

write.csv(gi, "data-raw/GroupInputs.csv", row.names = FALSE)
cat("CSV updated:", ncol(gi), "columns\n")
cat("New columns:", paste(tail(names(gi), 5), collapse = ", "), "\n")

# Rebuild package data
GroupInputs <- readr::read_csv("data-raw/GroupInputs.csv", show_col_types = FALSE)
usethis::use_data(GroupInputs, overwrite = TRUE)
cat("Package data rebuilt.\n")
