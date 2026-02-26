# ============================================================================
# PRE-COMPUTE ALL ZOOMSS MODELS FOR THE REPRODUCTION VIGNETTE
# ============================================================================
#
# Run this script ONCE from the R console (not inside knitr) before knitting
# the vignette. It parallelises all gradient models using a PSOCK cluster
# that works on Windows, Linux, and macOS.
#
# Models are saved to vignettes/cache/ as .rds files. The vignette Rmd
# loads them with readRDS() — no parallel code needed at knit time.
#
# NOTE: This script runs the full model set twice:
#       1. All fish groups with repro_eff = 0.1  (tag: "reproeff0.1")
#       2. All fish groups with repro_eff = 0.01 (tag: "reproeff0.01")
#
# Usage:
#   source("vignettes/run_gradient_models_DBPM_repro.R")
#
# To force re-run of all models, delete the cache directory first:
#   unlink("vignettes/cache", recursive = TRUE)
#
# Hardware: Tested on 16-core / 64 GB RAM Windows machine.
#           Adjust n_cores below if memory is a concern (~4 GB per worker).
# ============================================================================

library(parallel)
devtools::load_all(here::here(), quiet = TRUE)

# ---- Configuration --------------------------------------------------------
n_cores    <- max(1, detectCores() - 1)  # 15 on a 16-thread machine
cache_dir  <- here::here("vignettes", "cache")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

# ---- Helper: run a single model if not cached -----------------------------
run_single <- function(cache_name, sst, chl, Groups, cache_dir, isave = 2) {
  cache_file <- file.path(cache_dir, paste0(cache_name, ".rds"))
  if (file.exists(cache_file)) return(paste0("CACHED: ", cache_name))

  env <- createInputParams(
    time = seq(0, 1000, by = 0.1),
    sst  = sst,
    chl  = chl
  )
  mdl <- zoomss_model(input_params = env, Groups = Groups, isave = isave)
  tmp <- tempfile(tmpdir = cache_dir, fileext = ".rds")
  saveRDS(mdl, tmp)
  file.rename(tmp, cache_file)
  return(paste0("RAN: ", cache_name))
}

# ---- Repro_eff values to loop over ----------------------------------------
repro_eff_values <- c(0.1, 0.01)

for (re_val in repro_eff_values) {

  # Build a filename-safe tag: "reproeff0.1" or "reproeff0.01"
  re_tag <- paste0("reproeff", re_val)

  # Set repro_eff for all fish groups
  Groups <- getGroups()
  Groups$repro_eff[Groups$Type == "Fish"] <- re_val

  cat("==================================================\n")
  cat("ZooMSS Gradient Model Pre-Computation\n")
  cat("  Fish repro_eff set to", re_val, "\n")
  cat("==================================================\n")
  cat("Workers:   ", n_cores, "\n")
  cat("Cache dir: ", cache_dir, "\n")
  cat("Fish repro_eff values:\n")
  print(Groups[Groups$Type == "Fish", c("Species", "repro_eff")])
  cat("==================================================\n\n")

  # ---- 1. Sequential single-run models (constant, seasonal) ----------------
  cat("--- Part 1: Single-run models (sequential) ---\n")

  # Constant
  cache_file <- file.path(cache_dir, paste0("mdl_constant_", re_tag, "_1000yr.rds"))
  if (!file.exists(cache_file)) {
    cat("Running: mdl_constant_", re_tag, "_1000yr\n", sep = "")
    env_constant <- createInputParams(time = seq(0, 1000, by = 0.1), sst = 15, chl = 1.0)
    mdl <- zoomss_model(input_params = env_constant, Groups = Groups, isave = 2)
    saveRDS(mdl, cache_file)
    cat("  Saved.\n")
  } else {
    cat("Cached:  mdl_constant_", re_tag, "_1000yr\n", sep = "")
  }

  # Seasonal
  cache_file <- file.path(cache_dir, paste0("mdl_seasonal_", re_tag, "_1000yr.rds"))
  if (!file.exists(cache_file)) {
    cat("Running: mdl_seasonal_", re_tag, "_1000yr\n", sep = "")
    env_seasonal <- createEnviroData(n_years = 1000, dt = 0.1, base_sst = 15,
                                      base_chl = 1.0, seasonal = TRUE,
                                      sst_amplitude = 4, chl_amplitude = 0.5)
    env_seasonal_params <- createInputParams(time = env_seasonal$time,
                                              sst = env_seasonal$sst,
                                              chl = env_seasonal$chl)
    mdl <- zoomss_model(input_params = env_seasonal_params, Groups = Groups, isave = 2)
    saveRDS(mdl, cache_file)
    cat("  Saved.\n")
  } else {
    cat("Cached:  mdl_seasonal_", re_tag, "_1000yr\n", sep = "")
  }

  # ---- 2. Chlorophyll gradient × temperature (69 models) -------------------
  cat("\n--- Part 2: Chlorophyll gradient (69 models, parallel) ---\n")

  log10_chl_seq <- seq(-1.7, 0.5, by = 0.1)
  chl_levels    <- 10^log10_chl_seq
  sst_variants  <- c(10, 15, 20)

  task_grid <- expand.grid(
    chl_idx = seq_along(chl_levels),
    sst     = sst_variants
  )

  # Count how many actually need running
  n_to_run <- sum(!sapply(seq_len(nrow(task_grid)), function(row) {
    i <- task_grid$chl_idx[row]; temp <- task_grid$sst[row]
    file.exists(file.path(cache_dir,
      paste0(sprintf("mdl_chl_log10_%.1f_sst%d_%s_1000yr", log10(chl_levels[i]), temp, re_tag), ".rds")))
  }))
  cat("Models to run:", n_to_run, "of", nrow(task_grid), "\n")

  if (n_to_run > 0) {
    cl <- makeCluster(n_cores)
    pkg_path <- here::here()  # root of the zoomss package
    clusterExport(cl, "pkg_path")
    clusterEvalQ(cl, devtools::load_all(pkg_path, quiet = TRUE))
    clusterExport(cl, c("run_single", "task_grid", "chl_levels", "cache_dir", "Groups", "re_tag"))

    t0 <- Sys.time()
    results_chl <- parLapply(cl, seq_len(nrow(task_grid)), function(row) {
      i    <- task_grid$chl_idx[row]
      temp <- task_grid$sst[row]
      cache_name <- sprintf("mdl_chl_log10_%.1f_sst%d_%s_1000yr", log10(chl_levels[i]), temp, re_tag)
      run_single(cache_name, sst = temp, chl = chl_levels[i],
                 Groups = Groups, cache_dir = cache_dir)
    })
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    stopCluster(cl)

    cat("Chlorophyll gradient complete in", elapsed, "minutes\n")
    cat("  Ran:   ", sum(grepl("^RAN:", results_chl)), "\n")
    cat("  Cached:", sum(grepl("^CACHED:", results_chl)), "\n")
  } else {
    cat("All chlorophyll models already cached.\n")
  }

  # ---- 3. Temperature gradient (5 models) ----------------------------------
  cat("\n--- Part 3: Temperature gradient (5 models, parallel) ---\n")

  sst_levels <- c(5, 10, 15, 20, 25)

  n_to_run <- sum(!sapply(sst_levels, function(temp) {
    file.exists(file.path(cache_dir, paste0(sprintf("mdl_sst%d_chl1_%s_1000yr", temp, re_tag), ".rds")))
  }))
  cat("Models to run:", n_to_run, "of", length(sst_levels), "\n")

  if (n_to_run > 0) {
    cl <- makeCluster(n_cores)
    pkg_path <- here::here()
    clusterExport(cl, "pkg_path")
    clusterEvalQ(cl, devtools::load_all(pkg_path, quiet = TRUE))
    clusterExport(cl, c("run_single", "sst_levels", "cache_dir", "Groups", "re_tag"))

    t0 <- Sys.time()
    results_sst <- parLapply(cl, seq_along(sst_levels), function(i) {
      cache_name <- sprintf("mdl_sst%d_chl1_%s_1000yr", sst_levels[i], re_tag)
      run_single(cache_name, sst = sst_levels[i], chl = 1.0,
                 Groups = Groups, cache_dir = cache_dir)
    })
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    stopCluster(cl)

    cat("Temperature gradient complete in", elapsed, "minutes\n")
    cat("  Ran:   ", sum(grepl("^RAN:", results_sst)), "\n")
    cat("  Cached:", sum(grepl("^CACHED:", results_sst)), "\n")
  } else {
    cat("All temperature models already cached.\n")
  }

  cat("\n--- Finished repro_eff =", re_val, "---\n\n")

}  # end repro_eff loop

# ---- Summary --------------------------------------------------------------
cat("\n==================================================\n")
all_rds <- list.files(cache_dir, pattern = "\\.rds$")
cat("Done! Total cached models:", length(all_rds), "\n")
cat("Cache size:", round(sum(file.size(file.path(cache_dir, all_rds))) / 1e9, 2), "GB\n")
cat("==================================================\n")
cat("\nYou can now knit the vignette:\n")
cat('  rmarkdown::render("vignettes/reproduction_v6.Rmd")\n')
