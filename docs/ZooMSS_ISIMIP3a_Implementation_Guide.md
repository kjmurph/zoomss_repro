# ZooMSS-ISIMIP3a Implementation Guide (Phase B)

**Purpose:** Guide for implementing the `zoomss-isimip3a` repository following the ISIMIP3a protocol  
**Date:** January 2026  
**Dependencies:** Requires completed Phase A (zoomss package with fishing implementation)

---

## Executive Summary

This guide provides step-by-step instructions for creating the `zoomss-isimip3a` repository, which implements the ISIMIP3a protocol using the `zoomss` R package as the underlying ecological model. The workflow is based on the DBPM ISIMIP3a implementation but adapted for ZooMSS's functional group structure.

**Reference Repository:** https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a  
**Branch:** `new_features`

---

## Architecture

### Two-Repository Design

1. **`zoomss`** (existing, Phase A complete)
   - Generic size spectrum model package
   - Dynamic fishing mortality functionality
   - Standalone and reusable
   
2. **`zoomss-isimip3a`** (to be created, Phase B)
   - ISIMIP3a protocol implementation
   - Depends on `zoomss` package
   - Contains workflow scripts, data processing, and calibration

### Key Principle

**Do NOT add protocol-specific code to the `zoomss` package.** All ISIMIP3a-specific functionality (FAO regions, calibration routines, gridded output formatting) belongs in `zoomss-isimip3a`.

---

## Repository Structure

Create the following structure for `zoomss-isimip3a`:

```
zoomss-isimip3a/
├── README.md                    # Overview and quick start
├── DESCRIPTION                  # R package metadata (Depends: zoomss)
├── .gitignore
├── scripts/
│   ├── 01_climate_inputs.R      # Task B.1
│   ├── 02_regional_inputs.R     # Task B.2
│   ├── 03_fishing_inputs.R      # Task B.3
│   ├── 04_calibration.R         # Task B.4 (CRITICAL)
│   ├── 05_setup_grid.R          # Task B.5
│   ├── 06_run_model.R           # Task B.6
│   ├── 07_process_outputs.R     # Task B.7
│   └── utils/
│       ├── download_helpers.R   # Data download utilities
│       ├── netcdf_helpers.R     # NetCDF I/O
│       └── spatial_helpers.R    # Grid/region operations
├── data/
│   ├── climate/                 # ISIMIP climate forcing (gitignored)
│   ├── regions/                 # FAO/LME shapefiles
│   ├── fishing/                 # Watson catch data, effort
│   └── calibration/             # Calibrated parameters
├── outputs/
│   ├── calibration/             # Calibration diagnostics
│   ├── gridded/                 # Gridded model outputs (gitignored)
│   └── regional/                # Regional aggregations
├── tests/
│   └── testthat/
│       ├── test-climate.R
│       ├── test-calibration.R
│       └── test-outputs.R
├── docs/
│   ├── calibration_report.Rmd  # Calibration validation
│   └── protocol_compliance.md  # ISIMIP3a checklist
└── vignettes/
    └── zoomss_isimip3a_workflow.Rmd
```

---

## Phase B Tasks

### Task B.1: Climate Input Processing

**Objective:** Download and process ISIMIP climate forcing data for ZooMSS

**DBPM Reference:** Fetch and study `scripts/01_processing_dbpm_global_inputs.ipynb`

**Implementation Steps:**

1. **Identify ISIMIP Climate Data Sources**
   - ISIMIP3a Climate Input Data: https://data.isimip.org/
   - Required variables:
     - `tos` (Sea Surface Temperature, °C)
     - `chl` or `chlos` (Surface Chlorophyll, mg/m³)
     - Alternative: `phyc` (Phytoplankton Carbon) if chlorophyll unavailable — can be converted to chlorophyll values
   - Climate models: GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL
   - Scenarios: historical (1850-2014), ssp126, ssp370, ssp585 (2015-2100)

2. **Download Climate Data**
   ```r
   # scripts/01_climate_inputs.R
   
   library(ncdf4)
   library(terra)
   
   download_isimip_climate <- function(model, scenario, variable, output_dir) {
     # ISIMIP data structure:
     # https://data.isimip.org/datasets/[model]/[scenario]/[variable]/
     
     url <- sprintf(
       "https://data.isimip.org/api/v1/datasets/%s_%s_%s_global_monthly_2015_2100.nc4",
       model, scenario, variable
     )
     
     dest_file <- file.path(output_dir, basename(url))
     
     if (!file.exists(dest_file)) {
       download.file(url, dest_file, mode = "wb")
       message("Downloaded: ", dest_file)
     } else {
       message("Already exists: ", dest_file)
     }
     
     return(dest_file)
   }
   ```

3. **Process Climate Data for ZooMSS**
   ```r
   process_climate_for_zoomss <- function(nc_file, spatial_res = 1.0) {
     # Read NetCDF
     nc <- nc_open(nc_file)
     
     # Extract SST and chlorophyll
     sst <- ncvar_get(nc, "tos")
     chl <- ncvar_get(nc, "chl")
     lon <- ncvar_get(nc, "lon")
     lat <- ncvar_get(nc, "lat")
     time <- ncvar_get(nc, "time")
     
     nc_close(nc)
     
     # Regrid to target resolution if needed
     if (spatial_res != attr(lon, "resolution")) {
       # Use terra::resample() or similar
     }
     
     # Create data structure
     climate_data <- list(
       sst = sst,      # dims: lon × lat × time
       chl = chl,
       lon = lon,
       lat = lat,
       time = time,
       units = list(sst = "degC", chl = "mg/m3")
     )
     
     return(climate_data)
   }
   ```

4. **Validation**
   - Check spatial coverage (global ocean)
   - Verify temporal continuity
   - Compare to observations (e.g., WOA for SST)
   - Plot climatology maps

**Outputs:**
- `data/climate/[model]_[scenario]_processed.rds`
- Validation plots in `outputs/climate/`

**Tests:**
```r
# tests/testthat/test-climate.R
test_that("Climate data has correct dimensions", {
  climate <- readRDS("data/climate/GFDL-ESM4_historical_processed.rds")
  expect_true(all(c("sst", "chl", "lon", "lat", "time") %in% names(climate)))
  expect_equal(dim(climate$sst)[3], length(climate$time))
})

test_that("SST values are physically reasonable", {
  climate <- readRDS("data/climate/GFDL-ESM4_historical_processed.rds")
  expect_true(all(climate$sst >= -2, na.rm = TRUE))  # Freezing point
  expect_true(all(climate$sst <= 35, na.rm = TRUE))  # Max ocean temp
})
```

---

### Task B.2: Regional Input Processing

**Objective:** Define spatial regions (FAO, LME) and extract regional environmental conditions

**DBPM Reference:** Fetch and study `scripts/02_processing_dbpm_regional_inputs.py`

**Implementation Steps:**

1. **Download Region Shapefiles**
   ```r
   # scripts/02_regional_inputs.R
   
   library(sf)
   library(rnaturalearth)
   
   download_fao_regions <- function(output_dir) {
     # FAO Major Fishing Areas
     url <- "https://www.fao.org/fishery/geoserver/fifao/ows?service=WFS&request=GetFeature&typeName=fifao:FAO_AREAS&outputFormat=SHAPE-ZIP"
     
     dest <- file.path(output_dir, "fao_areas.zip")
     download.file(url, dest, mode = "wb")
     unzip(dest, exdir = output_dir)
     
     fao <- st_read(file.path(output_dir, "FAO_AREAS.shp"))
     return(fao)
   }
   
   download_lme_regions <- function(output_dir) {
     # Large Marine Ecosystems
     url <- "https://marineregions.org/download_file.php?name=lme_v1.zip"
     
     dest <- file.path(output_dir, "lme.zip")
     download.file(url, dest, mode = "wb")
     unzip(dest, exdir = output_dir)
     
     lme <- st_read(file.path(output_dir, "lme_v1.shp"))
     return(lme)
   }
   ```

2. **Assign Grid Cells to Regions**
   ```r
   assign_cells_to_regions <- function(grid_coords, regions_sf) {
     # Convert grid to spatial points
     grid_points <- st_as_sf(
       grid_coords,
       coords = c("lon", "lat"),
       crs = st_crs(regions_sf)
     )
     
     # Spatial join
     cell_regions <- st_join(grid_points, regions_sf)
     
     return(cell_regions)
   }
   ```

3. **Calculate Regional Climatology**
   ```r
   calculate_regional_means <- function(climate_data, cell_regions) {
     regional_means <- lapply(unique(cell_regions$region_id), function(rid) {
       # Get cells in this region
       cells <- which(cell_regions$region_id == rid)
       
       # Calculate spatial mean for each time step
       sst_mean <- apply(climate_data$sst[cells, , ], 2, mean, na.rm = TRUE)
       chl_mean <- apply(climate_data$chl[cells, , ], 2, mean, na.rm = TRUE)
       
       list(
         region_id = rid,
         sst = sst_mean,
         chl = chl_mean,
         n_cells = length(cells)
       )
     })
     
     return(regional_means)
   }
   ```

4. **Create Region-Specific Groups (Optional)**
   ```r
   # If functional group parameters vary by region
   create_regional_groups <- function(base_groups, region_id, regional_params) {
     Groups <- base_groups
     
     # Modify parameters based on region
     # Example: Different growth rates in tropical vs temperate
     if (regional_params$latitude > 30) {
       # Temperate: slower growth
       Groups$lmax <- Groups$lmax * 0.8
     }
     
     return(Groups)
   }
   ```

**Outputs:**
- `data/regions/fao_areas.rds`
- `data/regions/lme.rds`
- `data/regions/cell_region_assignments.rds`
- `data/regions/regional_climatology.rds`

**Tests:**
```r
test_that("All grid cells assigned to regions", {
  assignments <- readRDS("data/regions/cell_region_assignments.rds")
  expect_equal(sum(is.na(assignments$region_id)), 0)
})
```

---

### Task B.3: Fishing Effort Processing

**Objective:** Process fishing effort time series from Watson catch data or ISIMIP fishing input

**DBPM Reference:** Fetch and study `scripts/03_processing_effort_fishing_inputs.R`

**Implementation Steps:**

1. **Download Watson Catch Data**
   ```r
   # scripts/03_fishing_inputs.R
   
   library(data.table)
   
   download_watson_data <- function(output_dir) {
     # Watson Industrial Catch Reconstruction
     # http://www.seaaroundus.org/data/#/spatial-catch
     
     url <- "http://www.seaaroundus.org/data/watson/global_catch_v4.rds"
     dest <- file.path(output_dir, "watson_catch.rds")
     
     if (!file.exists(dest)) {
       download.file(url, dest, mode = "wb")
     }
     
     catch_data <- readRDS(dest)
     return(catch_data)
   }
   ```

2. **Extract Catch by Region and Taxonomic Group**
   ```r
   process_watson_catch <- function(watson_data, regions) {
     # Aggregate catch by:
     # - FAO region
     # - Taxonomic group (map to ZooMSS functional groups)
     # - Year
     
     catch_regional <- watson_data %>%
       filter(!is.na(fao_area)) %>%
       mutate(
         zoomss_group = map_taxonomy_to_zoomss(taxonomic_group)
       ) %>%
       group_by(fao_area, zoomss_group, year) %>%
       summarise(
         catch_tonnes = sum(catch_tonnes, na.rm = TRUE),
         .groups = "drop"
       )
     
     return(catch_regional)
   }
   
   # Mapping function: taxonomy → ZooMSS groups
   map_taxonomy_to_zoomss <- function(taxonomy) {
     # ZooMSS groups: 1-9 zooplankton, 10-12 fish
     # Watson groups: demersal fish, pelagic fish, cephalopods, etc.
     
     case_when(
       str_detect(taxonomy, "demersal|groundfish") ~ 10,
       str_detect(taxonomy, "pelagic|small pelagic") ~ 11,
       str_detect(taxonomy, "large pelagic|tuna") ~ 12,
       TRUE ~ NA_integer_
     )
   }
   ```

3. **Convert Catch to Effort or Use ISIMIP Effort Directly**
   ```r
   # Option A: Back-calculate effort from catch (if catchability known)
   catch_to_effort <- function(catch, catchability, biomass) {
     # Simple model: Catch = q × Effort × Biomass
     # Solve for Effort = Catch / (q × Biomass)
     
     # Use initial biomass estimate or assume equilibrium
     effort <- catch / (catchability * biomass)
     return(effort)
   }
   
   # Option B: Use ISIMIP fishing effort directly
   download_isimip_fishing <- function(scenario, output_dir) {
     # ISIMIP3a Fishing Input Data
     url <- sprintf(
       "https://data.isimip.org/datasets/fishing_effort_%s_global_monthly.nc4",
       scenario
     )
     
     dest <- file.path(output_dir, basename(url))
     download.file(url, dest, mode = "wb")
     
     return(dest)
   }
   ```

4. **Create Effort Time Series for ZooMSS**
   ```r
   create_effort_timeseries <- function(catch_data, start_year, end_year, dt = 1/12) {
     # Generate monthly effort from annual catch data
     
     years <- seq(start_year, end_year, by = 1)
     n_steps <- length(years) * 12
     
     # Expand annual to monthly (assuming uniform within year)
     effort_monthly <- rep(catch_data$effort_annual, each = 12)
     
     # Or interpolate for smooth transitions
     effort_smooth <- approx(
       x = years,
       y = catch_data$effort_annual,
       xout = seq(start_year, end_year, by = dt)
     )$y
     
     return(effort_smooth)
   }
   ```

**Outputs:**
- `data/fishing/watson_catch_processed.rds`
- `data/fishing/effort_timeseries_by_region.rds`
- `data/fishing/taxonomic_mapping.csv`

**Tests:**
```r
test_that("Effort time series matches model time vector", {
  effort <- readRDS("data/fishing/effort_timeseries_by_region.rds")
  env_data <- readRDS("data/climate/example_cell_env.rds")
  
  expect_equal(length(effort$region_01), nrow(env_data))
})
```

---

### Task B.4: Calibration (CRITICAL)

**Objective:** Calibrate catchability parameters to match historical catch observations

**DBPM Reference:** Fetch and study `scripts/04_calculating_dbpm_fishing_params.R` — **Study this script carefully**

**Implementation Steps:**

1. **Define Calibration Objective Function**
   ```r
   # scripts/04_calibration.R
   
   library(zoomss)
   library(optimx)
   
   calibration_objective <- function(q_params, 
                                     observed_catch, 
                                     env_data, 
                                     effort,
                                     Groups,
                                     selectivity) {
     # q_params: catchability values to optimize (length = n_groups or 1 for all fish)
     
     # Create fishing parameters
     fishing_params <- list(
       effort = effort,
       catchability = q_params,  # Being optimized
       selectivity = selectivity  # Fixed
     )
     
     # Run ZooMSS model
     mdl <- zoomss_model(
       env_data, 
       Groups, 
       isave = 12,  # Annual saves
       fishing_params = fishing_params
     )
     
     # Extract modeled catch (sum across fish groups)
     modeled_catch <- rowSums(mdl$catch[, 10:12])  # Fish groups only
     
     # Calculate residual sum of squares
     # Match temporal overlap with observed data
     n_years <- min(length(observed_catch), length(modeled_catch))
     residual <- sum((modeled_catch[1:n_years] - observed_catch[1:n_years])^2)
     
     # Add penalty for unrealistic q values
     penalty <- sum(q_params[q_params < 0 | q_params > 1]^2) * 1e6
     
     return(residual + penalty)
   }
   ```

2. **Calibration Workflow**
   ```r
   calibrate_region <- function(region_id, 
                                observed_catch_data,
                                climate_data,
                                effort_data,
                                initial_q = 0.01) {
     
     # Prepare data for this region
     region_env <- extract_regional_env(climate_data, region_id)
     region_effort <- extract_regional_effort(effort_data, region_id)
     region_catch <- observed_catch_data %>% filter(fao_area == region_id)
     
     # Create environmental input
     env_data <- createInputParams(
       time = seq(1950, 2014, by = 1/12),  # Historical period
       sst = region_env$sst,
       chl = region_env$chl
     )
     
     # Define selectivity (fixed during calibration)
     Groups <- getGroups()
     ngrps <- nrow(Groups)
     
     selectivity <- lapply(1:ngrps, function(g) {
       if (g %in% 10:12) {  # Fish groups
         list(
           type = "logistic",
           params = list(
             L50 = (Groups$W0[g] + Groups$Wmax[g]) / 2,
             L95 = Groups$W0[g] + 0.75 * (Groups$Wmax[g] - Groups$W0[g])
           ),
           log_scale = TRUE
         )
       } else {  # Non-fished groups
         list(
           type = "knife_edge",
           params = list(w_min = Groups$Wmax[g] + 1),  # Never selected
           log_scale = TRUE
         )
       }
     })
     
     # Initial parameters: one q for each fish group
     initial_params <- rep(initial_q, 3)  # Groups 10, 11, 12
     
     # Expand to full ngrps (only fish groups have non-zero q)
     q_full <- rep(0, ngrps)
     q_full[10:12] <- initial_params
     
     # Optimization
     cat("Calibrating region:", region_id, "\n")
     
     opt_result <- optim(
       par = initial_params,
       fn = function(q_fish) {
         q_all <- rep(0, ngrps)
         q_all[10:12] <- q_fish
         calibration_objective(
           q_params = q_all,
           observed_catch = region_catch$catch_tonnes,
           env_data = env_data,
           effort = region_effort,
           Groups = Groups,
           selectivity = selectivity
         )
       },
       method = "L-BFGS-B",
       lower = rep(0, 3),
       upper = rep(1, 3),
       control = list(trace = 1, maxit = 100)
     )
     
     # Store results
     calibration_result <- list(
       region_id = region_id,
       q_optimal = opt_result$par,
       convergence = opt_result$convergence,
       final_objective = opt_result$value,
       groups = 10:12,
       selectivity = selectivity
     )
     
     return(calibration_result)
   }
   ```

3. **Run Calibration for All Regions**
   ```r
   calibrate_all_regions <- function(regions, ...) {
     calibration_results <- lapply(regions, function(rid) {
       tryCatch({
         calibrate_region(rid, ...)
       }, error = function(e) {
         warning("Calibration failed for region ", rid, ": ", e$message)
         return(NULL)
       })
     })
     
     # Remove failed calibrations
     calibration_results <- calibration_results[!sapply(calibration_results, is.null)]
     
     # Save results
     saveRDS(calibration_results, "data/calibration/regional_calibration.rds")
     
     return(calibration_results)
   }
   ```

4. **Validation Plots**
   ```r
   validate_calibration <- function(calibration_result, observed_catch, env_data, effort) {
     # Re-run model with calibrated parameters
     q_full <- rep(0, nrow(Groups))
     q_full[calibration_result$groups] <- calibration_result$q_optimal
     
     fishing_params <- list(
       effort = effort,
       catchability = q_full,
       selectivity = calibration_result$selectivity
     )
     
     mdl <- zoomss_model(env_data, Groups, isave = 12, fishing_params = fishing_params)
     modeled_catch <- rowSums(mdl$catch[, calibration_result$groups])
     
     # Time series comparison plot
     plot_df <- data.frame(
       year = seq_along(modeled_catch),
       modeled = modeled_catch,
       observed = observed_catch[1:length(modeled_catch)]
     )
     
     p <- ggplot(plot_df, aes(x = year)) +
       geom_line(aes(y = modeled, color = "Modeled")) +
       geom_point(aes(y = observed, color = "Observed")) +
       labs(title = paste("Calibration Validation - Region", calibration_result$region_id),
            y = "Catch (tonnes)", x = "Year") +
       theme_minimal()
     
     # Goodness of fit metrics
     rmse <- sqrt(mean((plot_df$modeled - plot_df$observed)^2, na.rm = TRUE))
     r2 <- cor(plot_df$modeled, plot_df$observed, use = "complete.obs")^2
     
     cat("RMSE:", rmse, "\n")
     cat("R²:", r2, "\n")
     
     return(list(plot = p, rmse = rmse, r2 = r2))
   }
   ```

**Key Decision Points:**

1. **Catchability Structure:**
   - Option A: Single q for all fish groups (simple, 1 parameter per region)
   - Option B: Group-specific q (realistic, 3 parameters per region)
   - Option C: Size-dependent q (complex, but more mechanistic)
   
   **Recommendation:** Start with Option A, validate with Option B

2. **Calibration Period:**
   - Use historical period with good catch data (e.g., 1980-2014)
   - Longer period = more data but potential non-stationarity

3. **Objective Function:**
   - Sum of squared errors (simple)
   - Log-transformed errors (if catch spans orders of magnitude)
   - Weighted errors (emphasize recent data or high-catch periods)

**Outputs:**
- `data/calibration/regional_calibration.rds`
- `outputs/calibration/validation_plots/region_[ID].pdf`
- `outputs/calibration/calibration_summary.csv`

**Tests:**
```r
test_that("Calibrated catchability values are reasonable", {
  cal <- readRDS("data/calibration/regional_calibration.rds")
  
  for (region_cal in cal) {
    expect_true(all(region_cal$q_optimal >= 0))
    expect_true(all(region_cal$q_optimal <= 1))
    expect_true(region_cal$convergence == 0)  # Successful convergence
  }
})

test_that("Calibration improves fit compared to default", {
  # Run with default q = 0.01
  # Run with calibrated q
  # Compare RMSE
  expect_lt(rmse_calibrated, rmse_default)
})
```

---

### Task B.5: Gridded Model Setup

**Objective:** Set up spatial grid and prepare for parallel execution

**DBPM Reference:** Fetch and study `scripts/05_setup_gridded_DBPM.py`

**Implementation Steps:**

1. **Create Spatial Grid**
   ```r
   # scripts/05_setup_grid.R
   
   library(terra)
   library(ncdf4)
   
   create_spatial_grid <- function(resolution = 1.0, mask_land = TRUE) {
     # Global grid
     lon <- seq(-180 + resolution/2, 180 - resolution/2, by = resolution)
     lat <- seq(-90 + resolution/2, 90 - resolution/2, by = resolution)
     
     grid <- expand.grid(lon = lon, lat = lat)
     
     # Mask land cells
     if (mask_land) {
       ocean_mask <- get_ocean_mask(resolution)
       grid$ocean <- ocean_mask[cbind(grid$lon, grid$lat)]
       grid <- grid[grid$ocean == 1, ]
     }
     
     # Assign region IDs
     regions <- st_read("data/regions/fao_areas.shp")
     grid$region_id <- assign_grid_to_regions(grid, regions)
     
     return(grid)
   }
   
   get_ocean_mask <- function(resolution) {
     # Use global bathymetry (e.g., ETOPO1)
     # or climate model land-sea mask
     
     bathy_file <- "data/bathymetry/etopo1.nc"
     nc <- nc_open(bathy_file)
     depth <- ncvar_get(nc, "elevation")
     nc_close(nc)
     
     # Ocean = depth < 0
     ocean_mask <- (depth < 0) * 1
     
     return(ocean_mask)
   }
   ```

2. **Initialize NetCDF Output File**
   ```r
   initialize_output_netcdf <- function(grid, time, variables, output_file) {
     # Dimensions
     londim <- ncdim_def("lon", "degrees_east", unique(grid$lon))
     latdim <- ncdim_def("lat", "degrees_north", unique(grid$lat))
     timedim <- ncdim_def("time", "years since 1850-01-01", time)
     
     # Variables to save
     # Following ISIMIP protocol variable names
     
     tcb_var <- ncvar_def(
       name = "tcb",
       units = "g m-2",
       dim = list(londim, latdim, timedim),
       longname = "Total Consumer Biomass",
       missval = -9999
     )
     
     tc_var <- ncvar_def(
       name = "tc",
       units = "g m-2 yr-1",
       dim = list(londim, latdim, timedim),
       longname = "Total Catch",
       missval = -9999
     )
     
     # Create NetCDF file
     nc <- nc_create(output_file, list(tcb_var, tc_var))
     nc_close(nc)
     
     message("Initialized NetCDF output: ", output_file)
   }
   ```

3. **Setup Parallel Execution Framework**
   ```r
   setup_parallel_execution <- function(n_cores = NULL) {
     library(parallel)
     library(foreach)
     library(doParallel)
     
     if (is.null(n_cores)) {
       n_cores <- detectCores() - 1
     }
     
     cl <- makeCluster(n_cores)
     registerDoParallel(cl)
     
     # Export required packages to workers
     clusterEvalQ(cl, {
       library(zoomss)
       library(ncdf4)
     })
     
     return(cl)
   }
   ```

**Outputs:**
- `data/grid/spatial_grid_1deg.rds`
- `outputs/gridded/template_output.nc`

---

### Task B.6: Model Execution

**Objective:** Run ZooMSS for each grid cell with appropriate forcing and fishing

**DBPM Reference:** Fetch and study `scripts/06_running_gridded_DBPM.py`

**Implementation Steps:**

1. **Define Grid Cell Runner Function**
   ```r
   # scripts/06_run_model.R
   
   run_grid_cell <- function(cell_idx, grid, climate_data, effort_data, calibration) {
     # Extract cell coordinates
     lon <- grid$lon[cell_idx]
     lat <- grid$lat[cell_idx]
     region_id <- grid$region_id[cell_idx]
     
     # Extract environmental forcing for this cell
     cell_env <- extract_cell_climate(climate_data, lon, lat)
     
     # Create input parameters
     env_data <- createInputParams(
       time = climate_data$time,
       sst = cell_env$sst,
       chl = cell_env$chl
     )
     
     # Get calibrated fishing parameters for this region
     region_cal <- calibration[[as.character(region_id)]]
     
     if (is.null(region_cal)) {
       warning("No calibration for region ", region_id, " - using default")
       q_fish <- rep(0.01, 3)
     } else {
       q_fish <- region_cal$q_optimal
     }
     
     # Setup fishing parameters
     Groups <- getGroups()
     ngrps <- nrow(Groups)
     
     q_full <- rep(0, ngrps)
     q_full[10:12] <- q_fish
     
     fishing_params <- list(
       effort = effort_data[[as.character(region_id)]],
       catchability = q_full,
       selectivity = region_cal$selectivity
     )
     
     # Run ZooMSS
     tryCatch({
       mdl <- zoomss_model(
         env_data, 
         Groups, 
         isave = 12,  # Annual output
         fishing_params = fishing_params,
         save_catch_by_size = FALSE  # Memory efficient
       )
       
       # Extract key outputs
       result <- list(
         cell_idx = cell_idx,
         lon = lon,
         lat = lat,
         biomass = rowSums(mdl$biomass[, 10:12, ]),  # Total fish biomass
         catch = rowSums(mdl$catch[, 10:12]),        # Total catch
         success = TRUE
       )
       
     }, error = function(e) {
       warning("Failed to run cell ", cell_idx, ": ", e$message)
       result <- list(
         cell_idx = cell_idx,
         success = FALSE,
         error = e$message
       )
     })
     
     return(result)
   }
   ```

2. **Run All Grid Cells in Parallel**
   ```r
   run_gridded_model <- function(grid, climate_data, effort_data, calibration, 
                                  output_file, n_cores = NULL) {
     
     # Setup parallel backend
     cl <- setup_parallel_execution(n_cores)
     
     # Progress tracking
     n_cells <- nrow(grid)
     cat("Running ZooMSS for", n_cells, "grid cells...\n")
     
     # Parallel execution
     results <- foreach(
       i = 1:n_cells,
       .packages = c("zoomss"),
       .combine = rbind,
       .errorhandling = "pass"
     ) %dopar% {
       run_grid_cell(i, grid, climate_data, effort_data, calibration)
     }
     
     # Stop cluster
     stopCluster(cl)
     
     # Write to NetCDF
     write_results_to_netcdf(results, grid, output_file)
     
     return(results)
   }
   ```

3. **Write Results to NetCDF**
   ```r
   write_results_to_netcdf <- function(results, grid, output_file) {
     nc <- nc_open(output_file, write = TRUE)
     
     # Reshape results to lon × lat × time grid
     n_lon <- length(unique(grid$lon))
     n_lat <- length(unique(grid$lat))
     n_time <- ncol(results[[1]]$biomass)
     
     biomass_array <- array(NA, dim = c(n_lon, n_lat, n_time))
     catch_array <- array(NA, dim = c(n_lon, n_lat, n_time))
     
     for (i in 1:nrow(results)) {
       if (results[[i]]$success) {
         lon_idx <- which(unique(grid$lon) == results[[i]]$lon)
         lat_idx <- which(unique(grid$lat) == results[[i]]$lat)
         
         biomass_array[lon_idx, lat_idx, ] <- results[[i]]$biomass
         catch_array[lon_idx, lat_idx, ] <- results[[i]]$catch
       }
     }
     
     # Write variables
     ncvar_put(nc, "tcb", biomass_array)
     ncvar_put(nc, "tc", catch_array)
     
     nc_close(nc)
     
     message("Wrote results to: ", output_file)
   }
   ```

4. **Batch Processing for Large Grids**
   ```r
   # For very large grids, process in chunks
   run_gridded_model_batched <- function(grid, ..., batch_size = 1000) {
     n_cells <- nrow(grid)
     n_batches <- ceiling(n_cells / batch_size)
     
     for (b in 1:n_batches) {
       start_idx <- (b - 1) * batch_size + 1
       end_idx <- min(b * batch_size, n_cells)
       
       cat("Processing batch", b, "of", n_batches, "\n")
       
       batch_grid <- grid[start_idx:end_idx, ]
       batch_results <- run_gridded_model(batch_grid, ...)
       
       # Save intermediate results
       saveRDS(
         batch_results, 
         paste0("outputs/gridded/batch_", sprintf("%03d", b), ".rds")
       )
     }
   }
   ```

**Outputs:**
- `outputs/gridded/zoomss_isimip3a_[model]_[scenario]_tcb.nc`
- `outputs/gridded/zoomss_isimip3a_[model]_[scenario]_tc.nc`

**Tests:**
```r
test_that("Gridded output has correct dimensions", {
  nc <- nc_open("outputs/gridded/zoomss_isimip3a_historical_tcb.nc")
  tcb <- ncvar_get(nc, "tcb")
  nc_close(nc)
  
  expect_equal(length(dim(tcb)), 3)  # lon × lat × time
  expect_true(all(!is.na(tcb[!is.na(tcb)])))  # No NaN in ocean cells
})
```

---

### Task B.7: Output Processing and Aggregation

**Objective:** Aggregate gridded outputs to regional/global summaries and format for ISIMIP

**DBPM Reference:** Fetch and study `scripts/07_calculating_catches_DBPM.py`

**Implementation Steps:**

1. **Aggregate Outputs by Region**
   ```r
   # scripts/07_process_outputs.R
   
   aggregate_to_regions <- function(gridded_output, grid, regions) {
     nc <- nc_open(gridded_output)
     tcb <- ncvar_get(nc, "tcb")  # Total consumer biomass
     tc <- ncvar_get(nc, "tc")    # Total catch
     lon <- ncvar_get(nc, "lon")
     lat <- ncvar_get(nc, "lat")
     time <- ncvar_get(nc, "time")
     nc_close(nc)
     
     # Calculate regional totals
     regional_summary <- lapply(unique(grid$region_id), function(rid) {
       # Get cells in this region
       region_cells <- which(grid$region_id == rid)
       
       # Sum across cells (weighted by cell area if needed)
       regional_biomass <- apply(tcb[region_cells, , ], 2, sum, na.rm = TRUE)
       regional_catch <- apply(tc[region_cells, , ], 2, sum, na.rm = TRUE)
       
       data.frame(
         region_id = rid,
         year = time,
         biomass = regional_biomass,
         catch = regional_catch
       )
     })
     
     regional_df <- do.call(rbind, regional_summary)
     return(regional_df)
   }
   ```

2. **Calculate Global Totals**
   ```r
   calculate_global_totals <- function(gridded_output) {
     nc <- nc_open(gridded_output)
     tcb <- ncvar_get(nc, "tcb")
     tc <- ncvar_get(nc, "tc")
     time <- ncvar_get(nc, "time")
     nc_close(nc)
     
     # Sum across all ocean cells
     global_biomass <- apply(tcb, 3, sum, na.rm = TRUE)
     global_catch <- apply(tc, 3, sum, na.rm = TRUE)
     
     global_df <- data.frame(
       year = time,
       biomass_global = global_biomass,
       catch_global = global_catch
     )
     
     return(global_df)
   }
   ```

3. **Format for ISIMIP Protocol**
   ```r
   format_for_isimip <- function(output_nc, metadata) {
     # ISIMIP requires specific NetCDF attributes
     nc <- nc_open(output_nc, write = TRUE)
     
     # Global attributes
     ncatt_put(nc, 0, "title", "ZooMSS ISIMIP3a Output")
     ncatt_put(nc, 0, "institution", "Your Institution")
     ncatt_put(nc, 0, "source", "ZooMSS v1.0")
     ncatt_put(nc, 0, "contact", "your.email@institution.edu")
     ncatt_put(nc, 0, "climate_forcing", metadata$climate_model)
     ncatt_put(nc, 0, "scenario", metadata$scenario)
     ncatt_put(nc, 0, "creation_date", as.character(Sys.Date()))
     ncatt_put(nc, 0, "Conventions", "CF-1.6")
     
     # Variable attributes
     ncatt_put(nc, "tcb", "standard_name", "total_consumer_biomass_in_sea_water")
     ncatt_put(nc, "tc", "standard_name", "total_catch_of_consumers")
     
     nc_close(nc)
   }
   ```

4. **Create Summary Reports**
   ```r
   create_summary_report <- function(regional_outputs, global_outputs, output_dir) {
     # Time series plots
     p1 <- ggplot(global_outputs, aes(x = year, y = biomass_global)) +
       geom_line() +
       labs(title = "Global Fish Biomass", y = "Biomass (g)", x = "Year") +
       theme_minimal()
     
     p2 <- ggplot(global_outputs, aes(x = year, y = catch_global)) +
       geom_line() +
       labs(title = "Global Catch", y = "Catch (g/yr)", x = "Year") +
       theme_minimal()
     
     # Regional comparison
     p3 <- ggplot(regional_outputs, aes(x = year, y = catch, color = as.factor(region_id))) +
       geom_line() +
       labs(title = "Regional Catch Time Series", y = "Catch", x = "Year") +
       theme_minimal() +
       theme(legend.position = "none")
     
     # Save plots
     ggsave(file.path(output_dir, "global_biomass.pdf"), p1)
     ggsave(file.path(output_dir, "global_catch.pdf"), p2)
     ggsave(file.path(output_dir, "regional_catch.pdf"), p3, width = 12, height = 8)
     
     # Summary statistics
     summary_stats <- global_outputs %>%
       summarise(
         mean_biomass = mean(biomass_global),
         mean_catch = mean(catch_global),
         trend_biomass = coef(lm(biomass_global ~ year))[2],
         trend_catch = coef(lm(catch_global ~ year))[2]
       )
     
     write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"))
   }
   ```

**Outputs:**
- `outputs/regional/regional_summaries.csv`
- `outputs/regional/global_totals.csv`
- `outputs/regional/summary_plots.pdf`
- Formatted NetCDF files compliant with ISIMIP protocol

**Tests:**
```r
test_that("Regional aggregation conserves total", {
  regional <- readRDS("outputs/regional/regional_summaries.rds")
  global <- readRDS("outputs/regional/global_totals.rds")
  
  regional_total <- regional %>%
    group_by(year) %>%
    summarise(total = sum(catch))
  
  expect_equal(regional_total$total, global$catch_global, tolerance = 1e-6)
})
```

---

## Testing Strategy

### Unit Tests

Create tests for each major component:

```r
# tests/testthat/test-climate.R
test_that("Climate data processing works", { ... })

# tests/testthat/test-regions.R
test_that("Region assignment is correct", { ... })

# tests/testthat/test-fishing.R
test_that("Effort time series are valid", { ... })

# tests/testthat/test-calibration.R
test_that("Calibration converges", { ... })

# tests/testthat/test-gridded.R
test_that("Grid cell runner works", { ... })

# tests/testthat/test-outputs.R
test_that("Output aggregation is correct", { ... })
```

### Integration Tests

Test complete workflow:

```r
# tests/testthat/test-integration.R

test_that("Complete workflow runs end-to-end", {
  # Small test case: single region, short time series
  
  # 1. Process climate
  climate <- process_test_climate()
  expect_s3_class(climate, "list")
  
  # 2. Process regions
  regions <- process_test_regions()
  expect_gt(nrow(regions), 0)
  
  # 3. Process fishing
  effort <- process_test_effort()
  expect_true(all(effort >= 0))
  
  # 4. Calibrate
  cal <- calibrate_test_region()
  expect_equal(cal$convergence, 0)
  
  # 5. Run grid
  results <- run_test_grid()
  expect_true(all(sapply(results, function(x) x$success)))
  
  # 6. Aggregate
  summary <- aggregate_test_results(results)
  expect_s3_class(summary, "data.frame")
})
```

---

## Validation and Quality Control

### Validation Checklist

1. **Climate Data:**
   - [ ] Spatial coverage complete
   - [ ] Temporal continuity verified
   - [ ] Values within physical bounds
   - [ ] Compared to observations

2. **Regional Definitions:**
   - [ ] All ocean cells assigned
   - [ ] No gaps or overlaps
   - [ ] Matches official FAO/LME boundaries

3. **Fishing Data:**
   - [ ] Taxonomic mapping verified
   - [ ] Effort time series continuous
   - [ ] Matches Watson catch totals

4. **Calibration:**
   - [ ] All regions converged
   - [ ] Catchability values reasonable (0 < q < 1)
   - [ ] Modeled catch matches observed (R² > 0.7)
   - [ ] Residuals examined for bias

5. **Model Outputs:**
   - [ ] NetCDF format correct
   - [ ] ISIMIP metadata complete
   - [ ] No NaN in ocean cells
   - [ ] Mass balance verified
   - [ ] Compared to other models (e.g., DBPM, FishMIP)

---

## ISIMIP Protocol Compliance

### Required Output Variables

| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| Total Consumer Biomass | `tcb` | g m⁻² | Sum of all fish biomass |
| Total Catch | `tc` | g m⁻² yr⁻¹ | Sum of all catch |
| Total Consumer Biomass Density | `tcbd` | g m⁻³ | Biomass per volume |
| Biomass by Size Class | `b10cm`, `b30cm`, etc. | g m⁻² | Size-specific biomass |

### NetCDF Metadata Requirements

```r
# Required global attributes
ncatt_put(nc, 0, "Conventions", "CF-1.6")
ncatt_put(nc, 0, "title", "ZooMSS ISIMIP3a Output")
ncatt_put(nc, 0, "institution", "Your Institution")
ncatt_put(nc, 0, "source", "ZooMSS v1.0")
ncatt_put(nc, 0, "climate_forcing", "GFDL-ESM4")
ncatt_put(nc, 0, "scenario", "historical")
ncatt_put(nc, 0, "contact", "your.email@institution.edu")

# Variable attributes
ncatt_put(nc, "tcb", "standard_name", "total_consumer_biomass_in_sea_water")
ncatt_put(nc, "tcb", "long_name", "Total Consumer Biomass")
ncatt_put(nc, "tcb", "units", "g m-2")
```

### File Naming Convention

```
zoomss_[climate-model]_[bias-adjustment]_[climate-scenario]_[soc-scenario]_[sens-experiment]_[variable]_[region]_[time-step]_[time-period].nc4

Example:
zoomss_gfdl-esm4_none_historical_histsoc_default_tcb_global_monthly_1950_2014.nc4
```

---

## Performance Optimization

### Memory Management

```r
# For large grids, process in batches
run_gridded_model_batched(grid, batch_size = 1000)

# Use save_catch_by_size = FALSE unless needed
zoomss_model(..., save_catch_by_size = FALSE)

# Clear intermediate results
rm(large_object)
gc()
```

### Parallel Execution

```r
# Optimize number of cores
n_cores <- detectCores() - 1

# Use parallel backend
library(doParallel)
registerDoParallel(n_cores)

# Monitor memory per worker
# Each worker needs ~2GB for ZooMSS run
max_workers <- floor(available_memory_gb / 2)
```

---

## Troubleshooting Guide

### Common Issues

1. **Calibration doesn't converge**
   - Check initial parameter values
   - Try different optimization method (e.g., "Nelder-Mead")
   - Reduce time series length
   - Check for data quality issues

2. **Model runs fail for some cells**
   - Check environmental data for NaNs
   - Verify effort time series length matches
   - Examine error messages in log files

3. **NetCDF output too large**
   - Reduce temporal resolution (annual instead of monthly)
   - Aggregate spatially before writing
   - Use compression: `compression = 9` in `ncvar_def()`

4. **Parallel execution crashes**
   - Reduce number of workers
   - Check memory limits
   - Use batch processing

---

## References and Resources

### ISIMIP Resources
- ISIMIP Homepage: https://www.isimip.org/
- FishMIP: https://www.isimip.org/about/marine-ecosystems-fisheries/
- Data Portal: https://data.isimip.org/
- Protocol Documentation: https://www.isimip.org/protocol/

### Data Sources
- Watson Catch Data: http://www.seararoundus.org/data/#/spatial-catch
- FAO Fishing Areas: https://www.fao.org/fishery/area/search/en
- MARMAP LMEs: https://marineregions.org/
- ISIMIP Climate Data: https://data.isimip.org/

### ZooMSS Package
- GitHub: https://github.com/[your-org]/zoomss
- Documentation: See `vignettes/fishing.Rmd`
- Fishing Implementation: See `docs/ZooMSS_Fishing_Implementation_Assessment_v2.md`

---

## Appendix: Example Complete Workflow

```r
# Complete example workflow script

library(zoomss)
library(ncdf4)
library(sf)
library(dplyr)

# 1. Download and process climate data
climate_data <- download_isimip_climate("GFDL-ESM4", "historical", "data/climate/")
climate_processed <- process_climate_for_zoomss(climate_data, spatial_res = 1.0)

# 2. Setup regions
fao_regions <- download_fao_regions("data/regions/")
grid <- create_spatial_grid(resolution = 1.0, mask_land = TRUE)
grid$region_id <- assign_cells_to_regions(grid, fao_regions)

# 3. Process fishing data
watson_catch <- download_watson_data("data/fishing/")
catch_processed <- process_watson_catch(watson_catch, fao_regions)
effort_timeseries <- create_effort_timeseries(catch_processed, 1950, 2014)

# 4. Calibrate (CRITICAL STEP)
calibration_results <- calibrate_all_regions(
  regions = unique(grid$region_id),
  observed_catch_data = catch_processed,
  climate_data = climate_processed,
  effort_data = effort_timeseries
)

# Validate calibration
for (cal in calibration_results) {
  val <- validate_calibration(cal, catch_processed, climate_processed, effort_timeseries)
  print(val$plot)
}

# 5. Run gridded model
output_file <- "outputs/gridded/zoomss_gfdl-esm4_historical_tcb.nc"
initialize_output_netcdf(grid, climate_processed$time, c("tcb", "tc"), output_file)

results <- run_gridded_model(
  grid = grid,
  climate_data = climate_processed,
  effort_data = effort_timeseries,
  calibration = calibration_results,
  output_file = output_file,
  n_cores = 10
)

# 6. Process and aggregate outputs
regional_summary <- aggregate_to_regions(output_file, grid, fao_regions)
global_summary <- calculate_global_totals(output_file)

# 7. Create reports
create_summary_report(regional_summary, global_summary, "outputs/reports/")

# 8. Format for ISIMIP submission
format_for_isimip(output_file, list(climate_model = "GFDL-ESM4", scenario = "historical"))
```

---

## Success Criteria

Phase B is complete when:

- [ ] All 7 workflow scripts (B.1-B.7) are implemented and tested
- [ ] Calibration achieves R² > 0.7 for >80% of regions
- [ ] Gridded model runs successfully for full historical period
- [ ] Outputs pass ISIMIP protocol validation
- [ ] Documentation is complete (README, vignettes)
- [ ] At least one scenario simulation completed (e.g., historical + ssp370)
- [ ] Results compared to DBPM and other FishMIP models
- [ ] Code is version controlled and reproducible

---

**End of Implementation Guide**

For questions or issues during implementation, refer to:
1. ZooMSS package documentation
2. DBPM reference scripts
3. ISIMIP protocol documentation
4. ZooMSS_Fishing_Implementation_Assessment_v2.md (Phase A decisions)
