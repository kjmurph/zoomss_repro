# ZooMSS Revision Instructions

## Context

These instructions are for fixing documentation errors, a diffusion scaling bug, and establishing a calibration workflow in the ZooMSS R package. The package implements a size-structured marine ecosystem model with an explicit energy budget and fish reproduction system.

The energy budget follows two stages:
1. Ingestion → Defecation (D, prey-quality dependent) + Assimilation (1 - D)
2. Assimilation → Metabolism (f_M) + Growth (K_growth) + Reproduction (R_frac)

Where `f_M + K_growth + R_frac = 1` for all groups.

For zooplankton groups, `R_frac > 0` represents real reproductive costs that are currently "lost" from the system (zooplankton recruitment uses boundary conditions, not explicit reproduction). This is intentional — the energy loss is biologically real and the framework supports future explicit zooplankton reproduction.

---

## Task 1: Fix Documentation in `R/data.R`

### File: `R/data.R`

The roxygen documentation for `GroupInputs` contains an incorrect description of `def_low`. Find and update:

**Find this line:**
```
#'   \item{def_low}{Numeric. Defecation fraction for low-quality (low Carbon) prey (typically 0.50)}
```

**Replace with:**
```
#'   \item{def_low}{Numeric. Defecation fraction for low-quality (low Carbon) prey (typically 0.95, reflecting high water content of gelatinous prey)}
```

No other changes to `R/data.R`.

---

## Task 2: Fix Vignette `vignettes/reproduction.Rmd`

### 2.1 Fix zooplankton energy budget description (line 79)

**Find:**
```
- **Zooplankton** have `repro_on = 0` (all production goes to somatic growth)
```

**Replace with:**
```
- **Zooplankton** have `repro_on = 0` (reproduction is not explicitly modelled). Zooplankton groups with `K_growth < 0.50` have a non-zero `R_frac` representing real reproductive costs that are lost from the modelled system — zooplankton recruitment is instead maintained via boundary conditions. This energy loss is biologically realistic (copepods invest 30–60% of production in egg production) and the framework supports future explicit zooplankton reproduction.
```

### 2.2 Fix zooplankton growth rate comment (lines 188–189)

**Find:**
```
# Growth rates for zooplankton (higher K_growth = 0.50)
```

**Replace with:**
```
# Growth rates for zooplankton (K_growth ranges from 0.25 to 0.50 across groups)
```

### 2.3 Fix fish growth rate comment (line 193)

**Find:**
```
# Growth rates for fish (lower K_growth = 0.36, with reproduction)
```

**Replace with:**
```
# Growth rates for fish (K_growth = 0.30, with R_frac allocated to reproduction for mature individuals)
```

### 2.4 Fix the summary table (lines 812–818)

**Find:**
```
### Key Differences: Zooplankton vs Fish

| Aspect | Zooplankton | Fish |
|--------|-------------|------|
| K_growth | 0.25–0.50 (group-specific) | 0.30 (all fish groups) |
| R_frac | 0.00 (not used) | 0.20 (= 1 - f_M - K_growth) |
| repro_on | 0 (disabled) | 1 (enabled) |
| Recruitment | Boundary condition | SSB-driven |
```

**Replace with:**
```
### Key Differences: Zooplankton vs Fish

| Aspect | Zooplankton | Fish |
|--------|-------------|------|
| K_growth | 0.25–0.50 (group-specific) | 0.30 (all fish groups) |
| R_frac | 0.00–0.25 (group-specific; energy lost from system as implicit reproductive cost) | 0.20 (= 1 - f_M - K_growth; drives explicit reproduction) |
| repro_on | 0 (no explicit reproduction) | 1 (enabled) |
| Recruitment | Boundary condition (independent of energetics) | SSB-driven via R_frac and maturity ogive |

**Note on zooplankton R_frac:** Groups with K_growth < 0.50 have R_frac > 0, representing real energy allocated to reproduction that is not recycled back into the modelled population. This is biologically realistic — zooplankton genuinely invest energy in reproduction — but the current boundary-condition recruitment does not respond to this investment. The energy budget framework supports future implementation of explicit zooplankton reproduction.
```

### 2.5 Fix summary point 1 (lines 794–795)

**Find:**
```
1. **Energy Budget Components**: The model explicitly tracks non-nutritive fraction (defecation + water/structural losses), metabolism, growth, and reproduction fractions. Energy budget parameters (`f_M`, `K_growth`, `R_frac`) are group-specific, ranging from K_growth = 0.25 (Euphausiids) to 0.50 (Jellyfish).
```

**Replace with:**
```
1. **Energy Budget Components**: The model explicitly tracks non-nutritive fraction (defecation + water/structural losses), metabolism, growth, and reproduction fractions. Energy budget parameters (`f_M`, `K_growth`, `R_frac`) are group-specific, with K_growth ranging from 0.25 (Euphausiids) to 0.50 (Jellyfish). For zooplankton, R_frac represents an implicit reproductive cost lost from the system (0.00–0.25 depending on group); for fish, R_frac drives explicit SSB-based recruitment.
```

### 2.6 Add Scenario A/B documentation

Insert a new subsection after the "Energy Budget Allocation" plot (after line 120). Place this before the "Defecation by Prey Quality" subheading:

```markdown
#### Energy Budget Scenarios

The model supports two energy budget scenarios, controlled by `energy_budget_scenario` in `zoomss_params()`:

- **Scenario A (default):** Full energy budget for all groups. Zooplankton groups with `K_growth < 0.50` have `R_frac > 0`, representing reproductive energy costs that are lost from the modelled system. This is biologically realistic but means not all assimilated energy contributes to modelled growth.

- **Scenario B:** Zooplankton `R_frac` is redistributed to `K_growth` (i.e., `R_frac = 0` for all zooplankton). Fish groups are unaffected. This preserves more energy in the modelled system and may simplify calibration against observed biomass patterns.

The choice of scenario affects zooplankton biomass levels and should be documented when reporting model results.
```

---

## Task 3: Fix the Diffusion Scaling Bug in `R/zoomss_run.R`

### Problem

In the MvF-D framework, diffusion scales as growth². Currently, diffusion is uniformly scaled by `K_growth²` for all groups:

```r
diff <- sweep(diff, 1, K_growth^2, '*')
```

But for fish with `repro_on == 1`, immature individuals have an effective growth fraction higher than K_growth because R_frac energy is redirected to growth:

```r
gg[fg, ] <- K_growth[fg] * gg_total[fg, ] + R_frac[fg] * gg_total[fg, ] * (1 - mat_ogive[fg, ])
```

The effective growth fraction for immature fish is `K_growth + R_frac × (1 - mat_ogive)`, which varies across sizes. The diffusion term must match this size-dependent effective growth fraction.

### Fix

In `R/zoomss_run.R`, find the diffusion partitioning section. The current code is:

```r
    # Partition diffusion by K_growth^2 for consistency with growth partitioning
    # (diffusion scales as growth^2 in the MvF-D framework)
    diff <- sweep(diff, 1, K_growth^2, '*')
```

Replace this entire block with:

```r
    # ==========================================================================
    # DIFFUSION PARTITIONING
    # ==========================================================================
    # Diffusion scales as growth^2 in the MvF-D framework.
    # For most groups, the effective growth fraction is simply K_growth.
    # For fish with repro_on = 1, immature individuals redirect R_frac to
    # growth, so their effective growth fraction is size-dependent:
    #   eff_K = K_growth + R_frac * (1 - mat_ogive)
    # Diffusion must use eff_K^2 to remain consistent with the actual growth rate.

    # Start with K_growth^2 for all groups (correct for zooplankton and fish without reproduction)
    eff_K_sq <- matrix(rep(K_growth^2, each = ngrid), nrow = ngrps, ncol = ngrid, byrow = TRUE)

    # Adjust for fish with active reproduction: use size-dependent effective growth fraction
    for (f in 1:num_fish) {
      fg <- fish_grps[f]
      if (repro_on[fg] == 1 && R_frac[fg] > 0) {
        eff_K <- K_growth[fg] + R_frac[fg] * (1 - mat_ogive[fg, ])
        eff_K_sq[fg, ] <- eff_K^2
      }
    }

    # Apply size-dependent diffusion scaling
    diff <- diff * eff_K_sq
```

### Verification

After applying this fix, confirm that:
1. For zooplankton groups, `eff_K_sq` equals `K_growth^2` (unchanged behaviour)
2. For fish with `repro_on = 0`, `eff_K_sq` equals `K_growth^2` (unchanged behaviour)
3. For fish with `repro_on = 1`, `eff_K_sq` smoothly transitions from `(K_growth + R_frac)^2` at small sizes to `K_growth^2` at large sizes, following the maturity ogive
4. Model still reaches stable equilibrium under constant environmental forcing

---

## Task 4: Expose Energy Budget Scenario as Parameter

### File: `R/zoomss_params.R`

The energy budget scenario is currently hard-coded. Expose it as a parameter.

**In the function signature**, find:

```r
zoomss_params <- function(Groups, input_params, isave){
```

Replace with:

```r
zoomss_params <- function(Groups, input_params, isave, energy_budget_scenario = "A"){
```

**Then find the hard-coded scenario line:**

```r
  energy_budget_scenario <- "A"  # Toggle: "A" or "B"
```

Replace with:

```r
  # Validate energy budget scenario
  energy_budget_scenario <- match.arg(energy_budget_scenario, choices = c("A", "B"))
```

### File: `R/zoomss_model.R`

Pass the parameter through. Find:

```r
zoomss_model <- function(input_params, Groups = NULL, isave = 1){
```

Replace with:

```r
zoomss_model <- function(input_params, Groups = NULL, isave = 1, energy_budget_scenario = "A"){
```

And find:

```r
  param <- zoomss_params(Groups, input_params, isave) # Set up parameter list
```

Replace with:

```r
  param <- zoomss_params(Groups, input_params, isave, energy_budget_scenario) # Set up parameter list
```

---

## Task 5: Add Validation for repro_eff When repro_on = 1

### File: `R/zoomss_groups.R`

In the `validateGroups()` function, after the existing `repro_on` validation block (around the line checking `repro_on must be 0 or 1`), add:

```r
  # Check that fish with reproduction enabled have positive reproductive efficiency
  fish_repro_on <- groups$repro_on == 1
  if (any(fish_repro_on)) {
    assertthat::assert_that(all(groups$repro_eff[fish_repro_on] > 0),
                           msg = "repro_eff must be > 0 for groups with repro_on = 1 (otherwise recruitment will be zero)")
  }
```

---

## Task 6: Calibration and Benchmarking Strategy

### Overview

The goal is to calibrate the revised ZooMSS (with explicit energy budget and fish reproduction) against the original ZooMSS (https://github.com/MathMarEcol/zoomss) by matching zooplankton community composition across a productivity gradient. The original model's behaviour serves as the calibration target because it has been validated against global observational data (Heneghan et al. 2020).

**Workflow:**
1. Generate baseline outputs from the original ZooMSS (Section 6.1)
2. Run sensitivity analysis to identify the most influential parameters (Section 6.2)
3. Define the objective function using sensitivity-informed parameter selection (Section 6.3)
4. Run optimisation with parallelised model evaluations (Section 6.4)
5. Evaluate and visualise calibration results (Section 6.5)

**Key constraints:**
- All simulations run for **400 years minimum** to ensure oscillating steady state
- Steady-state metrics averaged over the **final 100 years**
- Full chlorophyll gradient (0.05–10 mg/m³) maintained throughout — do not reduce
- Chlorophyll levels within each objective evaluation are run in parallel via `future.apply`

### 6.1 Generate Baseline from Original ZooMSS

Write a script `calibration/01_generate_baseline.R` that:

1. Installs the original zoomss package from GitHub:
```r
# remotes::install_github("MathMarEcol/zoomss")
```

2. Runs the original model across a chlorophyll gradient at constant temperature (15°C):
```r
chl_levels <- c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)
```

3. For each chlorophyll level, runs to oscillating steady state (400 years minimum) and extracts time-averaged values over the final 100 years:
   - Time-averaged biomass per functional group
   - Zooplankton group proportions (fraction of total zooplankton biomass)
   - Fish biomass by group
   - Total community size spectrum slope
   - Zooplankton:Fish biomass ratio

4. Saves the baseline as an RDS file:
```r
baseline <- list(
  chl_levels = chl_levels,
  zoo_proportions = zoo_prop_matrix,   # matrix: n_chl x n_zoo_groups
  fish_biomass = fish_biomass_matrix,  # matrix: n_chl x n_fish_groups
  total_biomass = total_biomass_vec,   # vector: n_chl
  zoo_fish_ratio = zoo_fish_ratio_vec  # vector: n_chl
)
saveRDS(baseline, "calibration/baseline_original_zoomss.rds")
```

### 6.2 Sensitivity Analysis (Run First)

Before optimisation, run a one-at-a-time (OAT) sensitivity analysis to identify which parameters most strongly influence the calibration targets. This avoids wasting optimisation effort on parameters the model is insensitive to, and reveals potential parameter interactions.

Write `calibration/02_sensitivity_analysis.R`:

```r
library(zoomss)
library(ggplot2)
library(patchwork)
library(future.apply)

plan(multisession, workers = parallelly::availableCores() - 1)

# Load baseline
baseline <- readRDS("calibration/baseline_original_zoomss.rds")
chl_levels <- baseline$chl_levels

# ── Define candidate parameters and perturbation ranges ──
# Each parameter is varied across a range while all others are held at defaults.
# The range should span biologically plausible values.

sensitivity_params <- list(
  f_M = list(
    default = 0.50,
    range = seq(0.30, 0.70, by = 0.05),
    description = "Metabolic fraction of assimilated energy"
  ),
  K_growth_zoo_base = list(
    default = 0.35,
    range = seq(0.15, 0.49, by = 0.05),
    description = "Zooplankton growth fraction (base, scaled per group)"
  ),
  K_growth_fish = list(
    default = 0.30,
    range = seq(0.15, 0.45, by = 0.05),
    description = "Fish growth fraction"
  ),
  repro_eff = list(
    default = 0.002,
    range = c(1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2),
    description = "Fish reproductive efficiency (egg-to-recruit survival)"
  ),
  def_low = list(
    default = 0.95,
    range = seq(0.50, 0.95, by = 0.05),
    description = "Defecation fraction for low-Carbon prey"
  ),
  def_high = list(
    default = 0.30,
    range = seq(0.15, 0.45, by = 0.05),
    description = "Defecation fraction for high-Carbon prey"
  )
)

# ── Helper: run model across chl gradient and extract diagnostics ──
run_gradient <- function(Groups, chl_levels) {
  # Parallelise across chl levels
  results <- future_lapply(chl_levels, function(chl) {
    env <- createInputParams(
      time = seq(0, 400, by = 0.1),
      sst = 15,
      chl = chl
    )
    mdl <- zoomss_model(input_params = env, Groups = Groups, isave = 10)
    avg_biomass <- averageTimeSeries(mdl, var = "biomass", n_years = 100)
    group_biomass <- rowSums(avg_biomass)
    list(group_biomass = group_biomass)
  }, future.seed = TRUE)

  zoo_idx <- which(Groups$Type == "Zooplankton")
  fish_idx <- which(Groups$Type == "Fish")
  n_chl <- length(chl_levels)

  zoo_prop <- matrix(NA, n_chl, length(zoo_idx))
  fish_bm <- matrix(NA, n_chl, length(fish_idx))
  total_bm <- numeric(n_chl)

  for (i in seq_along(chl_levels)) {
    gb <- results[[i]]$group_biomass
    zoo_bm <- gb[zoo_idx]
    zoo_prop[i, ] <- zoo_bm / sum(zoo_bm)
    fish_bm[i, ] <- gb[fish_idx]
    total_bm[i] <- sum(gb)
  }

  list(zoo_proportions = zoo_prop, fish_biomass = fish_bm, total_biomass = total_bm)
}

# ── Helper: apply a parameter perturbation to Groups ──
apply_perturbation <- function(Groups_default, param_name, value) {
  Groups <- Groups_default
  zoo_idx <- which(Groups$Type == "Zooplankton")
  fish_idx <- which(Groups$Type == "Fish")

  switch(param_name,
    f_M = {
      Groups$f_M <- value
    },
    K_growth_zoo_base = {
      default_K_zoo <- Groups$K_growth[zoo_idx]
      K_relative <- default_K_zoo / mean(default_K_zoo)
      Groups$K_growth[zoo_idx] <- pmin(pmax(value * K_relative, 0.05), 0.49)
    },
    K_growth_fish = {
      Groups$K_growth[fish_idx] <- value
    },
    repro_eff = {
      Groups$repro_eff[fish_idx] <- value
    },
    def_low = {
      Groups$def_low <- value
    },
    def_high = {
      Groups$def_high <- value
    }
  )

  # Validate energy budget closure
  R_frac <- 1 - Groups$f_M - Groups$K_growth
  if (any(R_frac < 0)) return(NULL)  # Invalid combination

  Groups
}

# ── Run OAT sensitivity analysis ──
Groups_default <- getGroups()
sensitivity_results <- list()

for (param_name in names(sensitivity_params)) {
  param_info <- sensitivity_params[[param_name]]
  cat("\n=== Sensitivity analysis for:", param_name, "===\n")

  param_results <- list()
  for (val in param_info$range) {
    cat("  Running", param_name, "=", val, "...")

    Groups_mod <- apply_perturbation(Groups_default, param_name, val)
    if (is.null(Groups_mod)) {
      cat(" SKIPPED (invalid energy budget)\n")
      next
    }

    tryCatch({
      res <- run_gradient(Groups_mod, chl_levels)
      res$param_value <- val
      param_results <- c(param_results, list(res))
      cat(" done\n")
    }, error = function(e) {
      cat(" ERROR:", e$message, "\n")
    })
  }

  sensitivity_results[[param_name]] <- param_results
}

saveRDS(sensitivity_results, "calibration/sensitivity_results.rds")

# ── Compute sensitivity metrics ──
# For each parameter, compute the range of effect on key outputs:
# 1. Mean absolute change in zoo proportions across gradient
# 2. Mean absolute change in log10 total biomass
# 3. Mean absolute change in log10 zoo:fish ratio

sensitivity_summary <- data.frame(
  parameter = character(),
  delta_zoo_prop = numeric(),   # Sensitivity of zoo proportions
  delta_total_bm = numeric(),   # Sensitivity of total biomass
  delta_zf_ratio = numeric(),   # Sensitivity of zoo:fish ratio
  stringsAsFactors = FALSE
)

for (param_name in names(sensitivity_results)) {
  results_list <- sensitivity_results[[param_name]]
  if (length(results_list) < 2) next

  # Compare each perturbation to the default run
  default_idx <- which(sapply(results_list, function(x) {
    abs(x$param_value - sensitivity_params[[param_name]]$default) < 1e-8
  }))

  if (length(default_idx) == 0) {
    # Use midpoint as reference
    default_idx <- ceiling(length(results_list) / 2)
  }

  ref <- results_list[[default_idx]]
  deltas_prop <- numeric(length(results_list))
  deltas_bm <- numeric(length(results_list))
  deltas_ratio <- numeric(length(results_list))

  for (j in seq_along(results_list)) {
    res <- results_list[[j]]
    deltas_prop[j] <- mean(abs(res$zoo_proportions - ref$zoo_proportions), na.rm = TRUE)
    deltas_bm[j] <- mean(abs(log10(res$total_biomass) - log10(ref$total_biomass)), na.rm = TRUE)

    ref_ratio <- rowSums(ref$zoo_proportions) / rowSums(ref$fish_biomass)
    res_ratio <- rowSums(res$zoo_proportions) / rowSums(res$fish_biomass)
    deltas_ratio[j] <- mean(abs(log10(res_ratio) - log10(ref_ratio)), na.rm = TRUE)
  }

  sensitivity_summary <- rbind(sensitivity_summary, data.frame(
    parameter = param_name,
    delta_zoo_prop = max(deltas_prop),
    delta_total_bm = max(deltas_bm),
    delta_zf_ratio = max(deltas_ratio)
  ))
}

cat("\n=== SENSITIVITY SUMMARY ===\n")
cat("(larger values = higher sensitivity = better calibration parameter)\n\n")
print(sensitivity_summary[order(-sensitivity_summary$delta_zoo_prop), ])

# ── Save summary ──
saveRDS(sensitivity_summary, "calibration/sensitivity_summary.rds")
```

**How to interpret**: Parameters with large `delta_zoo_prop` are strong candidates for the optimisation. Parameters with negligible effect can be fixed at their defaults to reduce dimensionality. If two parameters show similar sensitivity patterns, they may be correlated (non-identifiable) — check by examining their individual response curves.

**Expected outcomes**: Based on model structure, `f_M` and `K_growth_zoo_base` should dominate zooplankton proportions (they directly control zooplankton growth rates). `repro_eff` should strongly affect fish biomass and zoo:fish ratio but less so zooplankton proportions. `def_low` and `def_high` affect all groups through prey-quality scaling — if sensitivity is low, they can be excluded from the optimisation.

Use the results of this analysis to select the final set of tuning parameters for optimisation. Only parameters showing meaningful sensitivity to the calibration targets should be included.

### 6.3 Define the Objective Function

Write `calibration/03_objective_function.R`:

The objective function takes a parameter vector, runs the revised model across the same chlorophyll gradient, and returns the sum-of-squared differences in zooplankton proportions compared to the baseline. Chlorophyll levels are evaluated in parallel via `future.apply` — ensure `plan(multisession)` is set before calling this function.

```r
#' Calibration Objective Function for Revised ZooMSS
#'
#' @param par Named numeric vector of parameters to optimise:
#'   - f_M: Metabolic fraction (single value applied to all groups, or vector)
#'   - K_growth_zoo: Zooplankton growth fractions (vector, one per zoo group)
#'   - K_growth_fish: Fish growth fraction (single value for all fish)
#'   - repro_eff: Fish reproductive efficiency (single value for all fish)
#' @param baseline List from 01_generate_baseline.R
#' @param chl_levels Chlorophyll gradient to evaluate
#' @param Groups Default Groups data frame (modified by par)
#' @param verbose Logical, print progress
#' @return Scalar objective value (lower = better fit)

calibration_objective <- function(par, baseline, chl_levels, Groups, verbose = FALSE) {

  # ── Unpack and apply parameters ──
  # f_M: shared across all groups (single value)
  Groups$f_M <- par["f_M"]

  # K_growth for zooplankton: group-specific
  zoo_idx <- which(Groups$Type == "Zooplankton")
  n_zoo <- length(zoo_idx)

  # Option A: Single K_growth_zoo applied to all zooplankton
  # Groups$K_growth[zoo_idx] <- par["K_growth_zoo"]

  # Option B: Group-specific K_growth for zooplankton (preferred, more degrees of freedom)
  # Use a scaling approach: par contains a base value and relative scalars
  K_base_zoo <- par["K_growth_zoo_base"]
  # Scale individual groups relative to base (maintain rank order from defaults)
  default_K_zoo <- Groups$K_growth[zoo_idx]
  K_relative <- default_K_zoo / mean(default_K_zoo)  # preserve relative pattern
  Groups$K_growth[zoo_idx] <- K_base_zoo * K_relative
  # Clamp to valid range
  Groups$K_growth[zoo_idx] <- pmin(pmax(Groups$K_growth[zoo_idx], 0.05), 0.49)

  # K_growth for fish
  fish_idx <- which(Groups$Type == "Fish")
  Groups$K_growth[fish_idx] <- par["K_growth_fish"]

  # Reproductive efficiency for fish
  Groups$repro_eff[fish_idx] <- par["repro_eff"]

  # ── Validate energy budget closure ──
  R_frac <- 1 - Groups$f_M - Groups$K_growth
  if (any(R_frac < 0) || any(R_frac > 1)) {
    return(1e6)  # Penalty for invalid parameter combinations
  }

  # ── Run model across chlorophyll gradient (parallelised) ──
  n_chl <- length(chl_levels)

  # Each chl level is independent — run in parallel for speed
  run_single_chl <- function(chl, Groups) {
    tryCatch({
      env <- createInputParams(
        time = seq(0, 400, by = 0.1),
        sst = 15,
        chl = chl
      )

      mdl <- zoomss_model(
        input_params = env,
        Groups = Groups,
        isave = 10  # coarser saving for speed during calibration
      )

      # Extract steady-state biomass (average final 100 years)
      avg_biomass <- averageTimeSeries(mdl, var = "biomass", n_years = 100)
      group_biomass <- rowSums(avg_biomass)
      list(group_biomass = group_biomass, success = TRUE)
    }, error = function(e) {
      list(group_biomass = NULL, success = FALSE, error = e$message)
    })
  }

  # Run in parallel using future.apply (set up plan externally)
  if (requireNamespace("future.apply", quietly = TRUE)) {
    chl_results <- future.apply::future_lapply(
      chl_levels, run_single_chl, Groups = Groups,
      future.seed = TRUE
    )
  } else {
    chl_results <- lapply(chl_levels, run_single_chl, Groups = Groups)
  }

  # Unpack results
  zoo_proportions <- matrix(NA, nrow = n_chl, ncol = n_zoo)
  fish_biomass <- matrix(NA, nrow = n_chl, ncol = length(fish_idx))
  total_biomass <- numeric(n_chl)

  for (i in seq_along(chl_levels)) {
    res <- chl_results[[i]]
    if (res$success) {
      gb <- res$group_biomass

      zoo_biomass <- gb[zoo_idx]
      zoo_proportions[i, ] <- zoo_biomass / sum(zoo_biomass)
      fish_biomass[i, ] <- gb[fish_idx]
      total_biomass[i] <- sum(gb)
    } else {
      if (verbose) cat("Error at chl =", chl_levels[i], ":", res$error, "\n")
    }
  }

  # ── Calculate objective ──
  # If any runs failed, return large penalty
  if (any(is.na(zoo_proportions))) {
    return(1e6 - sum(!is.na(zoo_proportions)))  # Slightly less penalty for partial success
  }

  # Component 1: Zooplankton proportions across gradient (primary target)
  # Weight = 1.0 (this is the main calibration target)
  ss_zoo_prop <- sum((zoo_proportions - baseline$zoo_proportions)^2)

  # Component 2: Zoo:Fish biomass ratio across gradient (secondary target)
  # This constrains the overall energy transfer from zoo to fish
  zoo_fish_ratio <- rowSums(zoo_proportions * total_biomass) /
                    rowSums(fish_biomass)
  # Use log-ratio to handle scale differences
  ss_ratio <- sum((log10(zoo_fish_ratio) - log10(baseline$zoo_fish_ratio))^2,
                  na.rm = TRUE)

  # Component 3: Total biomass magnitude (tertiary target)
  # Prevents solutions where proportions match but absolute values are wrong
  ss_total <- sum((log10(total_biomass) - log10(baseline$total_biomass))^2,
                  na.rm = TRUE)

  # Weighted objective
  w_prop <- 1.0    # Primary: zooplankton community composition
  w_ratio <- 0.3   # Secondary: zoo:fish ratio
  w_total <- 0.1   # Tertiary: absolute biomass magnitude

  objective <- w_prop * ss_zoo_prop + w_ratio * ss_ratio + w_total * ss_total

  if (verbose) {
    cat(sprintf("f_M=%.3f K_zoo=%.3f K_fish=%.3f repro_eff=%.4f | obj=%.4f (prop=%.4f ratio=%.4f total=%.4f)\n",
                par["f_M"], par["K_growth_zoo_base"], par["K_growth_fish"],
                par["repro_eff"], objective, ss_zoo_prop, ss_ratio, ss_total))
  }

  return(objective)
}
```

### 6.4 Run the Optimisation

Write `calibration/04_run_calibration.R`:

```r
library(zoomss)
library(future.apply)

# Set up parallel execution across chlorophyll levels
plan(multisession, workers = parallelly::availableCores() - 1)

# Load baseline from original model
baseline <- readRDS("calibration/baseline_original_zoomss.rds")
chl_levels <- baseline$chl_levels

# Load sensitivity analysis results to inform parameter selection
sensitivity <- readRDS("calibration/sensitivity_summary.rds")
cat("Sensitivity ranking (by zoo proportions):\n")
print(sensitivity[order(-sensitivity$delta_zoo_prop), ])
# Use this to decide which parameters to include in optimisation.
# Parameters with negligible delta_zoo_prop can be fixed at defaults.

# Load default Groups as starting point
Groups <- getGroups()

# Source the objective function
source("calibration/03_objective_function.R")

# ── Define parameter bounds ──
# These bounds enforce the energy budget constraint f_M + K_growth <= 1
par_lower <- c(
  f_M = 0.30,              # Minimum metabolic fraction
  K_growth_zoo_base = 0.15, # Minimum zoo growth fraction
  K_growth_fish = 0.15,     # Minimum fish growth fraction
  repro_eff = 1e-6          # Minimum reproductive efficiency
)

par_upper <- c(
  f_M = 0.70,              # Maximum metabolic fraction
  K_growth_zoo_base = 0.50, # Maximum zoo growth fraction (f_M + K = 1 means R_frac = 0)
  K_growth_fish = 0.45,     # Maximum fish growth fraction
  repro_eff = 0.01          # Maximum reproductive efficiency
)

# Starting values (from current defaults)
par_start <- c(
  f_M = 0.50,
  K_growth_zoo_base = 0.35,
  K_growth_fish = 0.30,
  repro_eff = 0.002
)

# ── Option A: Nelder-Mead (derivative-free, good for noisy objectives) ──
result_nm <- optim(
  par = par_start,
  fn = calibration_objective,
  method = "Nelder-Mead",
  baseline = baseline,
  chl_levels = chl_levels,
  Groups = Groups,
  verbose = TRUE,
  control = list(
    maxit = 200,
    trace = 1,
    reltol = 1e-4
  )
)

# ── Option B: L-BFGS-B (box-constrained, more efficient if objective is smooth) ──
result_bfgs <- optim(
  par = par_start,
  fn = calibration_objective,
  method = "L-BFGS-B",
  lower = par_lower,
  upper = par_upper,
  baseline = baseline,
  chl_levels = chl_levels,
  Groups = Groups,
  verbose = TRUE,
  control = list(
    maxit = 100,
    trace = 1
  )
)

# ── Option C: DEoptim (global optimiser, slower but avoids local minima) ──
# Recommended for initial exploration because the objective landscape may be
# multimodal (multiple parameter combinations produce similar patterns)
if (requireNamespace("DEoptim", quietly = TRUE)) {
  result_de <- DEoptim::DEoptim(
    fn = calibration_objective,
    lower = par_lower,
    upper = par_upper,
    baseline = baseline,
    chl_levels = chl_levels,
    Groups = Groups,
    verbose = TRUE,
    control = DEoptim::DEoptim.control(
      NP = 40,         # Population size (10x number of parameters)
      itermax = 100,    # Maximum iterations
      trace = 10,       # Print every 10 iterations
      parallelType = 1  # Parallel evaluation (if available)
    )
  )
}

# ── Save results ──
saveRDS(list(
  nelder_mead = result_nm,
  lbfgsb = result_bfgs,
  deoptim = if (exists("result_de")) result_de else NULL,
  baseline = baseline,
  par_bounds = list(lower = par_lower, upper = par_upper)
), "calibration/calibration_results.rds")
```

### 6.5 Evaluate and Visualise Calibration Results

Write `calibration/05_evaluate_calibration.R`:

```r
library(zoomss)
library(ggplot2)
library(patchwork)
library(future.apply)

plan(multisession, workers = parallelly::availableCores() - 1)

# Load results
results <- readRDS("calibration/calibration_results.rds")
baseline <- results$baseline

# Pick best result (compare across methods)
best <- results$deoptim  # or results$nelder_mead / results$lbfgsb
best_par <- best$optim$bestmem  # DEoptim syntax; use best$par for optim()

cat("Best parameters:\n")
print(best_par)

# ── Run revised model with best parameters ──
Groups <- getGroups()
# Apply best parameters to Groups (same logic as in objective function)
Groups$f_M <- best_par["f_M"]
zoo_idx <- which(Groups$Type == "Zooplankton")
fish_idx <- which(Groups$Type == "Fish")
K_base_zoo <- best_par["K_growth_zoo_base"]
default_K_zoo <- Groups$K_growth[zoo_idx]
K_relative <- default_K_zoo / mean(default_K_zoo)
Groups$K_growth[zoo_idx] <- pmin(pmax(K_base_zoo * K_relative, 0.05), 0.49)
Groups$K_growth[fish_idx] <- best_par["K_growth_fish"]
Groups$repro_eff[fish_idx] <- best_par["repro_eff"]

# Print the calibrated energy budget
cat("\nCalibrated energy budget:\n")
R_frac <- 1 - Groups$f_M - Groups$K_growth
print(data.frame(
  Species = Groups$Species,
  f_M = Groups$f_M,
  K_growth = Groups$K_growth,
  R_frac = round(R_frac, 3),
  repro_eff = Groups$repro_eff
))

# Run across gradient (parallelised)
chl_levels <- baseline$chl_levels
revised_results <- future_lapply(chl_levels, function(chl) {
  env <- createInputParams(time = seq(0, 400, by = 0.1), sst = 15, chl = chl)
  zoomss_model(input_params = env, Groups = Groups, isave = 10)
}, future.seed = TRUE)

# ── Extract and compare ──
# ... (extract zoo proportions, fish biomass, total biomass from revised_results
#      and plot side-by-side with baseline)

# Key diagnostic plots:
# 1. Zooplankton proportions: baseline vs calibrated (stacked bars at each chl level)
# 2. Zoo:Fish ratio: baseline vs calibrated (line plot across gradient)
# 3. Size spectra comparison at 3 representative chl levels (0.1, 1.0, 10.0)
# 4. Fish SSB and recruitment across gradient (new model only, no baseline equivalent)
```

### 6.6 Calibration Strategy Notes

**Why these parameters (confirm with sensitivity analysis):**

The table below lists candidate parameters and their expected sensitivity. The sensitivity analysis (Section 6.2) should be used to confirm which parameters warrant inclusion in optimisation — only include parameters that the OAT analysis shows have meaningful influence on the calibration targets.

| Parameter | Role | Expected Sensitivity |
|-----------|------|-------------|
| `f_M` | Controls total energy available for growth + reproduction. Higher f_M = less growth. | High — directly scales all growth rates |
| `K_growth_zoo_base` | Controls how much assimilated energy goes to zooplankton somatic growth. Indirectly sets R_frac (the "lost" energy). | High — determines zooplankton biomass accumulation |
| `K_growth_fish` | Controls fish somatic growth rate. With reproduction, lower K_growth means more R_frac for reproduction. | Moderate — interacts with repro_eff |
| `repro_eff` | Egg-to-recruit survival. Controls how efficiently reproductive investment translates to new fish. | High for fish dynamics — tiny values (1e-4 to 1e-2) make large differences |
| `def_low` | Defecation for low-Carbon prey. Controls energy loss when eating gelatinous prey. | Unknown — test with sensitivity analysis |
| `def_high` | Defecation for high-Carbon prey. Controls baseline energy loss. | Unknown — test with sensitivity analysis |

**Practical tips:**

- **Sensitivity first**: Always run the sensitivity analysis (Section 6.2) before optimisation. This identifies which parameters actually matter and avoids wasting compute on insensitive parameters. If the sensitivity analysis shows that (e.g.) `def_low` has negligible impact on zooplankton proportions, exclude it from optimisation entirely.

- **Selecting tuning parameters**: Only include parameters in the optimisation that the sensitivity analysis identifies as having meaningful influence on the calibration targets. Fewer tuning parameters = faster convergence and more identifiable results.

- **Simulation length**: All simulations must run for **400 years minimum** to ensure the model reaches an oscillating steady state. Average over the **final 100 years** to smooth out any residual oscillations. Do not reduce simulation length — under-converged simulations produce misleading objective values that will derail the optimiser.

- **Full chlorophyll gradient**: Always use the complete gradient `c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)`. The extremes (oligotrophic and eutrophic) are where community composition shifts most dramatically and provide the strongest constraints on parameters. Reducing the gradient loses critical information.

- **Parallelise across chl levels**: Each chlorophyll level is independent. Set up `future.apply` at the top of your session and all 8 chl levels run simultaneously:
  ```r
  library(future.apply)
  plan(multisession, workers = parallelly::availableCores() - 1)
  ```
  With 8 cores, wall-clock time per objective evaluation drops from ~8× single-run time to ~1× single-run time. At 400 years with dt=0.1 and isave=10, each individual run takes ~30–90 seconds depending on hardware, so a parallelised objective evaluation takes ~1–2 minutes. Budget 2–4 hours for DEoptim with 100 iterations on a multicore machine.

- **Two-phase optimisation**: Based on sensitivity analysis results, consider splitting:
  - **Phase 1**: Optimise zooplankton parameters (`f_M`, `K_growth_zoo_base`) with fish parameters fixed. The zooplankton community structure is largely set by these two parameters.
  - **Phase 2**: Fix zooplankton parameters at Phase 1 optima, then optimise fish parameters (`K_growth_fish`, `repro_eff`). Fish dynamics respond to zooplankton structure, so this ordering is natural.
  This halves the dimensionality at each phase, dramatically improving convergence.

- **Sensitivity analysis post-calibration**: After optimisation, re-run the sensitivity analysis centred on the optimised parameter values (±20%). This checks for parameter identifiability (flat objective surface = non-identifiable) and correlations between parameters. Plot 2D objective surfaces for parameter pairs to visualise interactions.

---

## Task 7: Comparison Document Structure

Create `vignettes/model_comparison.Rmd` with this structure:

1. **Introduction**: Purpose of comparison, what changed between versions
2. **Parameter Mapping Table**: Old GGE → new f_M, K_growth, R_frac, def_high, def_low
3. **Null Test**: Run revised model with Scenario B (R_frac=0 for zoo) and def_high = def_low = (1 - old_GGE) to verify it approximates the original. Any remaining differences are due to the prey-quality-dependent defecation
4. **Productivity Gradient Comparison**: Side-by-side zooplankton proportions, total biomass, zoo:fish ratios across chl = 0.05 to 10 mg/m³
5. **Temperature Gradient Comparison**: Same diagnostics across SST = 5 to 25°C
6. **Fish Dynamics**: SSB, recruitment, biomass stability in new model vs old boundary condition
7. **Calibrated vs Original**: Results from the optimisation in Task 6
8. **Discussion**: What the explicit energy budget adds (prey quality effects, reproduction, future extensibility) vs what it costs (more parameters, calibration complexity)
