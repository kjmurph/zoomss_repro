# ZooMSS Fishing Implementation: Revised Assessment and Implementation Roadmap

**Prepared by:** Claude (Opus 4.5)  
**For:** Kieran Murphy  
**Date:** January 2026  
**Purpose:** Pre-implementation analysis with modular two-repository architecture

---

## Architecture Overview

Following the DBPM framework, the implementation is split across two repositories:

### Repository 1: `zoomss` (R package)
**Purpose:** Core ZooMSS model with generic fishing mortality functionality  
**Analogy:** Like `mizer` — standalone package for size spectrum modelling with fishing

Adds to the existing package:
- Fishing mortality functions (selectivity, F calculation, catch calculation)
- Integration of fishing mortality into model dynamics
- Generic interfaces for effort, catchability, selectivity parameters
- Theoretical/scenario-based fishing applications

### Repository 2: `zoomss-isimip3a` (workflow scripts)
**Purpose:** ISIMIP3a protocol implementation using ZooMSS  
**Template:** `lme_scale_calibration_ISMIP3a` (DBPM repo)

Contains:
- Data processing scripts (climate, fishing effort, catch)
- Regional calibration workflows
- Gridded model execution
- Output processing and validation
- Protocol-specific configurations

---

## Part 1: Critical Assessment (Revised)

### 1.1 Issues Affecting `zoomss` Package Design

These must be resolved before package implementation:

#### A. Fishing Mortality Integration Point (RESOLVED ✓)
**Finding:** Fishing mortality already exists in ZooMSS and the integration point is clear.

**Current implementation (from code review):**

1. **Groups data frame** contains fishing parameters:
   - `Fmort`: Fishing mortality rate (yr⁻¹)
   - `Fmort_W0`: Minimum vulnerable size (log10 g)
   - `Fmort_Wmax`: Maximum vulnerable size (log10 g)

2. **`zoomss_setup.R` (lines 130-133)** applies static fishing mortality:
   ```r
   for(g in 1:param$ngrps){
     model$fish_mort[g, match(param$Groups$Fmort_W0[g], param$w_log10):
                        match(param$Groups$Fmort_Wmax[g], param$w_log10)] <- param$Groups$Fmort[g]
   }
   ```

3. **`zoomss_run.R` (line 180)** combines mortalities:
   ```r
   Z <- sweep(model$M_sb + model$fish_mort, 2, M2, '+')
   ```
   Where: `Z = senescence + fishing + predation`

**What this means for implementation:**
- Fishing mortality is **already knife-edge selectivity** (between Fmort_W0 and Fmort_Wmax)
- The architecture supports fishing; we just need to make it **dynamic over time**
- Follow the pattern used for temperature effects: pre-calculate in `zoomss_params`, update each timestep in `zoomss_run`

**Recommended approach:**
```r
# In zoomss_run, inside time loop (similar to temp_eff update):
if(dynamic_fishing) {
  current_effort <- effort_ts[itime]
  for(g in fish_grps) {
    model$fish_mort[g, selected_sizes] <- current_effort * q[g] * selectivity[g, ]
  }
}
```

---

#### B. Size Class Structure (RESOLVED ✓)
**Finding:** Size classes are log10-spaced with 0.1 increments.

From `zoomss_params.R`:
```r
param$dx = 0.1  # log10 weight step
param2$w_log10 <- round(seq(from = min(Groups$W0), to = max(Groups$Wmax), param$dx), digits = 2)
param2$w <- 10^(seq(from = min(Groups$W0), to = max(Groups$Wmax), param$dx))
```

**Key values:**
- `dx = 0.1` (log10 step size)
- `w`: Body mass in grams
- `w_log10`: Log10 of body mass
- `ngrid`: Number of size classes

This means selectivity thresholds (w_min) should be specified in **log10 grams** to match existing conventions (Fmort_W0, Fmort_Wmax).

---

#### C. Environmental Forcing Pattern (RESOLVED ✓)
**Finding:** Dynamic forcing uses pre-calculated time series arrays.

From `zoomss_params.R` and `zoomss_run.R`:
1. Environmental data (SST, Chl) provided as columns in `input_params`
2. Derived parameters (phyto_int, phyto_slope, temp_eff) pre-calculated for all timesteps
3. In run loop, current timestep values accessed by index: `temp_eff_zoo[itime, i]`

**Pattern to follow for fishing effort:**
```r
# In zoomss_params:
param2$effort_ts <- input_params$effort  # Add effort time series
param2$fishing_q <- fishing_params$q      # Catchability per group
param2$fishing_wmin <- fishing_params$w_min  # Size threshold per group

# In zoomss_run:
current_effort <- param$effort_ts[itime]
# Update fish_mort based on current_effort and parameters
```

---

#### B. Selectivity Function Flexibility (RESOLVED ✓)
**Finding:** Current implementation is knife-edge. Package should support multiple types for flexibility.

**Current implementation:** `Fmort_W0` to `Fmort_Wmax` defines vulnerable size range (knife-edge).

**Recommendation for `zoomss` package:**
- Implement selectivity as a **generic interface** (Task A.2)
- Provide built-in options: `knife_edge`, `logistic`
- Allow user-defined selectivity functions
- Default to knife-edge for backward compatibility
- ISIMIP3a workflow uses knife-edge (as specified in implementation guide)

---

#### C. Effort vs. Direct F Specification (RESOLVED ✓)
**Finding:** Current implementation uses direct F specification via `Groups$Fmort`.

**Recommendation for `zoomss` package:**
- Support both modes:
  1. **Static F (current):** Use `Groups$Fmort`, `Fmort_W0`, `Fmort_Wmax`
  2. **Dynamic effort-based:** Add `effort` column to `input_params` + fishing parameters
- Backward compatible: if no `fishing` argument provided, use current static behaviour
- ISIMIP3a workflow uses effort-based approach

---

#### D. Size Class Structure (RESOLVED ✓)
**Finding:** Log10-spaced with dx = 0.1 increments.

**Key details from code:**
- `w`: Body mass in grams
- `w_log10`: Log10(body mass), rounded to 2 decimal places
- `dx = 0.1`: Step size in log10 space
- Range: `min(Groups$W0)` to `max(Groups$Wmax)`

**Implication:** All size parameters (w_min, Fmort_W0, etc.) should be in **log10 grams** to match existing conventions.

---

### 1.2 Issues Affecting `zoomss-isimip3a` Workflow

These are specific to the protocol implementation:

#### A. Data Format Alignment with DBPM Repo
**Issue:** The DBPM repo uses zarr for climate data. Does ZooMSS currently expect NetCDF, CSV, or other formats?

**Recommendation:** Adopt the DBPM approach:
- Steps 1-2: Python scripts process GFDL-MOM6-COBALT2 to zarr
- Step 3: R script loads zarr (or converted format) and fishing data
- This maximises code reuse from the DBPM repo

---

#### B. Regional Unit: LME vs FAO
**Issue:** Your implementation guide references LMEs. The DBPM repo uses FAO Major Fishing Areas.

**Decision required:** Which regional unit to use?
- **LME:** 66 regions, coastal focus, commonly used in FishMIP
- **FAO:** 19 major areas, includes open ocean, used in DBPM repo

**Recommendation:** Use FAO areas to align with DBPM repo structure, which will ease code adaptation. Document mapping to LMEs if needed for comparison.

---

#### C. Functional Group Mapping
**Issue:** DBPM has different functional groups than ZooMSS.

| DBPM | ZooMSS |
|------|--------|
| Benthic detritivores | ? |
| Pelagic predators | Small / Medium / Large fish |

**Required:** Define explicit mapping from FishMIP 6 groups → ZooMSS 3 groups (as per your implementation guide Section 3.1).

---

### 1.3 Numerical and Computational Concerns

(These remain from original assessment, now clearly assigned to appropriate repo)

| Concern | Affects | Mitigation |
|---------|---------|------------|
| Parameter identifiability (q vs w_min) | `zoomss-isimip3a` | Regularisation, fix w_min from literature |
| Objective function scaling | `zoomss-isimip3a` | Use log-transformed catch |
| Knife-edge discontinuity | `zoomss-isimip3a` | Use Nelder-Mead optimizer |
| Initial conditions / spin-up | `zoomss-isimip3a` | Include 1841-1959 spin-up period (per DBPM) |
| Negative biomass safeguards | `zoomss` | Floor at zero, log warnings |

---

## Part 2: Decision Log

### Decisions for `zoomss` Package (CONFIRMED)

| ID | Decision | Resolution | Rationale |
|----|----------|------------|-----------|
| P1 | Selectivity interface | **Multiple types** (knife_edge, logistic) | Flexibility for different use cases |
| P2 | Fishing input mode | **Both** static and dynamic | Backward compatibility + new functionality |
| P3 | Default selectivity | **Knife-edge** | Matches current implementation |
| P4 | Fishing output | **F + Catch** | Comprehensive but not excessive |
| P5 | Size parameter units | **Log10 grams** | Matches existing ZooMSS conventions |

### Decisions for `zoomss-isimip3a` Workflow (CONFIRMED)

| ID | Decision | Resolution | Rationale |
|----|----------|------------|-----------|
| W1 | Regional unit | **FAO Major Areas** | Aligns with DBPM repo (confirmed) |
| W2 | Effort aggregation | **Sum** | Additive fishing pressure |
| W3 | Objective function | **Log-SSE** | Handles magnitude differences |
| W4 | Optimizer | **Nelder-Mead** | Robust to discontinuities |
| W5 | Spin-up period | **1841-1959** | Per DBPM approach |
| W6 | Climate data format | **zarr** | Reuse DBPM processing |
| W7 | Model execution | **Adapt DBPM pattern** | May need R↔Python bridge |

---

## Part 3: Implementation Roadmap

### Phase A: `zoomss` Package — Fishing Module

These tasks extend the existing `zoomss` R package with dynamic fishing functionality.

**Key finding from code review:** ZooMSS already has static fishing mortality. The goal is to make it **dynamic** (time-varying based on effort) while maintaining backward compatibility.

---

#### Task A.1: Document existing ZooMSS mortality structure (COMPLETED)

**Summary of findings:**

| Component | Location | Description |
|-----------|----------|-------------|
| Fishing parameters | `Groups` data frame | `Fmort`, `Fmort_W0`, `Fmort_Wmax` |
| Static F setup | `zoomss_setup.R:130-133` | Applies constant F from Groups |
| F storage | `model$fish_mort` | Matrix (ngrps × ngrid) |
| Mortality combination | `zoomss_run.R:180` | `Z = M_sb + fish_mort + M2` |
| Size classes | `zoomss_params.R` | Log10-spaced, dx=0.1, stored in `w` and `w_log10` |

**Existing selectivity:** Already knife-edge (F applied between Fmort_W0 and Fmort_Wmax).

---

#### Task A.2: Implement selectivity functions
**File:** `R/selectivity.R` (new file — justified as distinct functional module)

```r
#' Calculate size selectivity
#'
#' @param w Numeric vector of body masses (g) OR log10 body masses
#' @param type Selectivity type: "knife_edge", "logistic", or "custom"
#' @param params List of parameters:
#'   - knife_edge: w_min (minimum vulnerable size)
#'   - logistic: L50, L95 (sizes at 50% and 95% selectivity)
#' @param log_scale Logical, TRUE if w and params are in log10 scale (default TRUE to match ZooMSS convention)
#' @return Numeric vector of selectivity values (0-1)
#' @export
calc_selectivity <- function(w, type = "knife_edge", params = list(), log_scale = TRUE) {
  # Implementation
}

#' @keywords internal
selectivity_knife_edge <- function(w, w_min) {
  as.numeric(w >= w_min)
}

#' @keywords internal
selectivity_logistic <- function(w, L50, L95) {
  1 / (1 + exp(-log(19) * (w - L50) / (L95 - L50)))
}
```

**Tests:**
```r
# Knife-edge (log10 scale, matching ZooMSS convention)
calc_selectivity(c(-1, 0, 1, 2), type = "knife_edge", params = list(w_min = 0))
# Expected: c(0, 1, 1, 1)

# Logistic
calc_selectivity(c(0, 1, 2, 3), type = "logistic", params = list(L50 = 1.5, L95 = 2.5))
# Expected: values between 0 and 1, ~0.5 at w=1.5
```

---

#### Task A.3: Implement fishing mortality calculation
**File:** `R/fishing.R` (new file)

```r
#' Calculate fishing mortality from effort
#'
#' @param effort Numeric, fishing effort at time t
#' @param q Numeric vector, catchability coefficient per group
#' @param selectivity Matrix (ngrps × ngrid), selectivity at each size for each group
#' @return Matrix (ngrps × ngrid), F(w) for each group and size class
#' @export
calc_fishing_mortality <- function(effort, q, selectivity) {
  # F(w,g) = effort * q[g] * selectivity[g, w]
  sweep(selectivity, 1, effort * q, "*")
}
```

---

#### Task A.4: Implement catch calculation
**File:** `R/fishing.R` (add to same file)

```r
#' Calculate catch from fishing mortality and abundance
#'
#' @param F_w Matrix (ngrps × ngrid), fishing mortality by group and size
#' @param N_w Matrix (ngrps × ngrid), abundance by group and size
#' @param w Numeric vector, body mass of each size class (g)
#' @param dw Numeric, width of size classes in log10 space (default 0.1)
#' @return Numeric vector, catch per group (g), or total if sum = TRUE
#' @export
calc_catch <- function(F_w, N_w, w, dw = 0.1) {
  # Catch per group = sum over sizes of F * N * w * dw
  # Note: dw in log space, so actual width = w * log(10) * dw
  rowSums(F_w * N_w * w) * log(10) * dw
}
```

---

#### Task A.5: Add dynamic fishing to ZooMSS
**Files to modify:** `R/zoomss_params.R`, `R/zoomss_run.R`, `R/zoomss_model.R`

**A.5.1: Extend `input_params` to accept effort time series**

In `zoomss_params.R`, add after line 151 (after temp_eff calculation):

```r
# Process fishing effort time series if provided
if ("effort" %in% names(input_params)) {
  param2$effort_ts <- input_params$effort
  param2$dynamic_fishing <- TRUE
} else {
  param2$dynamic_fishing <- FALSE
}
```

**A.5.2: Add fishing parameters argument to `zoomss_model`**

```r
zoomss_model <- function(input_params, Groups = NULL, isave = 1, fishing = NULL) {
  # fishing = NULL: use static Fmort from Groups (current behaviour)
  # fishing = list(q = c(...), w_min = c(...), w_max = c(...)):
  #   use effort from input_params with these parameters
}
```

**A.5.3: Modify `zoomss_setup.R` to handle dynamic fishing**

Replace lines 130-133 with:

```r
# Fishing mortality setup
if (param$dynamic_fishing && !is.null(param$fishing_params)) {
  # Dynamic fishing: initialize fish_mort as zeros, will be updated each timestep
  model$fish_mort <- matrix(0, nrow = param$ngrps, ncol = param$ngrid)
  
  # Pre-calculate selectivity matrix (constant over time)
  model$fishing_selectivity <- matrix(0, nrow = param$ngrps, ncol = param$ngrid)
  for (g in param$fish_grps) {
    w_min <- param$fishing_params$w_min[g]
    w_max <- param$fishing_params$w_max[g]
    selected <- param$w_log10 >= w_min & param$w_log10 <= w_max
    model$fishing_selectivity[g, selected] <- 1
  }
  model$fishing_q <- param$fishing_params$q
  
} else {
  # Static fishing (original behaviour)
  for(g in 1:param$ngrps){
    model$fish_mort[g, match(param$Groups$Fmort_W0[g], param$w_log10):
                       match(param$Groups$Fmort_Wmax[g], param$w_log10)] <- param$Groups$Fmort[g]
  }
}
```

**A.5.4: Modify `zoomss_run.R` to update fishing mortality each timestep**

Add after line 159 (after M_sb update), inside the time loop:

```r
# Update fishing mortality with current effort (if dynamic fishing enabled)
if (param$dynamic_fishing) {
  current_effort <- param$effort_ts[itime]
  model$fish_mort <- calc_fishing_mortality(
    effort = current_effort,
    q = model$fishing_q,
    selectivity = model$fishing_selectivity
  )
}
```

---

#### Task A.6: Add fishing outputs to model results
**File:** Modify `zoomss_run.R`

Add to the output saving block (after line 261):

```r
# Save fishing outputs if dynamic fishing
if (param$dynamic_fishing) {
  model$F[isav,,] <- model$fish_mort  # Fishing mortality
  model$catch[isav,] <- calc_catch(model$fish_mort, N, param$w, param$dx)  # Catch by group
}
```

Add array initialization in `zoomss_setup.R`:

```r
if (param$dynamic_fishing) {
  model$F <- array(NA, dim = c(param$nsave, param$ngrps, param$ngrid))
  model$catch <- array(NA, dim = c(param$nsave, param$ngrps))
}
```

---

#### Task A.7: Create fishing vignette
**File:** `vignettes/fishing.Rmd`

**Content outline:**
1. Introduction to fishing in ZooMSS
2. Static fishing (current approach using Groups$Fmort)
3. Dynamic fishing with effort time series
4. Example: comparing fished vs unfished scenarios
5. Example: time-varying effort simulation
6. Selectivity options (knife-edge vs logistic)

---

### Phase B: `zoomss-isimip3a` Repository Setup

Create a new repository mirroring the DBPM structure.

**CRITICAL:** For each task in Phase B, fetch the corresponding DBPM script from  
`https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a/blob/new_features/scripts/`  
and use it as the implementation template.

---

#### Task B.1: Initialise repository structure
**Based on:** `lme_scale_calibration_ISMIP3a`

```
zoomss-isimip3a/
├── scripts/
│   ├── 01_processing_zoomss_global_inputs.ipynb
│   ├── 02_processing_zoomss_regional_inputs.py
│   ├── 03_processing_effort_fishing_inputs.R
│   ├── 04_calibrating_zoomss_fishing_params.R
│   ├── 05_setup_gridded_zoomss.py
│   ├── 06_run_zoomss_gridded.R  # or .py depending on ZooMSS interface
│   ├── 07_merging_outputs.py
│   ├── 08_calculating_catches_zoomss.R
│   └── 09_plotting_outputs.Rmd
├── data/
│   └── (processed regional inputs)
├── outputs/
│   └── (model outputs by region)
├── docs/
│   ├── ZooMSS_Fishing_Implementation_Guide.md
│   └── ZooMSS_Fishing_Implementation_Assessment.md
├── README.md
├── LICENSE
└── zoomss-isimip3a.Rproj
```

---

#### Task B.2: Adapt Step 1 — Global climate input processing
**Fetch first:** `scripts/01_processing_dbpm_global_inputs.ipynb`  
**Adaptation:** Minimal — same GFDL-MOM6-COBALT2 processing

**Changes needed:**
- Update variable selection if ZooMSS requires different inputs than DBPM
- Verify SST and Chl extraction matches ZooMSS requirements

---

#### Task B.3: Adapt Step 2 — Regional climate input processing
**Fetch first:** `scripts/02_processing_dbpm_regional_inputs.py`  
**Adaptation:** Change region definitions if using different areas

---

#### Task B.4: Adapt Step 3 — Fishing effort processing
**Fetch first:** `scripts/03_processing_effort_fishing_inputs.R`  
**Adaptation required:**

- Update functional group mapping:
  ```r
  # DBPM mapping (for reference)
  # Benthic detritivores, Pelagic predators
  
  # ZooMSS mapping (to implement)
  fishmip_to_zoomss <- list(
    Small = c("Small pelagics", "Bathypelagics"),
    Medium = c("Medium pelagics", "Benthopelagics"),
    Large = c("Large pelagics", "Demersals")
  )
  ```
- Aggregate effort by summing within ZooMSS groups
- Output format compatible with ZooMSS

---

#### Task B.5: Create Step 4 — ZooMSS calibration script
**Fetch first:** `scripts/04_calculating_dbpm_fishing_params.R`  
**Major adaptation required — study DBPM implementation carefully**

This is the core calibration script. Understand how DBPM structures:
- Parameter definitions and bounds
- Objective function construction
- Optimisation loop
- Diagnostic output generation

Then adapt for ZooMSS's 6-parameter structure (q and w_min for Small/Medium/Large).

Structure:

```r
# 04_calibrating_zoomss_fishing_params.R

library(zoomss)
library(optimx)  # or optim

# 1. Load regional data (climate + fishing + observed catch)
# 2. Define objective function
objective_function <- function(params, data, zoomss_args) {
  # Unpack params: q_small, w_min_small, q_medium, w_min_medium, q_large, w_min_large
  # Run ZooMSS with fishing
  # Calculate predicted catch
  # Return log-SSE against observed catch
}

# 3. Set parameter bounds
bounds <- list(
  lower = c(q_small = 0.001, w_min_small = 1, ...),
  upper = c(q_small = 0.05, w_min_small = 10, ...)
)

# 4. Run optimisation (per region)
results <- optim(
  par = starting_values,
  fn = objective_function,
  method = "Nelder-Mead",
  ...
)

# 5. Save calibrated parameters
# 6. Generate diagnostic plots
```

---

#### Task B.6: Adapt Steps 5-6 — Gridded model setup and execution
**Fetch first:** `scripts/05_setup_gridded_DBPM.py` and `scripts/06_running_gridded_DBPM.py`

**Decision needed:** Is ZooMSS run from R or Python?
- If R: Convert these to R scripts calling `zoomss` package
- If Python: Ensure ZooMSS has Python bindings or use rpy2

---

#### Task B.7: Adapt Steps 7-9 — Output processing and plotting
**Fetch first:** `scripts/07_calculating_catches_DBPM.py` and `scripts/08_plotting_gridded_DBPM_outputs.ipynb`  
**Adaptation:** Update for ZooMSS output structure

---

## Part 4: Suggested Implementation Order

### Immediate (Weeks 1-2): `zoomss` Package Core

1. **Task A.1** — Document existing mortality structure
2. **Task A.2** — Implement selectivity functions
3. **Task A.3** — Implement fishing mortality calculation
4. **Task A.4** — Implement catch calculation
5. **Task A.5** — Integrate into ZooMSS dynamics (this is the critical integration)

### Short-term (Weeks 3-4): `zoomss` Package Polish

6. **Task A.6** — Add fishing outputs
7. **Task A.7** — Create vignette with examples
8. Run `devtools::check()`, ensure CRAN compatibility

### Medium-term (Weeks 5-8): `zoomss-isimip3a` Workflow

9. **Task B.1** — Initialise repo structure
10. **Task B.2-B.3** — Adapt climate processing (likely minimal changes)
11. **Task B.4** — Adapt fishing effort processing
12. **Task B.5** — Create calibration script (most complex adaptation)
13. **Test calibration** on single FAO region before scaling

### Later (Weeks 9+): Full Implementation

14. **Task B.6** — Gridded model execution
15. **Task B.7** — Output processing
16. Run full calibration across all regions
17. Validation and documentation

---

## Part 5: Instructions for Copilot (Sonnet 4.5)

### General Workflow

Before implementing any task:
1. Read the task specification in this document
2. If it's a Phase B task, **fetch the corresponding DBPM script** from the reference repo
3. Study the DBPM implementation to understand patterns and data structures
4. Plan your approach and share it before coding
5. Implement, following existing code conventions

### DBPM Reference Repository

**URL:** `https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a`  
**Branch:** `new_features`

Fetch scripts from: `scripts/` directory

### For `zoomss` package work (Phase A):

```markdown
## Context
You are adding fishing mortality functionality to the `zoomss` R package.
This should be generic/reusable (like mizer), not protocol-specific.

## Key files to create
- `R/selectivity.R` — Selectivity functions
- `R/fishing.R` — Fishing mortality and catch calculations

## Design principles
- Support multiple selectivity types (knife_edge, logistic, custom)
- Support both direct F specification and effort-based calculation
- Maintain backward compatibility (fishing = NULL runs without fishing)
- Follow existing zoomss code style and conventions
```

### For `zoomss-isimip3a` workflow work (Phase B):

```markdown
## Context
You are adapting the DBPM ISIMIP3a workflow for ZooMSS.

## CRITICAL: Fetch-first approach
Before implementing ANY Phase B task:
1. Fetch the corresponding DBPM script from:
   https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a/blob/new_features/scripts/
2. Read and understand the full script
3. Document what needs to change for ZooMSS vs what can be reused
4. Then implement

## Key adaptations
- Functional group mapping: FishMIP 6 groups → ZooMSS 3 groups (Small/Medium/Large)
- Model calls: Replace DBPM functions with zoomss equivalents
- Parameters: ZooMSS uses q + w_min per group (6 total), DBPM structure differs

## What to preserve from DBPM
- Overall workflow structure (Steps 1-8)
- Data loading and formatting patterns
- Output file structures
- Plotting conventions
- NCI/Gadi infrastructure setup

## What requires adaptation
- Functional group definitions and mapping
- Calibration objective function (Step 4)
- Model execution interface (Step 6)
- Catch calculation specifics (Step 7)
```

---

## Part 6: Status and Next Steps

### Resolved Questions

| Question | Resolution |
|----------|------------|
| Regional unit | **FAO Major Areas** (confirmed) |
| Mortality implementation | **Reviewed** — clear integration point at `zoomss_run.R:180` |
| Size class structure | **Reviewed** — log10-spaced, dx=0.1, units in log10 grams |
| Fishing parameters | **Existing** — Groups has Fmort, Fmort_W0, Fmort_Wmax |

### Remaining Questions

1. **NCI/Gadi access:** Will this run on Gadi? (Affects infrastructure setup reuse from DBPM)

2. **Model execution for gridded runs:** DBPM uses Python. ZooMSS is R-based. Options:
   - Keep R (adapt DBPM Python scripts to R)
   - Use rpy2 to call ZooMSS from Python
   - Other approach?

3. **Timeline priority:** Package-first or parallel development?

### Suggested First Steps

1. **Start with `zoomss` package (Phase A)**
   - Task A.2: Implement `R/selectivity.R`
   - Task A.3-A.4: Implement `R/fishing.R`
   - Task A.5: Modify zoomss_params, zoomss_setup, zoomss_run
   - Test with synthetic effort time series before moving to calibration

2. **Then `zoomss-isimip3a` (Phase B)**
   - Task B.1: Initialize repo structure
   - Have Sonnet fetch DBPM scripts and adapt

---

**Document Version:** 2.1  
**Status:** Ready for implementation — package design confirmed
