# ZooMSS Fishing Implementation Guide
## Comprehensive Technical Summary for Historical Fishing Mortality Implementation

**Author:** Kieran Murphy  
**Date:** December 2024  
**Context:** ISIMIP3a Protocol - Historical Fishing Implementation (1960-2010)  
**Model:** Zooplankton Model of Size Spectra (ZooMSS)

---

## Executive Summary

This document provides a complete technical specification for implementing historical fishing mortality in ZooMSS, a global marine ecosystem model. The approach uses a catchability-based fishing mortality framework calibrated to observed catch data from FishMIP ISIMIP3a, with regional calibration at the Large Marine Ecosystem (LME) scale.

---

## 1. Overview of Approach

### 1.1 Model Context
- **Model:** ZooMSS (Zooplankton Model of Size Spectra)
- **Type:** Global marine ecosystem model with size-structured food web representation
- **Resolution:** 1° and 0.25° global grids
- **Fish Groups:** Three size-based fish groups (Small, Medium, Large)
- **Time Period:** 1960-2010 (calibration and historical simulations)

### 1.2 Fishing Implementation Philosophy
- **Approach:** Catchability-based fishing mortality with knife-edge selectivity
- **Calibration Scale:** Large Marine Ecosystem (LME) or regional level
- **Parameters per Group:** 2 parameters × 3 fish groups = 6 total parameters
  - Catchability coefficient (q)
  - Size threshold for vulnerability (w_min)
- **Data Source:** FishMIP ISIMIP3a fishing effort data (Rousseau et al. 2024, Nature Scientific Data)
- **Environmental Forcing:** SST(x,y,t) and Chl(x,y,t) from GFDL-MOM6-COBALT2

---

## 2. Fishing Mortality Equation

### 2.1 Core Equation
The fishing mortality rate at size w and time t is given by:

```
F(w,t) = Effort(t) × q × Selectivity(w)
```

Where:
- **F(w,t):** Fishing mortality rate (per time unit) for organisms of body mass w at time t
- **Effort(t):** Observed fishing effort at time t from FishMIP data (units: kW-days or equivalent)
- **q:** Catchability coefficient (dimensionless, to be calibrated)
- **Selectivity(w):** Size selectivity function (dimensionless, 0 or 1)

### 2.2 Selectivity Function (Knife-Edge)
A simplified knife-edge selectivity function is used:

```
Selectivity(w) = {
    0,  if w < w_min
    1,  if w ≥ w_min
}
```

**Biological Interpretation:**
- Fish below w_min are not vulnerable to fishing (selectivity = 0)
- Fish at or above w_min are fully vulnerable to fishing (selectivity = 1)
- This represents a sharp transition at the threshold size
- Simpler than logistic selectivity curves but computationally efficient

### 2.3 Predicted Catch Calculation
The predicted catch at time t is calculated by integrating fishing mortality over the size spectrum:

```
Catch(t) = ∫[w_min to ∞] F(w,t) × N(w,t) × w dw
```

Where:
- **N(w,t):** Abundance (number of individuals) at size w and time t
- **w:** Body mass (individual weight)
- The integration is performed over all vulnerable sizes (w ≥ w_min)

---

## 3. Functional Group Mapping

### 3.1 FishMIP to ZooMSS Mapping
FishMIP provides fishing effort data for 6 functional groups, which must be mapped to ZooMSS's 3 size-based fish groups:

**FishMIP Functional Groups → ZooMSS Fish Groups:**

| FishMIP Groups | ZooMSS Group | Size Range |
|----------------|--------------|------------|
| Small pelagics, Bathypelagics | **Small** | Smallest fish |
| Medium pelagics, Benthopelagics | **Medium** | Medium-sized fish |
| Large pelagics, Demersals | **Large** | Largest fish |

**Implementation Note:**
- Each ZooMSS fish group receives its own set of parameters (q, w_min)
- Fishing effort from multiple FishMIP groups may need to be aggregated or weighted
- The mapping ensures that pelagic and demersal components are represented across size classes

---

## 4. Calibration Workflow

### 4.1 Overview
The calibration is an iterative optimization process that estimates catchability (q) and size thresholds (w_min) for each fish group by minimizing the difference between predicted and observed catch.

### 4.2 Four-Step Calibration Loop

#### Step 1: Map FishMIP Functional Groups to Fish Groups in ZooMSS
- Establish correspondence between FishMIP's 6 groups and ZooMSS's 3 size-based groups
- Aggregate or weight fishing effort data as needed
- This mapping is fixed throughout calibration

#### Step 2: Calculate Fishing Mortality F(w,t)
- For each fish group and time step:
  - Apply current parameter estimates (q, w_min)
  - Calculate F(w,t) = Effort(t) × q × Selectivity(w)
- Generate size- and time-resolved fishing mortality fields

#### Step 3: Run Dynamic ZooMSS with Fishing Mortality
- Integrate F(w,t) into ZooMSS's size-spectrum dynamics
- Run the model forward in time with:
  - Environmental forcing: SST(x,y,t), Chl(x,y,t) from GFDL-MOM6-COBALT2
  - Fishing mortality: F(w,t) calculated in Step 2
- Generate predictions of:
  - Size distributions N(w,t)
  - Predicted catch from each fish group

#### Step 4: Optimize Parameters
- **Objective:** Minimize sum of squared errors (SSE) between predicted and observed catch
- **Objective Function:**
  ```
  SSE = Σ[t=1960 to 2010] (Catch_predicted(t) - Catch_observed(t))²
  ```
- **Optimization Method:** e.g., L-BFGS-B (or similar bounded optimization)
- **Parameters:** Adjust q and w_min for each fish group within bounds
- **Iteration:** Return to Step 2 with updated parameters until convergence

### 4.3 Iteration Scope
The **iterate arrow** in the calibration loop connects:
- **From:** Step 4 (Optimize Parameters)
- **To:** Step 2 (Calculate Fishing Mortality F(w,t))
- **Not included in iteration:** Step 1 (functional group mapping is static)

---

## 5. Parameter Specifications

### 5.1 Parameters to Estimate (6 Total)

Each of the 3 fish groups (Small, Medium, Large) requires 2 parameters:

#### Catchability (q)
- **Definition:** Proportionality constant relating fishing effort to fishing mortality
- **Units:** Dimensionless
- **Interpretation:** Higher q means more effective fishing (higher catch per unit effort)

#### Size Threshold (w_min)
- **Definition:** Minimum body mass at which fish become vulnerable to fishing
- **Units:** Grams (g)
- **Interpretation:** Represents gear selectivity and minimum catchable size

### 5.2 Parameter Bounds

Parameter bounds are based on biological realism and numerical stability:

| Fish Group | Catchability (q) | Size Threshold (w_min) |
|------------|------------------|------------------------|
| **Small Fish** | 0.001 - 0.05 | 1 - 10 g |
| **Medium Fish** | 0.005 - 0.08 | 10 - 100 g |
| **Large Fish** | 0.01 - 0.1 | 200 - 500 g |

**Rationale for Catchability Bounds:**
- **Large fish have higher catchability:** Larger fish are bigger targets, easier to locate with fishing gear, and fishing operations typically target larger individuals
- **Small fish have lower catchability:** Smaller fish can escape through nets more easily and are less targeted by commercial fisheries
- **Progressive increase:** q increases with fish size (Small < Medium < Large)

**Rationale for Size Threshold Bounds:**
- Reflect typical gear selectivity and market preferences
- Larger fish groups have higher w_min (larger minimum catchable size)
- Bounds prevent biologically unrealistic selectivity patterns

---

## 6. Optimization Framework

### 6.1 Objective Function
Minimize the sum of squared errors (SSE) between predicted and observed catch:

```
Minimize: SSE = Σ[t=1960 to 2010] Σ[groups] (Catch_predicted(group,t) - Catch_observed(group,t))²
```

**Considerations:**
- May include weighting by fish group if data quality varies
- May use relative errors or log-transformed catch to handle large magnitude differences
- Temporal weighting may be applied if recent data is more reliable

### 6.2 Optimization Method
**Recommended:** L-BFGS-B (Limited-memory Broyden-Fletcher-Goldfarb-Shanno with bounds)

**Characteristics:**
- Quasi-Newton method suitable for bounded optimization
- Handles ~6 parameters efficiently
- Gradient-based, so requires smooth objective function
- Well-suited for ecological parameter estimation

**Alternative Methods:**
- Particle Swarm Optimization (PSO) - gradient-free, more robust to local minima
- Differential Evolution - evolutionary algorithm, good for complex landscapes
- Bayesian optimization - if uncertainty quantification is needed

### 6.3 Calibration Period
- **Period:** 1960-2010
- **Justification:** FishMIP protocol specifies this historical period
- **Data Availability:** Fishing effort and catch observations available for this period

### 6.4 Regional Calibration Strategy
**Key Principle:** Calibrate separately for each Large Marine Ecosystem (LME) or region with distinct fishery characteristics

**Rationale:**
- Different regions have different fishing practices, gear types, and target species
- Fleet composition varies geographically
- Environmental conditions affect catchability
- Regional calibration captures local fishery dynamics better than global parameters

**Implementation:**
- Run optimization independently for each LME
- Parameters (q, w_min) are LME-specific
- Computational cost scales with number of LMEs
- May require hundreds to thousands of model runs per LME

---

## 7. Input Data Requirements

### 7.1 Fishing Effort Data
- **Source:** FishMIP ISIMIP3a Historical Fishing Effort Dataset
- **Reference:** Rousseau et al. 2024, Nature Scientific Data
- **Format:** Gridded (spatial) and time-resolved
- **Functional Groups:** 6 groups (small pelagics, medium pelagics, large pelagics, bathypelagics, benthopelagics, demersals)
- **Units:** Typically kW-days or equivalent fishing effort metric
- **Temporal Resolution:** Annual or monthly
- **Spatial Resolution:** Should match or be interpolated to ZooMSS grid resolution

### 7.2 Environmental Forcing
- **Source:** GFDL-MOM6-COBALT2 Earth System Model
- **Variables:**
  - **SST(x,y,t):** Sea surface temperature
  - **Chl(x,y,t):** Chlorophyll-a concentration (proxy for primary productivity)
- **Purpose:** Drive ZooMSS ecosystem dynamics (growth, reproduction, predation)
- **Period:** 1960-2010 (matching calibration period)

### 7.3 Observed Catch Data
- **Source:** FishMIP ISIMIP3a protocol (likely FAO catch statistics or similar)
- **Purpose:** Target data for calibration (observed catch to match)
- **Resolution:** Annual totals by functional group and region/LME
- **Units:** Tonnes or kg per year
- **Quality:** Critical for calibration success; data gaps or errors will affect parameter estimates

---

## 8. Model Integration Details

### 8.1 ZooMSS Model Structure
- **Framework:** Size-spectrum model representing marine food webs
- **Size Classes:** Continuous size spectrum or discrete size bins
- **Functional Groups:** Multiple groups including zooplankton and fish
- **Dynamics:** Growth, predation, natural mortality, and (with this implementation) fishing mortality

### 8.2 Fishing Mortality Integration
Fishing mortality F(w,t) must be integrated into ZooMSS's existing mortality terms:

```
dN(w,t)/dt = [Growth] + [Reproduction] - [Predation] - [Natural Mortality] - [Fishing Mortality]
```

Where the fishing mortality term is:
```
Fishing Loss(w,t) = F(w,t) × N(w,t)
```

**Implementation Notes:**
- F(w,t) is calculated externally (Step 2) and passed to ZooMSS
- Fishing mortality is applied additively to existing mortality sources
- Size-specific application ensures only vulnerable sizes experience fishing

### 8.3 Catch Calculation in Post-Processing
After each model time step, calculate predicted catch:

```
Catch(t) = ∫ F(w,t) × N(w,t) × w dw
```

**Numerical Implementation:**
- For discrete size bins: `Catch(t) = Σ[i] F(w_i,t) × N(w_i,t) × w_i × Δw_i`
- Accumulate catch over each time step (e.g., annual totals)
- Store for comparison with observed catch in optimization

---

## 9. Implementation Steps Summary

### Phase 1: Data Preparation
1. Obtain FishMIP ISIMIP3a fishing effort data (1960-2010)
2. Obtain GFDL-MOM6-COBALT2 environmental forcing (SST, Chl)
3. Obtain observed catch data by functional group and region
4. Define LME boundaries or regional domains for calibration
5. Map FishMIP functional groups to ZooMSS fish groups

### Phase 2: Model Modification
1. Implement fishing mortality calculation: F(w,t) = Effort(t) × q × Selectivity(w)
2. Implement knife-edge selectivity function based on w_min
3. Integrate F(w,t) into ZooMSS mortality equations
4. Implement catch calculation: Catch(t) = ∫ F(w,t) × N(w,t) × w dw
5. Set up parameter input/output for q and w_min (6 parameters)

### Phase 3: Calibration Setup
1. Define parameter bounds for each fish group
2. Set up objective function (SSE between predicted and observed catch)
3. Choose optimization algorithm (e.g., L-BFGS-B)
4. Define convergence criteria
5. Set up regional/LME-specific calibration loops

### Phase 4: Calibration Execution
1. For each LME/region:
   a. Initialize parameters (q, w_min) - use reasonable starting values
   b. Run optimization loop:
      - Calculate F(w,t) with current parameters
      - Run ZooMSS with fishing mortality
      - Calculate predicted catch
      - Compute SSE against observed catch
      - Update parameters via optimizer
      - Repeat until convergence
   c. Store optimized parameters for this LME
2. Validate calibrated parameters (check biological realism, temporal patterns)

### Phase 5: Historical Simulations
1. Apply calibrated parameters (q, w_min) by LME
2. Run ZooMSS for full historical period (1960-2010) with:
   - Environmental forcing from GFDL-MOM6-COBALT2
   - Fishing mortality using calibrated parameters
3. Generate outputs on 1° and 0.25° global grids:
   - Fish biomass by size and group
   - Fishing mortality fields F(w,t)
   - Predicted catch by group and region
   - Size distributions over time

### Phase 6: Validation and Analysis
1. Compare predicted vs. observed catch across LMEs
2. Analyze spatial and temporal patterns of fishing impacts
3. Examine size structure changes due to fishing
4. Document parameter estimates and their uncertainties
5. Assess model performance and identify areas for improvement

---

## 10. Key Assumptions and Limitations

### 10.1 Assumptions
1. **Catchability is constant over time** (within each calibration period and LME)
   - In reality, catchability may change with technology, fleet behavior, and fish availability
2. **Knife-edge selectivity** is an adequate approximation
   - Real selectivity curves are typically logistic or dome-shaped
3. **Fishing effort data accurately represents actual fishing pressure**
   - Effort data may have biases, reporting errors, or gaps
4. **FishMIP functional groups map cleanly to ZooMSS size groups**
   - Some species may not fit neatly into size categories
5. **Regional calibration (LME-scale) captures key spatial heterogeneity**
   - Within-LME variation is not captured
6. **Environmental forcing from GFDL-MOM6-COBALT2 adequately represents historical conditions**
   - Model biases or errors in forcing data will propagate to ecosystem predictions

### 10.2 Limitations
1. **Parameter Identifiability:** 
   - With only catch data, q and w_min may be partially correlated
   - Additional data (e.g., size composition of catch) would improve estimates
2. **Aggregation:** 
   - Mapping 6 FishMIP groups to 3 ZooMSS groups loses some functional diversity
3. **Size Threshold Simplicity:** 
   - Knife-edge selectivity ignores gradual transitions in vulnerability
4. **Temporal Resolution:** 
   - Annual or monthly data may miss seasonal or sub-annual dynamics
5. **Data Quality:** 
   - Catch and effort data have uncertainties, especially for historical periods
6. **Ecological Interactions:** 
   - Changes in predation, competition, and food web dynamics may complicate calibration

---

## 11. Expected Outputs

### 11.1 Calibrated Parameters
- 6 parameters (q and w_min) for each LME or region:
  - Small fish: q_small, w_min_small
  - Medium fish: q_medium, w_min_medium  
  - Large fish: q_large, w_min_large

### 11.2 Model Outputs (1960-2010)
- **Fishing Mortality Fields:** F(w,t) for each grid cell and time step
- **Fish Biomass:** Total and size-structured biomass by fish group
- **Predicted Catch:** Annual catch by fish group and region
- **Size Distributions:** Abundance N(w,t) over time showing fishing impacts
- **Ecosystem Indicators:** Metrics like mean fish size, biomass ratios, etc.

### 11.3 Validation Metrics
- **Fit Quality:** R², RMSE, bias between predicted and observed catch
- **Parameter Uncertainty:** Confidence intervals or posterior distributions (if Bayesian)
- **Residual Analysis:** Temporal and spatial patterns in prediction errors
- **Cross-Validation:** Performance on held-out data or time periods

---

## 12. Technical Considerations

### 12.1 Computational Requirements
- **Optimization:** Each optimization run requires 10s-100s of ZooMSS simulations
- **Parallelization:** Regional calibrations are independent and can be parallelized
- **Storage:** Time-resolved spatial outputs can be large (GB-TB for global grids)
- **Runtime:** Depends on ZooMSS computational cost; budget hours-days per LME

### 12.2 Numerical Stability
- **Parameter Bounds:** Essential to prevent unrealistic values
- **Initial Conditions:** Choose starting values near expected optima to speed convergence
- **Regularization:** Consider adding penalties for extreme parameter values if needed
- **Scaling:** Normalize objective function or parameters if magnitudes differ greatly

### 12.3 Uncertainty Quantification
- **Parameter Uncertainty:** Use Hessian matrix from L-BFGS-B or MCMC sampling
- **Propagation:** Run ensemble simulations with parameter uncertainty
- **Sensitivity Analysis:** Test robustness to parameter variations
- **Scenario Analysis:** Examine impacts of assumption violations

---

## 13. Extensions and Future Work

### 13.1 Potential Enhancements
1. **Logistic Selectivity:** Replace knife-edge with smoother selectivity curves
2. **Time-Varying Catchability:** Allow q to change over time (technology improvements)
3. **Multi-Objective Optimization:** Include size composition data alongside catch
4. **Fleet Disaggregation:** Model different gear types or fleets explicitly
5. **Dynamic Effort Allocation:** Predict effort distribution based on fish abundance
6. **Climate-Fishery Interactions:** Examine how climate change affects catchability

### 13.2 Model Intercomparison
- Compare ZooMSS fishing implementation with other FishMIP models
- Participate in FishMIP ensemble analyses
- Benchmark against models with more complex fishing representations

---

## 14. References and Resources

### Key Papers
1. **Rousseau et al. 2024.** "FishMIP ISIMIP3a Fishing Effort Data." *Nature Scientific Data*.
   - Source of historical fishing effort data
   
2. **FishMIP Protocol Documentation** (ISIMIP website)
   - Detailed specifications for model participation

3. **ZooMSS Model Description Paper** (if available)
   - Core model structure and equations

### Data Sources
- **FishMIP:** https://www.isimip.org/gettingstarted/marine-ecosystems-fisheries/
- **GFDL-MOM6-COBALT2:** CMIP6 or ISIMIP data repositories
- **FAO Fisheries Statistics:** http://www.fao.org/fishery/statistics/en

---

## 15. Contact and Collaboration

**Primary Contact:** Kieran Murphy  
**Affiliation:** University of Queensland  
**Role:** ZooMSS Development and FishMIP Coordination

For questions about:
- ZooMSS model specifics
- Parameter interpretation
- Calibration strategies
- FishMIP integration

---

## Appendix A: Notation Summary

| Symbol | Description | Units |
|--------|-------------|-------|
| F(w,t) | Fishing mortality rate at size w and time t | per time |
| Effort(t) | Observed fishing effort at time t | kW-days |
| q | Catchability coefficient | dimensionless |
| w | Body mass (individual weight) | grams (g) |
| w_min | Minimum vulnerable body mass | grams (g) |
| N(w,t) | Abundance at size w and time t | numbers |
| Catch(t) | Total catch at time t | tonnes/year |
| SST(x,y,t) | Sea surface temperature | °C |
| Chl(x,y,t) | Chlorophyll-a concentration | mg/m³ |
| SSE | Sum of squared errors | (tonnes)² |

---

## Appendix B: Quick Reference Checklist

- [ ] FishMIP fishing effort data downloaded and processed
- [ ] Environmental forcing (SST, Chl) from GFDL-MOM6-COBALT2 obtained
- [ ] Observed catch data compiled by functional group and region
- [ ] LME boundaries defined for regional calibration
- [ ] Functional group mapping (FishMIP → ZooMSS) established
- [ ] Fishing mortality equation F(w,t) implemented in ZooMSS
- [ ] Knife-edge selectivity function coded
- [ ] Catch calculation integrated into model output
- [ ] Parameter bounds set for all 6 parameters
- [ ] Optimization algorithm (e.g., L-BFGS-B) configured
- [ ] Calibration loop tested on single LME
- [ ] Parallel computing infrastructure set up for multiple LMEs
- [ ] Validation metrics defined and coded
- [ ] Output storage and post-processing scripts ready
- [ ] Documentation of methods prepared for publication

---

**Document Version:** 1.0  
**Last Updated:** December 2024  
**Status:** Ready for Implementation
