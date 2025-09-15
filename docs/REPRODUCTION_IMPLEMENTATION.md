# ZooMSS Fish Reproduction Implementation

## Overview
This implementation adds dynamic fish reproduction to the ZooMSS model based on spawning stock biomass (SSB) and Beverton-Holt stock-recruitment relationships.

## Files Generated
- `R/zoomss_reproduction.R` - Core reproduction functions
- `R/update_groups_for_reproduction.R` - Parameter update scripts
- `docs/zoomss_run_modifications.txt` - Modifications for main loop
- `docs/zoomss_setup_modifications.txt` - Modifications for setup
- `tests/test_reproduction.R` - Test script

## Implementation Steps

1. **Update GroupInputs parameters**:
   ```r
   source("R/update_groups_for_reproduction.R")
   Groups <- update_groups_csv("data-raw/GroupInputs.csv")
   ```

2. **Modify existing files**:
   - Follow instructions in `docs/zoomss_run_modifications.txt`
   - Follow instructions in `docs/zoomss_setup_modifications.txt`

3. **Test the implementation**:
   ```r
   source("tests/test_reproduction.R")
   ```

## Git Workflow

```bash
# Create new branch for reproduction feature
git checkout -b dev-reproduction

# Add new files
git add R/zoomss_reproduction.R
git add R/update_groups_for_reproduction.R
git commit -m "feat: Add reproduction function modules"

# After modifying existing files
git add R/zoomss_run.R R/zoomss_setup.R
git commit -m "feat: Integrate reproduction into simulation"

# Push to dev branch
git push origin dev-reproduction

# Create pull request to dev branch
```

## Key Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| ReproInvest | Investment at maturity | 0.06-0.10 |
| ReproExp | Size scaling | -0.25 |
| h (steepness) | Stock-recruitment | 0.75 |

## Testing
Run the test script to verify:
- Energy balance is maintained
- SSB responds to population changes
- Recruitment follows Beverton-Holt curve
- Size spectra remain stable

