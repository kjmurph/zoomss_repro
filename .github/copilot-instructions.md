# Copilot Agent Instructions for ZooMSS

You are assisting with development of the `zoomss` R package. Follow these principles:

---

## Architecture Context

Development follows a **two-repository architecture**:

1. **`zoomss` (this package)** — Core model with generic fishing mortality functionality. Standalone and reusable, similar to how `mizer` operates.

2. **`zoomss-isimip3a` (separate repo)** — ISIMIP3a protocol implementation. Uses `zoomss` as a dependency. Based on the DBPM framework.

When working in this package, implement **generic, reusable functionality** — not protocol-specific code.

---

## DBPM Reference Repository

The DBPM ISIMIP3a implementation serves as the template for the `zoomss-isimip3a` workflow:

**Repository:** `https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a`  
**Branch:** `new_features`

**When working on `zoomss-isimip3a` tasks, fetch and study the corresponding DBPM script BEFORE implementing:**

| ZooMSS Task | DBPM Script to Fetch |
|-------------|---------------------|
| Climate input processing | `scripts/01_processing_dbpm_global_inputs.ipynb` |
| Regional input processing | `scripts/02_processing_dbpm_regional_inputs.py` |
| Fishing effort processing | `scripts/03_processing_effort_fishing_inputs.R` |
| **Calibration (critical)** | `scripts/04_calculating_dbpm_fishing_params.R` |
| Gridded model setup | `scripts/05_setup_gridded_DBPM.py` |
| Model execution | `scripts/06_running_gridded_DBPM.py` |
| Catch calculation | `scripts/07_calculating_catches_DBPM.py` |
| Output plotting | `scripts/08_plotting_gridded_DBPM_outputs.ipynb` |

**How to use DBPM scripts:**
1. Fetch the script from the repo
2. Understand its structure, inputs, outputs, and data flow
3. Identify what must change for ZooMSS (functional groups, model interface, parameters)
4. Identify what can be reused directly (data loading, output formatting, plotting)
5. Implement the ZooMSS equivalent, preserving the workflow structure

---

## Reference Documentation

### ZooMSS Fishing Implementation Guide
The file `ZooMSS_Fishing_Implementation_Guide.md` provides the technical specification for implementing historical fishing mortality. This guide outlines the scientific approach and equations.

### ZooMSS Fishing Implementation Assessment
The file `ZooMSS_Fishing_Implementation_Assessment_v2.md` contains:
- **Part 1:** Critical assessment (ambiguities, numerical concerns)
- **Part 2:** Decision log with pre-made design choices
- **Part 3:** Implementation roadmap with task specifications
  - **Phase A (Tasks A.1–A.7):** `zoomss` package fishing module
  - **Phase B (Tasks B.1–B.7):** `zoomss-isimip3a` workflow (separate repo)
- **Part 5:** Specific instructions for Copilot

**When implementing fishing functionality in this package:**

1. **Follow Phase A tasks (A.1–A.7)** — These are the package-level implementations.
2. **Implement generic interfaces** — Support multiple selectivity types, both direct F and effort-based modes.
3. **Maintain backward compatibility** — Existing code must work unchanged when fishing is not specified.
4. **Reference the Decision Log (Part 2)** — Decisions P1–P4 apply to package design.
5. **Do not implement protocol-specific code** — FAO regions, ISIMIP calibration, etc. belong in `zoomss-isimip3a`.

---

## Code Completeness

- Never use placeholder code, `# ...`, ellipsis, or `TODO` comments to skip implementation. Write complete, functional code.
- Plan your approach before coding and share it with me first.
- Preserve existing functionality when editing files unless explicitly asked to remove it.
- This is scientific modeling code where numerical accuracy matters. Do not simplify equations or use approximations without explicit approval.

## File Management

- Minimise creation of new files. Prefer adding functions to existing, logically appropriate files.
- If a new file is genuinely needed, justify why before creating it.

## Coding Conventions

- Follow the existing coding style in the `zoomss` package and standard R package conventions.
- Use roxygen2 documentation for all exported functions.
- Write documentation that is clear and informative but concise—avoid verbose or redundant comments.
- Use meaningful variable and function names; let the code be self-documenting where possible.

## R Package Standards

- Respect the package structure (`R/`, `man/`, `data/`, `tests/`, etc.).
- Ensure any new dependencies are justified and added to DESCRIPTION.
- Write code compatible with `devtools::check()` and CRAN standards where applicable.

## Verification

- After writing code, explain how to test it or provide example usage.
- Ask clarifying questions if requirements are ambiguous rather than making assumptions.
