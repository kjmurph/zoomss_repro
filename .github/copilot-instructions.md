# Copilot Agent Instructions for ZooMSS

You are assisting with development of the `zoomss` R package. Follow these principles:

---

## Reference Documentation

### ZooMSS Fishing Implementation Guide
The file `docs/ZooMSS_Fishing_Implementation_Guide.md` provides the technical specification for implementing historical fishing mortality in ZooMSS. This guide outlines the ISIMIP3a protocol implementation (1960-2010) using a catchability-based fishing mortality framework.

**Before implementing any changes related to fishing mortality:**

1. **Read the guide thoroughly** — Understand the full context before writing code.
2. **Critically assess the specification** — Identify any:
   - Ambiguities or underspecified details
   - Potential numerical or computational issues
   - Assumptions that may need validation
   - Edge cases not addressed in the guide
   - Inconsistencies between sections
3. **Raise concerns first** — If you identify issues with the guide, discuss them before implementing. Do not silently "fix" perceived problems.
4. **Propose before implementing** — For significant additions or interpretations beyond what the guide specifies, outline your approach and rationale for approval.

The guide is a working document. Critical feedback that improves the implementation is encouraged.

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
