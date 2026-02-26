# =============================================================================
# probe_zoomss_structure.R
#
# Run this FIRST to understand the exact dimensions of your ZooMSS output.
# Paste the printed output back so the analysis script can be finalised.
# =============================================================================

cache_dir <- here::here("vignettes", "cache", "repro_eff_cache")

results_repro <- readRDS(file.path(cache_dir, "results_uniform_repro.rds"))

res <- results_repro[[1]][[2]]   # repro_eff=1, effort=0.025 (non-zero effort run)

cat("=== Top-level names ===\n")
cat(paste(names(res), collapse = "\n"), "\n\n")

probe <- function(name, obj) {
  cat(sprintf("%-25s  class: %-12s", name, paste(class(obj), collapse="/")))
  if (is.array(obj) || is.matrix(obj)) {
    cat("  dim:", paste(dim(obj), collapse=" x "))
  } else if (is.vector(obj) || is.numeric(obj)) {
    cat("  length:", length(obj))
    if (length(obj) <= 6) cat("  values:", paste(round(obj, 4), collapse=", "))
  } else if (is.data.frame(obj)) {
    cat("  nrow:", nrow(obj), " ncol:", ncol(obj))
  } else if (is.list(obj)) {
    cat("  list length:", length(obj))
  }
  cat("\n")
}

for (nm in names(res)) probe(nm, res[[nm]])

cat("\n=== Key array details ===\n")

# time
cat("time: length =", length(res$time),
    " range =", min(res$time), "to", max(res$time), "\n")

# abundance / biomass — show dim names if any
for (nm in c("abundance", "biomass", "biomassC", "catch",
             "SSB", "recruitment", "Fmort_ts", "repro_rate",
             "total_repro_output", "growth", "mortality")) {
  if (!is.null(res[[nm]])) {
    obj <- res[[nm]]
    cat(sprintf("\n%s:\n", nm))
    cat("  class:", paste(class(obj), collapse="/"), "\n")
    if (!is.null(dim(obj))) {
      cat("  dim:", paste(dim(obj), collapse=" x "), "\n")
      dn <- dimnames(obj)
      if (!is.null(dn)) {
        for (d in seq_along(dn)) {
          if (!is.null(dn[[d]]))
            cat(sprintf("  dimnames[[%d]]: %s\n", d,
                        paste(head(dn[[d]], 10), collapse=", ")))
        }
      }
    } else {
      cat("  length:", length(obj), "\n")
      if (length(obj) <= 10) cat("  values:", paste(round(obj, 6), collapse=", "), "\n")
    }
    # Show a small slice of values
    if (is.numeric(obj)) {
      cat("  range of values:", round(min(obj, na.rm=TRUE), 6),
          "to", round(max(obj, na.rm=TRUE), 6), "\n")
    }
  }
}

cat("\n=== param slot (Groups info) ===\n")
if (!is.null(res$param)) {
  cat("param class:", class(res$param), "\n")
  if (is.data.frame(res$param)) {
    cat("param cols:", paste(names(res$param), collapse=", "), "\n")
    if ("Type" %in% names(res$param))
      cat("Types:", paste(res$param$Type, collapse=", "), "\n")
    if ("Species" %in% names(res$param))
      cat("Species:", paste(res$param$Species, collapse=", "), "\n")
  } else if (is.list(res$param)) {
    cat("param names:", paste(names(res$param), collapse=", "), "\n")
  }
}
