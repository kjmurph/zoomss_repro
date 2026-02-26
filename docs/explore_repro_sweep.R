# =============================================================================
# explore_repro_sweep.R
#
# Diagnostic and exploratory analysis for the reproductive efficiency sweep.
# Handles a potentially truncated Part 2 RDS file and explores Part 1 fully.
#
# Usage:
#   source("explore_repro_sweep.R")
#   # or: Rscript explore_repro_sweep.R
#
# Assumptions:
#   - ZooMSS package is loadable via devtools::load_all() from pkg_path
#   - Each model result is a list/data.frame with columns including time,
#     and biomass/yield columns per fish group (Small, Medium, Large fish)
#   - cache_dir points to the folder with the .rds files
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# ── Adjust these paths as needed ─────────────────────────────────────────────
pkg_path  <- here::here()
cache_dir <- here::here("vignettes", "cache", "repro_eff_cache")
plot_dir  <- here::here("vignettes", "cache", "repro_eff_cache", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ── Recreate sweep config (mirrors run_repro_eff_sweep.R) ────────────────────
# Load from saved config if available, otherwise reconstruct inline
config_file <- file.path(cache_dir, "sweep_config.rds")
if (file.exists(config_file)) {
  cfg <- readRDS(config_file)
  repro_eff_values   <- cfg$repro_eff_values
  effort_levels      <- cfg$effort_levels
  repro_scenarios    <- cfg$repro_scenarios
  scenario_names     <- cfg$scenario_names
  fish_names         <- cfg$fish_names
  cat("Loaded sweep config from file.\n")
} else {
  # Reconstruct from script defaults
  repro_eff_values <- c(1, 0.5, 0.1)
  effort_levels    <- c(0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3,
                        0.4, 0.5, 0.75, 1.0, 1.5, 2.0)
  repro_scenarios  <- list(
    "Uniform 1.0"   = c(1,    1,    1),
    "Uniform 0.1"   = c(0.1,  0.1,  0.1),
    "Uniform 0.01"  = c(0.01, 0.01, 0.01),
    "Gradient 10x"  = c(1,    0.1,  0.01),
    "Gradient 100x" = c(1,    0.01, 0.0001),
    "Mild gradient" = c(1,    0.5,  0.1),
    "Steep gradient"= c(1,    0.05, 0.001)
  )
  scenario_names   <- names(repro_scenarios)
  # fish_names must come from the package; load if possible
  tryCatch({
    devtools::load_all(pkg_path, quiet = TRUE)
    Groups_base <- getGroups()
    fish_names  <- Groups_base$Species[Groups_base$Type == "Fish"]
  }, error = function(e) {
    fish_names <<- c("Small fish", "Medium fish", "Large fish")
    warning("Could not load package; using generic fish names.")
  })
  cat("Reconstructed sweep config from defaults.\n")
}

n_repro   <- length(repro_eff_values)
n_effort  <- length(effort_levels)
n_scen    <- length(repro_scenarios)

cat(sprintf("Part 1: %d repro_eff x %d effort = %d models\n",
            n_repro, n_effort, n_repro * n_effort))
cat(sprintf("Part 2: %d scenarios x %d effort = %d models\n",
            n_scen, n_effort, n_scen * n_effort))


# =============================================================================
# SECTION 1 — Diagnostic: test whether the Part 2 file is readable
# =============================================================================

cat("\n─── Part 2 RDS diagnostic ───\n")

p2_path <- file.path(cache_dir, "results_group_repro.rds")
p2_ok   <- FALSE
results_group_repro <- NULL

if (!file.exists(p2_path)) {
  cat("  ✗ File not found:", p2_path, "\n")
} else {
  fsize_gb <- file.info(p2_path)$size / 1e9
  cat(sprintf("  File size: %.2f GB\n", fsize_gb))

  # A fully-written Part 2 file should be ≥ Part 1 file in rough proportion
  p1_path   <- file.path(cache_dir, "results_uniform_repro.rds")
  p1_size   <- if (file.exists(p1_path)) file.info(p1_path)$size else NA
  if (!is.na(p1_size)) {
    expected_min <- p1_size * (n_scen * n_effort) / (n_repro * n_effort) * 0.7
    cat(sprintf("  Part 1 size: %.2f GB  |  Rough expected min for Part 2: %.2f GB\n",
                p1_size / 1e9, expected_min / 1e9))
    if (file.info(p2_path)$size < expected_min) {
      cat("  ⚠ Part 2 file is SMALLER than expected — likely truncated.\n")
    }
  }

  cat("  Attempting readRDS ... ")
  results_group_repro <- tryCatch({
    obj <- readRDS(p2_path)
    p2_ok <- TRUE
    cat("SUCCESS ✓\n")
    obj
  }, error = function(e) {
    cat("FAILED ✗\n")
    cat("  Error:", conditionMessage(e), "\n")
    NULL
  })
}

if (p2_ok) {
  # Count how many list slots are non-NULL
  total_slots  <- n_scen * n_effort
  filled_slots <- sum(sapply(seq_len(n_scen), function(i)
    sum(sapply(seq_len(n_effort), function(j)
      !is.null(results_group_repro[[i]][[j]])))))
  cat(sprintf("  Slots populated: %d / %d\n", filled_slots, total_slots))
  if (filled_slots < total_slots) {
    cat("  ⚠ Some slots are NULL — file may be partially written.\n")
  } else {
    cat("  All slots populated ✓\n")
  }
} else {
  cat("\n  Part 2 file is not usable. Options:\n")
  cat("  1. Re-run Part 2 only (see re_run_part2.R scaffold below).\n")
  cat("  2. Proceed with Part 1 analysis only.\n")
}


# =============================================================================
# SECTION 2 — Helper functions to extract summaries from a model result
# =============================================================================
# Adjust column names here to match your actual zoomss_model() output format.
# The helpers are written defensively so they fail informatively.

# Return the last `tail_years` of simulation as the quasi-steady state
get_tail <- function(res, tail_years = 50) {
  if (is.null(res)) return(NULL)
  # Support both data.frame and list-with-$output formats
  if (is.list(res) && !is.data.frame(res) && !is.null(res$output))
    res <- res$output
  if (!is.data.frame(res)) {
    warning("Unexpected model result format"); return(NULL)
  }
  if (!"time" %in% names(res)) {
    warning("No 'time' column found"); return(NULL)
  }
  t_max  <- max(res$time)
  res[res$time >= (t_max - tail_years), ]
}

# Detect fish biomass columns (adjust pattern if your naming differs)
fish_biomass_cols <- function(df) {
  grep("biomass.*fish|fish.*biomass|BioFish|bio_fish|^B_",
       names(df), value = TRUE, ignore.case = TRUE)
}

# Detect fish yield columns
fish_yield_cols <- function(df) {
  grep("yield.*fish|fish.*yield|catch|Yield",
       names(df), value = TRUE, ignore.case = TRUE)
}

# Summarise a single model run: mean biomass and yield at quasi-steady state
summarise_run <- function(res, tail_years = 50) {
  tail <- get_tail(res, tail_years)
  if (is.null(tail)) return(NULL)

  bio_cols   <- fish_biomass_cols(tail)
  yield_cols <- fish_yield_cols(tail)

  bio_means   <- if (length(bio_cols)   > 0) colMeans(tail[bio_cols],   na.rm = TRUE) else NULL
  yield_means <- if (length(yield_cols) > 0) colMeans(tail[yield_cols], na.rm = TRUE) else NULL
  total_yield <- if (!is.null(yield_means)) sum(yield_means) else NA

  list(
    bio_means   = bio_means,
    yield_means = yield_means,
    total_yield = total_yield,
    n_rows      = nrow(tail)
  )
}

# Coexistence metric: are all fish groups above a minimum biomass threshold?
coexistence_check <- function(res, threshold = 1e-6, tail_years = 50) {
  tail <- get_tail(res, tail_years)
  if (is.null(tail)) return(NA)
  bio_cols <- fish_biomass_cols(tail)
  if (length(bio_cols) == 0) return(NA)
  means <- colMeans(tail[bio_cols], na.rm = TRUE)
  all(means > threshold)
}


# =============================================================================
# SECTION 3 — Part 1 exploration: yield curves
# =============================================================================

cat("\n─── Loading Part 1 results ───\n")
results_repro <- readRDS(file.path(cache_dir, "results_uniform_repro.rds"))
cat("Loaded successfully.\n")

# ── 3a. Build summary data.frame ─────────────────────────────────────────────
cat("Extracting Part 1 summaries ...\n")

# First, peek at the structure of one result to identify columns
sample_res <- results_repro[[1]][[2]]  # [re_idx=1, eff_idx=2] (non-zero effort)
cat("\n--- Structure of one model result (head) ---\n")
if (is.data.frame(sample_res)) {
  cat("Data frame with", nrow(sample_res), "rows and", ncol(sample_res), "cols\n")
  cat("Columns:", paste(names(sample_res), collapse = ", "), "\n")
  print(head(sample_res, 3))
} else if (is.list(sample_res)) {
  cat("List with elements:", paste(names(sample_res), collapse = ", "), "\n")
  if (!is.null(sample_res$output)) {
    cat("$output columns:", paste(names(sample_res$output), collapse = ", "), "\n")
    print(head(sample_res$output, 3))
  }
}
cat("--------------------------------------------\n\n")

# Build summary rows — YOU MAY NEED TO ADJUST COLUMN NAMES BELOW
# after seeing the structure printout above.
p1_rows <- lapply(seq_len(n_repro), function(i) {
  lapply(seq_len(n_effort), function(j) {
    res <- results_repro[[i]][[j]]
    s   <- summarise_run(res)
    coex <- coexistence_check(res)

    data.frame(
      repro_eff   = repro_eff_values[i],
      effort      = effort_levels[j],
      total_yield = if (!is.null(s)) s$total_yield else NA,
      coexistence = coex,
      stringsAsFactors = FALSE
    )
  })
})

p1_df <- bind_rows(unlist(p1_rows, recursive = FALSE))
p1_df$repro_eff_label <- factor(paste0("repro_eff = ", p1_df$repro_eff))

cat("Part 1 summary (first 10 rows):\n")
print(head(p1_df, 10))

# Also try to extract per-group biomass and yield
p1_group_rows <- lapply(seq_len(n_repro), function(i) {
  lapply(seq_len(n_effort), function(j) {
    res  <- results_repro[[i]][[j]]
    tail <- get_tail(res)
    if (is.null(tail)) return(NULL)

    bio_cols   <- fish_biomass_cols(tail)
    yield_cols <- fish_yield_cols(tail)

    out <- data.frame(repro_eff = repro_eff_values[i],
                      effort    = effort_levels[j])
    if (length(bio_cols) > 0) {
      bm <- as.data.frame(t(colMeans(tail[bio_cols], na.rm = TRUE)))
      names(bm) <- paste0("bio_", seq_along(bio_cols))
      out <- cbind(out, bm)
    }
    if (length(yield_cols) > 0) {
      ym <- as.data.frame(t(colMeans(tail[yield_cols], na.rm = TRUE)))
      names(ym) <- paste0("yield_", seq_along(yield_cols))
      out <- cbind(out, ym)
    }
    out
  })
})
p1_group_df <- bind_rows(unlist(p1_group_rows, recursive = FALSE))

# ── 3b. Yield curves ──────────────────────────────────────────────────────────
if (!all(is.na(p1_df$total_yield))) {
  p1_yield_plot <- ggplot(p1_df, aes(x = effort, y = total_yield,
                                     colour = repro_eff_label,
                                     group  = repro_eff_label)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_viridis_d(option = "plasma", end = 0.85) +
    labs(title   = "Part 1: Total yield vs fishing effort",
         subtitle = "Uniform repro_eff across all fish groups",
         x       = "Fishing effort",
         y       = "Total yield (quasi-steady state mean)",
         colour  = NULL) +
    theme_bw(base_size = 13)

  ggsave(file.path(plot_dir, "p1_yield_curves.png"),
         p1_yield_plot, width = 8, height = 5, dpi = 150)
  cat("Saved: p1_yield_curves.png\n")
}

# ── 3c. Per-group biomass at steady state (no fishing) ───────────────────────
# Focus on effort = 0 to see how repro_eff alone changes biomass distribution
if (ncol(p1_group_df) > 2) {
  bio_cols_in_df <- grep("^bio_", names(p1_group_df), value = TRUE)
  if (length(bio_cols_in_df) > 0) {
    p1_bio_long <- p1_group_df %>%
      filter(effort == 0) %>%
      pivot_longer(all_of(bio_cols_in_df),
                   names_to  = "group",
                   values_to = "biomass") %>%
      mutate(group = gsub("bio_", "Fish group ", group))

    p1_bio_plot <- ggplot(p1_bio_long,
                          aes(x = factor(repro_eff), y = biomass, fill = group)) +
      geom_col(position = "dodge") +
      scale_fill_viridis_d(option = "turbo", end = 0.85) +
      labs(title   = "Part 1: Fish group biomass at zero effort",
           subtitle = "Quasi-steady state mean (last 50 yr)",
           x       = "repro_eff (uniform)",
           y       = "Biomass",
           fill    = "Group") +
      theme_bw(base_size = 13)

    ggsave(file.path(plot_dir, "p1_biomass_zero_effort.png"),
           p1_bio_plot, width = 7, height = 5, dpi = 150)
    cat("Saved: p1_biomass_zero_effort.png\n")
  }
}

# ── 3d. Coexistence heatmap ───────────────────────────────────────────────────
if (!all(is.na(p1_df$coexistence))) {
  p1_coex_plot <- ggplot(p1_df, aes(x = effort, y = factor(repro_eff),
                                    fill = coexistence)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = c("TRUE" = "#2ecc71", "FALSE" = "#e74c3c",
                                 "NA"   = "grey80"),
                      na.value = "grey80",
                      labels  = c("TRUE" = "Co-exist", "FALSE" = "Collapse")) +
    labs(title  = "Part 1: Coexistence of all fish groups",
         x      = "Fishing effort",
         y      = "repro_eff (uniform)",
         fill   = NULL) +
    theme_bw(base_size = 13)

  ggsave(file.path(plot_dir, "p1_coexistence_heatmap.png"),
         p1_coex_plot, width = 8, height = 4, dpi = 150)
  cat("Saved: p1_coexistence_heatmap.png\n")
}

# ── 3e. Time series for representative runs ───────────────────────────────────
# Plot time series for re=1.0 at low / medium / high effort
ts_efforts <- c(0, 0.1, 0.5, 1.0)
ts_effort_indices <- sapply(ts_efforts, function(e)
  which.min(abs(effort_levels - e)))

repro_idx <- 1  # repro_eff = 1.0

ts_list <- lapply(ts_effort_indices, function(j) {
  res <- results_repro[[repro_idx]][[j]]
  # Normalise to data.frame
  if (is.list(res) && !is.data.frame(res) && !is.null(res$output)) res <- res$output
  if (!is.data.frame(res)) return(NULL)
  bio_cols <- fish_biomass_cols(res)
  if (!"time" %in% names(res) || length(bio_cols) == 0) return(NULL)
  res %>%
    select(time, all_of(bio_cols)) %>%
    pivot_longer(-time, names_to = "group", values_to = "biomass") %>%
    mutate(effort_label = paste0("effort = ", effort_levels[j]))
})
ts_df <- bind_rows(Filter(Negate(is.null), ts_list))

if (nrow(ts_df) > 0) {
  ts_plot <- ggplot(ts_df, aes(x = time, y = biomass,
                               colour = group, group = group)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    facet_wrap(~effort_label, scales = "free_y", ncol = 2) +
    scale_colour_viridis_d(option = "turbo", end = 0.85) +
    labs(title   = "Part 1: Fish biomass time series (repro_eff = 1.0)",
         x       = "Time (years)",
         y       = "Biomass",
         colour  = "Group") +
    theme_bw(base_size = 12)

  ggsave(file.path(plot_dir, "p1_timeseries_repro1.png"),
         ts_plot, width = 10, height = 6, dpi = 150)
  cat("Saved: p1_timeseries_repro1.png\n")
}


# =============================================================================
# SECTION 4 — Part 2 exploration (only if file loaded successfully)
# =============================================================================

if (!p2_ok || is.null(results_group_repro)) {
  cat("\n─── Part 2: Skipped (file not readable) ───\n")
} else {
  cat("\n─── Part 2 exploration ───\n")

  # ── 4a. Summary data.frame ─────────────────────────────────────────────────
  p2_rows <- lapply(seq_len(n_scen), function(i) {
    lapply(seq_len(n_effort), function(j) {
      res  <- results_group_repro[[i]][[j]]
      if (is.null(res)) {
        return(data.frame(scenario = scenario_names[i],
                          effort   = effort_levels[j],
                          total_yield = NA, coexistence = NA,
                          stringsAsFactors = FALSE))
      }
      s    <- summarise_run(res)
      coex <- coexistence_check(res)
      data.frame(scenario    = scenario_names[i],
                 effort      = effort_levels[j],
                 total_yield = if (!is.null(s)) s$total_yield else NA,
                 coexistence = coex,
                 stringsAsFactors = FALSE)
    })
  })
  p2_df <- bind_rows(unlist(p2_rows, recursive = FALSE))

  cat("Part 2 summary (first 10 rows):\n")
  print(head(p2_df, 10))

  # Report completeness
  n_null <- sum(is.na(p2_df$coexistence))
  cat(sprintf("  NULL slots: %d / %d\n", n_null, nrow(p2_df)))

  # ── 4b. Yield curves by scenario ──────────────────────────────────────────
  if (!all(is.na(p2_df$total_yield))) {
    # Separate baselines vs gradients for clarity
    baseline_scens <- grep("Uniform", scenario_names, value = TRUE)
    gradient_scens <- grep("gradient|Gradient", scenario_names, value = TRUE)

    p2_yield_plot <- ggplot(p2_df %>% filter(!is.na(total_yield)),
                            aes(x = effort, y = total_yield,
                                colour = scenario, group = scenario)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2) +
      scale_colour_viridis_d(option = "plasma", end = 0.9) +
      labs(title   = "Part 2: Total yield vs effort by repro_eff scenario",
           x       = "Fishing effort", y = "Total yield",
           colour  = "Scenario") +
      theme_bw(base_size = 12) +
      theme(legend.position = "right")

    ggsave(file.path(plot_dir, "p2_yield_curves.png"),
           p2_yield_plot, width = 10, height = 5, dpi = 150)
    cat("Saved: p2_yield_curves.png\n")
  }

  # ── 4c. Coexistence heatmap across scenarios and effort ───────────────────
  if (!all(is.na(p2_df$coexistence))) {
    p2_coex_plot <- ggplot(p2_df %>% filter(!is.na(coexistence)),
                           aes(x = effort, y = scenario, fill = coexistence)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_manual(values = c("TRUE" = "#2ecc71", "FALSE" = "#e74c3c"),
                        labels = c("TRUE" = "Co-exist", "FALSE" = "Collapse")) +
      labs(title  = "Part 2: Coexistence heatmap by scenario",
           x      = "Fishing effort", y = NULL, fill = NULL) +
      theme_bw(base_size = 12)

    ggsave(file.path(plot_dir, "p2_coexistence_heatmap.png"),
           p2_coex_plot, width = 9, height = 5, dpi = 150)
    cat("Saved: p2_coexistence_heatmap.png\n")
  }

  # ── 4d. Group biomass at zero effort across scenarios ────────────────────
  p2_bio_rows <- lapply(seq_len(n_scen), function(i) {
    res  <- results_group_repro[[i]][[1]]  # effort_idx=1 → effort=0
    if (is.null(res)) return(NULL)
    tail <- get_tail(res)
    if (is.null(tail)) return(NULL)
    bio_cols <- fish_biomass_cols(tail)
    if (length(bio_cols) == 0) return(NULL)
    bm <- colMeans(tail[bio_cols], na.rm = TRUE)
    data.frame(scenario = scenario_names[i],
               group    = names(bm),
               biomass  = as.numeric(bm))
  })
  p2_bio_df <- bind_rows(Filter(Negate(is.null), p2_bio_rows))

  if (nrow(p2_bio_df) > 0) {
    # Order scenarios sensibly
    p2_bio_df$scenario <- factor(p2_bio_df$scenario,
                                 levels = rev(scenario_names))
    p2_bio_plot <- ggplot(p2_bio_df,
                          aes(x = scenario, y = biomass, fill = group)) +
      geom_col(position = "dodge") +
      coord_flip() +
      scale_fill_viridis_d(option = "turbo", end = 0.85) +
      labs(title   = "Part 2: Fish biomass at zero effort by scenario",
           x       = NULL, y = "Biomass (quasi-steady state)", fill = "Group") +
      theme_bw(base_size = 12)

    ggsave(file.path(plot_dir, "p2_biomass_zero_effort.png"),
           p2_bio_plot, width = 9, height = 5, dpi = 150)
    cat("Saved: p2_biomass_zero_effort.png\n")
  }
}


# =============================================================================
# SECTION 5 — Cross-part comparison: Uniform 1.0 should match Part 1 re=1
# =============================================================================

if (p2_ok && !is.null(results_group_repro)) {
  cat("\n─── Cross-validation: Part 1 re=1.0 vs Part 2 'Uniform 1.0' ───\n")
  idx_u1 <- which(scenario_names == "Uniform 1.0")
  if (length(idx_u1) == 1) {
    p1_yields <- sapply(seq_len(n_effort), function(j) {
      s <- summarise_run(results_repro[[1]][[j]])
      if (!is.null(s)) s$total_yield else NA
    })
    p2_yields <- sapply(seq_len(n_effort), function(j) {
      s <- summarise_run(results_group_repro[[idx_u1]][[j]])
      if (!is.null(s)) s$total_yield else NA
    })
    cv_df <- data.frame(effort = effort_levels,
                        Part1  = p1_yields,
                        Part2  = p2_yields) %>%
      pivot_longer(c(Part1, Part2), names_to = "source", values_to = "yield")

    cv_plot <- ggplot(cv_df %>% filter(!is.na(yield)),
                      aes(x = effort, y = yield,
                          colour = source, linetype = source)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      labs(title = "Cross-check: Part 1 repro_eff=1 vs Part 2 'Uniform 1.0'",
           x = "Effort", y = "Total yield") +
      theme_bw(base_size = 13)

    ggsave(file.path(plot_dir, "crosscheck_p1_vs_p2_uniform.png"),
           cv_plot, width = 7, height = 4, dpi = 150)
    cat("Saved: crosscheck_p1_vs_p2_uniform.png\n")

    # Numeric comparison
    diff_df <- data.frame(effort = effort_levels,
                          p1 = p1_yields, p2 = p2_yields,
                          rel_diff_pct = 100 * abs(p1_yields - p2_yields) / p1_yields)
    cat("Relative differences (should be ~0 if deterministic):\n")
    print(diff_df)
  }
}


# =============================================================================
# SECTION 6 — Console summary & key findings
# =============================================================================

cat("\n==================================================\n")
cat("SUMMARY\n")
cat("==================================================\n")
cat(sprintf("Part 1 file:  %s (%s)\n",
            file.path(cache_dir, "results_uniform_repro.rds"),
            if (file.exists(file.path(cache_dir, "results_uniform_repro.rds"))) "OK" else "MISSING"))
cat(sprintf("Part 2 file:  %s (%s)\n",
            p2_path,
            if (p2_ok) "READABLE" else "UNREADABLE/TRUNCATED"))

if (!all(is.na(p1_df$total_yield))) {
  cat("\nPart 1 key results:\n")
  # Max yield and effort-at-max-yield by repro_eff
  p1_summary <- p1_df %>%
    filter(!is.na(total_yield)) %>%
    group_by(repro_eff) %>%
    summarise(max_yield  = max(total_yield),
              opt_effort = effort[which.max(total_yield)],
              .groups    = "drop")
  print(p1_summary)
}

if (p2_ok && !all(is.na(p2_df$total_yield))) {
  cat("\nPart 2 key results (max yield per scenario):\n")
  p2_summary <- p2_df %>%
    filter(!is.na(total_yield)) %>%
    group_by(scenario) %>%
    summarise(max_yield  = max(total_yield),
              opt_effort = effort[which.max(total_yield)],
              .groups    = "drop") %>%
    arrange(desc(max_yield))
  print(p2_summary)
}

cat("\nPlots saved to:", normalizePath(plot_dir), "\n")

# =============================================================================
# SECTION 7 — Scaffold: re-run Part 2 only (if needed)
# =============================================================================

cat("\n─── Scaffold for re-running Part 2 (if needed) ───\n")
cat("If Part 2 was unreadable, save the following as re_run_part2.R:\n\n")
cat(
'# re_run_part2.R  — paste into your project and run
devtools::load_all()
library(parallel)
source("run_repro_eff_sweep.R", local = TRUE)  # loads helpers & config

# Only re-run Part 2 ----------------------------------------------------------
group_results_flat <- parLapply(cl, seq_len(nrow(group_jobs)), function(j) {
  sc_idx  <- group_jobs$sc_idx[j]
  eff_idx <- group_jobs$eff_idx[j]
  repro_vec  <- repro_scenarios[[sc_idx]]
  effort_vec <- rep(effort_levels[eff_idx], num_fish)
  q_vec      <- rep(q_fixed, num_fish)
  run_single(repro_vec, effort_vec, Groups_base, fish_rows, q_vec,
             sim_time, sst_const, chl_const, isave, effort_col_names)
})
results_group_repro <- vector("list", length(repro_scenarios))
for (i in seq_along(results_group_repro))
  results_group_repro[[i]] <- vector("list", length(effort_levels))
for (j in seq_len(nrow(group_jobs))) {
  results_group_repro[[group_jobs$sc_idx[j]]][[group_jobs$eff_idx[j]]] <-
    group_results_flat[[j]]
}
saveRDS(results_group_repro,
        file.path(cache_dir, "results_group_repro.rds"))
')
