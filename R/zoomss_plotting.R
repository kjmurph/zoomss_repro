#' Plot Predator-Prey Mass Ratio (PPMR)
#'
#' @title Visualize predator-prey mass ratio patterns in ZooMSS results
#' @description Creates a plot showing the distribution of predator-prey mass ratios (PPMR)
#'   across functional groups, providing insights into the trophic structure of the ecosystem.
#' @details This function calculates and visualizes PPMR patterns by:
#'   - Computing theoretical PPMR values for each functional group and size class
#'   - Weighting by biomass to show realized community patterns
#'   - Creating a density plot of PPMR distribution across the community
#'   - Overlaying species-specific PPMR values as points
#'
#'   PPMR is a key ecological metric that describes the size relationship between
#'   predators and their prey, providing insight into food web structure and
#'   energy transfer efficiency in marine ecosystems.
#'
#' @param mdl ZooMSS results object containing model outputs and parameters
#' @param idx The time index to plot
#'
#' @return ggplot object showing PPMR distribution with species-specific overlays
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' ppmr_plot <- plotPPMR(results)
#' print(ppmr_plot)
#' }
#'
plotPPMR <- function(mdl, idx){

  out <- extractPPMR(mdl)

  out <- out[[idx]] # Subset to the required timestep

  ggplot2::ggplot() +
    ggplot2::geom_line(data = out[[2]], mapping = ggplot2::aes(x = .data$Betas, y = .data$y, colour = .data$Species), linewidth = 1) +
    ggplot2::geom_line(data = out[[1]], mapping = ggplot2::aes(x = .data$x, y = .data$y), linewidth = 1.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    ggplot2::labs(x = expression('log' [10] * PPMR),
                  y = "Zoop. Biomass Proportion", subtitle = "PPMR") +
    ggplot2::geom_vline(data = out[[1]], mapping = ggplot2::aes(xintercept = .data$mn_beta), colour = 'black') +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = c("Flagellates" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Flagellates"],
                                            "Ciliates" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Ciliates"],
                                            "Larvaceans" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Larvaceans"],
                                            "Salps" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Salps"],
                                            "Jellyfish" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Jellyfish"],
                                            "CarnCopepods" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="CarnCopepods"],
                                            "Chaetognaths" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Chaetognaths"],
                                            "Euphausiids" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Euphausiids"],
                                            "OmniCopepods" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="OmniCopepods"]))
}

#' Plot Size Spectra for ZooMSS Results
#'
#' @title Visualize abundance size spectra across functional groups
#' @description Creates a log-log plot of abundance versus body size for all functional groups,
#'   showing the classic size spectrum pattern in marine ecosystems.
#' @details This function visualizes the abundance size spectrum by:
#'   - Converting abundance data to long format with body weights
#'   - Filtering out zero abundances to focus on active size classes
#'   - Creating log-log plots colored by functional group
#'   - Using species-specific colors defined in the Groups parameter table
#'
#'   Size spectra are fundamental patterns in marine ecology, typically showing
#'   declining abundance with increasing body size. This visualization helps
#'   assess model realism and identify dominant size classes within each
#'   functional group.
#'
#' @param mdl ZooMSS results object containing model outputs and parameters
#' @param by Character string specifying the metric to plot. Options: "abundance", "biomass", "mortality", "growth" (default: "abundance")
#' @param n_years The number of years (from the end) over which to average the size spectra
#'
#' @return ggplot object showing log abundance vs log body weight by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' size_plot <- plotSizeSpectra(results)
#' print(size_plot)
#' }
#'
plotSizeSpectra <- function(mdl, by = "abundance", n_years = 10) {

  species <- averageTimeSeries(mdl, var = by, n_years = n_years)

  rownames(species) <- mdl$param$Groups$Species
  species <- tibble::as_tibble(t(species))

  species <- species %>%
    tibble::add_column("Weight" = mdl$param$w) %>%
    tidyr::pivot_longer(-"Weight", names_to = "Species", values_to = by) %>%
    dplyr::filter(.data[[by]] > 0) %>%
    dplyr::mutate(Species = factor(.data$Species, levels = mdl$param$Groups$Species))

  ggplot2::ggplot(data = species,
                        mapping = ggplot2::aes(x = log10(.data$Weight),
                                               y = log10(.data[[by]]),
                                               colour = .data$Species)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = mdl$param$Groups$PlotColour) +
    ggplot2::theme_bw() +
    ggplot2::labs(subtitle = paste(stringr::str_to_title(by), "Spectrum"))

}

#' Plot Time Series Data for ZooMSS Results
#'
#' @title Unified function to visualize time series changes for different metrics
#' @description Creates time series plots showing how abundance, biomass, mortality, or growth
#'   rates of functional groups change throughout the ZooMSS simulation period.
#' @details This function creates time series visualizations by:
#'   - **Abundance**: Summing abundances across size classes, log-transformed y-axis
#'   - **Biomass**: Calculating biomass (abundance × weight), with optional stacking and proportional scaling
#'   - **Mortality**: Averaging predation mortality rates across size classes
#'   - **Growth**: Averaging growth rates across size classes, log-transformed y-axis
#'
#'   All plots use species-specific colors and filter out zero values. Time series plots help identify:
#'   - Equilibration time for model runs
#'   - Seasonal or cyclical patterns in ecological metrics
#'   - Relative patterns between functional groups
#'   - Model stability and convergence behavior
#'
#' @param mdl ZooMSS results object containing model outputs with time series data
#' @param by Character string specifying the metric to plot. Options: "abundance", "biomass", "mortality", "growth" (default: "abundance")
#' @param species Character vector of species names to include. If NULL, all species included (default: NULL, applies to all metrics)
#' @param type Character vector of plot type. Use `line` for the default line plot, `stack` or `fill` (as per geom_area) for stacked or proportional plots. (default: "line")
#' @param transform Character vector of the required y-axis transformation. Options from `scale_*_continuous` (Default: "identity).
#'
#' @return ggplot object showing the requested time series by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#'
#' # Plot different metrics
#' abundance_plot <- plotTimeSeries(results, by = "abundance", transform = "log10")
#' biomass_plot <- plotTimeSeries(results, by = "biomass", transform = "log10")
#' mortality_plot <- plotTimeSeries(results, by = "mortality")
#' growth_plot <- plotTimeSeries(results, by = "growth")
#'
#' stacked_plot <- plotTimeSeries(results, by = "biomass", type = "stack")
#' prop_plot <- plotTimeSeries(results, by = "biomass", type = "fill")
#'
#' # Focus on specific species (works for all metrics)
#' copepod_plot <- plotTimeSeries(results, by = "biomass",
#'                               species = c("OmniCopepods", "CarnCopepods"))
#' abundance_copepods <- plotTimeSeries(results, by = "abundance",
#'                                     species = c("OmniCopepods", "CarnCopepods"))
#' mortality_copepods <- plotTimeSeries(results, by = "mortality",
#'                                     species = c("OmniCopepods", "CarnCopepods"))
#' growth_copepods <- plotTimeSeries(results, by = "growth",
#'                                  species = c("OmniCopepods", "CarnCopepods"))
#' }
#'
plotTimeSeries <- function(mdl, by = "abundance", type = "line",
                           transform = "identity", species = NULL) {

  # Validate inputs
  assertthat::assert_that(
    is.list(mdl),
    msg = "mdl must be a list."
  )

  assertthat::assert_that(
    "abundance" %in% names(mdl),
    msg = "Time series data not available. Model may not have been run correctly."
  )
  assertthat::assert_that(
    by %in% c("abundance", "biomass", "mortality", "growth"),
    msg = paste0("'by' must be one of: ", paste(c("abundance", "biomass", "mortality", "growth"), collapse = ", "))
  )

  assertthat::assert_that(
    type %in% c("line", "stack", "fill"),
    msg = paste0("'type' must be one of: ", paste(c("line", "stack", "fill"), collapse = ", "))
  )

  assertthat::assert_that(
    transform %in% c("identity", "log", "log10", "log1p", "log2", "logit"),
    msg = paste0("'transform' must be one of: ", paste(c(c("identity", "log", "log10", "log1p", "log2", "logit")), collapse = ", "))
  )

  assertthat::assert_that(
    is.null(species) || all(species %in% mdl$param$Groups$Species),
    msg = paste0("All elements of 'species' must be in: ", paste(mdl$param$Groups$Species, collapse = ", "))
  )

  # Calculate data based on requested variable
  data_matrix <- switch( by,
                         "abundance" = mdl$abundance %>% reduceSize(),
                         "biomass" = mdl$biomass %>% reduceSize(),
                         "mortality" = mdl$mortality %>% reduceSize(method = "mean"),
                         "growth" = mdl$growth %>% reduceSize(method = "mean")
  )

  # Set up data frame
  colnames(data_matrix) <- mdl$param$Groups$Species
  data_df <- tibble::as_tibble(data_matrix)
  data_df$Time <- mdl$time

  # Filter species if specified
  if (!is.null(species)) {
    data_df <- data_df[, c("Time", species)]
  }

  # Convert to long format
  data_long <- data_df %>%
    tidyr::pivot_longer(-"Time", names_to = "Species", values_to = by) %>%
    dplyr::filter(!!rlang::sym(by) > 0) %>%
    dplyr::mutate(Species = factor(.data$Species, levels = mdl$param$Groups$Species))

  # Get colors for plotting
  if (!is.null(species)) {
    species_indices <- match(intersect(species, mdl$param$Groups$Species), mdl$param$Groups$Species)
    plot_colors <- mdl$param$Groups$PlotColour[species_indices]
    names(plot_colors) <- mdl$param$Groups$Species[species_indices]
  } else {
    plot_colors <- mdl$param$Groups$PlotColour
    names(plot_colors) <- mdl$param$Groups$Species
  }


  # Do plotting
  if (type %in% c("fill", "stack")) {

    # Stacked and proportion area plot for any metric
    gg <- ggplot2::ggplot(data = data_long,
                          mapping = ggplot2::aes(x = .data$Time,
                                                 y = !!rlang::sym(by),
                                                 fill = .data$Species)) +
      ggplot2::geom_area(position = type, alpha = 0.7) +
      ggplot2::scale_fill_manual(values = plot_colors) +
      ggplot2::scale_y_continuous(transform = transform, expand = c(0, 0))

  } else if (type == "line") {

    gg <- ggplot2::ggplot(data = data_long,
                          mapping = ggplot2::aes(x = .data$Time,
                                                 y = !!rlang::sym(by),
                                                 colour = .data$Species)) +
      ggplot2::scale_y_continuous(transform = transform) +
      ggplot2::geom_line(linewidth = 1) +
      # ggplot2::geom_point(size = 1.2) +
      ggplot2::scale_color_manual(values = plot_colors)

  }

  gg <- gg +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(y = stringr::str_to_sentence(by)) +
    ggplot2::xlab("Time (years)")
  return(gg)
}






#' Plot Environmental Time Series
#'
#' @title Plot environmental forcing data
#' @description Creates plots of sea surface temperature and chlorophyll time series
#'   for visualizing environmental forcing data used in ZooMSS model runs.
#' @details This function creates two separate plots with different y-axes scales:
#'   - SST plot (red line) with temperature in deg C
#'   - Chlorophyll plot (green line) with concentration in mg/m^3
#'
#'   The plots can be combined using the patchwork package if available, otherwise
#'   separate plots are returned as a list. This helps users visualize the
#'   environmental forcing that drives ZooMSS model dynamics.
#'
#' @param env_data Environmental data frame with time, sst, chl columns
#'
#' @return ggplot object (if patchwork available) or list of two ggplot objects
#' @export
#'
#' @examples
#' # Create sample data and plot
#' env_data <- data.frame(
#'   time = 1:100,
#'   dt = 0.01,
#'   sst = 15 + 3*sin(2*pi*(1:100)/50),
#'   chl = 0.5 + 0.2*cos(2*pi*(1:100)/50)
#' )
#' plots <- plotEnvironment(env_data)
#'
plotEnvironment <- function(env_data) {

  # Convert to long format for plotting
  env_long <- tidyr::pivot_longer(env_data,
                                  cols = c("sst", "chl"),
                                  names_to = "variable",
                                  values_to = "value")

  # Create separate y-axes for SST and chlorophyll
  p1 <- ggplot2::ggplot(data = dplyr::filter(env_long, .data$variable == "sst"),
                        ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_line(color = "red", linewidth = 1) +
    ggplot2::labs(y = "SST (deg C)", title = "Environmental Forcing") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  p2 <- ggplot2::ggplot(data = dplyr::filter(env_long, .data$variable == "chl"),
                        ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_line(color = "green", linewidth = 1) +
    ggplot2::labs(x = "Time (years)", y = "Chlorophyll (mg/m^3)") +
    ggplot2::theme_bw()

  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(p1, p2, ncol = 1))
  } else {
    cat("Install patchwork package to combine plots\n")
    return(list(sst_plot = p1, chl_plot = p2))
  }
}


#' Plot Reproduction Time Series
#'
#' @title Visualize fish reproduction metrics through time
#' @description Creates time series plots showing spawning stock biomass (SSB),
#'   recruitment, and total reproductive output for fish functional groups.
#' @details This function visualizes key fish reproduction metrics including:
#'   - **SSB (Spawning Stock Biomass)**: Total biomass of mature individuals
#'   - **Recruitment**: Number of new recruits entering the smallest size class
#'   - **Reproductive Output**: Total reproductive investment (biomass allocated to reproduction)
#'
#'   These metrics are essential for understanding population dynamics driven by
#'   the explicit energy budget and reproduction system.
#'
#' @param mdl ZooMSS results object containing model outputs with reproduction data
#' @param by Character string specifying the metric to plot. Options:
#'   "SSB" (Spawning Stock Biomass), "recruitment", "repro_output" (default: "SSB")
#' @param transform Character vector of the required y-axis transformation.
#'   Options from `scale_*_continuous` (Default: "identity")
#'
#' @return ggplot object showing the requested reproduction time series by fish group
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#'
#' # Plot different reproduction metrics
#' ssb_plot <- plotReproduction(results, by = "SSB")
#' recruitment_plot <- plotReproduction(results, by = "recruitment")
#' repro_plot <- plotReproduction(results, by = "repro_output")
#' }
#'
plotReproduction <- function(mdl, by = "SSB", transform = "identity") {

  # Validate inputs
  assertthat::assert_that(
    is.list(mdl),
    msg = "mdl must be a list."
  )

  assertthat::assert_that(
    "SSB" %in% names(mdl),
    msg = "Reproduction data not available. Ensure model was run with fish groups that have repro_on = 1."
  )

  assertthat::assert_that(
    by %in% c("SSB", "recruitment", "repro_output"),
    msg = paste0("'by' must be one of: ", paste(c("SSB", "recruitment", "repro_output"), collapse = ", "))
  )

  assertthat::assert_that(
    transform %in% c("identity", "log", "log10", "log1p", "log2"),
    msg = paste0("'transform' must be one of: ", paste(c("identity", "log", "log10", "log1p", "log2"), collapse = ", "))
  )

  # Get fish group names
  fish_grps <- mdl$param$fish_grps
  fish_names <- mdl$param$Groups$Species[fish_grps]
  fish_colors <- mdl$param$Groups$PlotColour[fish_grps]
  names(fish_colors) <- fish_names

  # Select data based on metric
  data_matrix <- switch(by,
                        "SSB" = mdl$SSB,
                        "recruitment" = mdl$recruitment,
                        "repro_output" = mdl$total_repro_output)

  # Create data frame
  colnames(data_matrix) <- fish_names
  data_df <- tibble::as_tibble(data_matrix)
  data_df$Time <- mdl$time

  # Convert to long format
  data_long <- data_df %>%
    tidyr::pivot_longer(-"Time", names_to = "Species", values_to = by) %>%
    dplyr::filter(!is.na(!!rlang::sym(by))) %>%
    dplyr::mutate(Species = factor(.data$Species, levels = fish_names))

  # Create y-axis label
  y_label <- switch(by,
                    "SSB" = "Spawning Stock Biomass",
                    "recruitment" = "Recruitment (individuals)",
                    "repro_output" = "Total Reproductive Output")

  # Create plot
  gg <- ggplot2::ggplot(data = data_long,
                        mapping = ggplot2::aes(x = .data$Time,
                                               y = !!rlang::sym(by),
                                               colour = .data$Species)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = fish_colors) +
    ggplot2::scale_y_continuous(transform = transform) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = y_label, x = "Time (years)",
                  subtitle = paste("Fish", y_label))

  return(gg)
}


#' Plot Reproduction Size Spectra
#'
#' @title Visualize reproductive investment across body sizes
#' @description Creates plots showing how reproductive investment is distributed
#'   across body sizes for fish functional groups.
#' @details This function visualizes:
#'   - Reproductive rate as a function of body size
#'   - The effect of the maturity ogive on reproduction allocation
#'   - Size-dependent patterns in reproductive investment
#'
#'   The plot shows the specific reproductive rate (g reproduction / g body weight / time)
#'   at each body size, averaged over the specified number of years.
#'
#' @param mdl ZooMSS results object containing model outputs with reproduction data
#' @param n_years The number of years (from the end) over which to average (default: 10)
#'
#' @return ggplot object showing reproductive rate vs body size by fish group
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' repro_size_plot <- plotReproSizeSpectra(results, n_years = 10)
#' }
#'
plotReproSizeSpectra <- function(mdl, n_years = 10) {

  # Validate inputs
  assertthat::assert_that(
    is.list(mdl),
    msg = "mdl must be a list."
  )

  assertthat::assert_that(
    "repro_rate" %in% names(mdl),
    msg = "Reproduction rate data not available."
  )

  # Get fish groups
  fish_grps <- mdl$param$fish_grps
  fish_names <- mdl$param$Groups$Species[fish_grps]
  fish_colors <- mdl$param$Groups$PlotColour[fish_grps]
  names(fish_colors) <- fish_names

  # Average reproductive rate over final n_years
  repro_avg <- averageTimeSeries(mdl, var = "repro_rate", n_years = n_years)

  # Subset to fish groups only
  repro_fish <- repro_avg[fish_grps, , drop = FALSE]
  rownames(repro_fish) <- fish_names

  # Create data frame
  species_df <- tibble::as_tibble(t(repro_fish))
  species_df <- species_df %>%
    tibble::add_column("Weight" = mdl$param$w) %>%
    tidyr::pivot_longer(-"Weight", names_to = "Species", values_to = "repro_rate") %>%
    dplyr::filter(.data$repro_rate > 0) %>%
    dplyr::mutate(Species = factor(.data$Species, levels = fish_names))

  # Create plot
  gg <- ggplot2::ggplot(data = species_df,
                        mapping = ggplot2::aes(x = log10(.data$Weight),
                                               y = .data$repro_rate,
                                               colour = .data$Species)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(values = fish_colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = expression(log[10] ~ "Body Weight (g)"),
                  y = expression("Reproductive Rate (" * yr^-1 * ")"),
                  subtitle = "Fish Reproductive Investment by Size")

  return(gg)
}


#' Plot Growth Rate Comparison
#'
#' @title Compare growth rates across functional groups
#' @description Creates comparative plots of growth rates for fish and zooplankton
#'   groups, useful for analyzing energy budget effects on growth.
#' @details This function visualizes growth rate patterns by:
#'   - Showing size-dependent growth rates for each group
#'   - Allowing comparison between zooplankton and fish groups
#'   - Highlighting differences due to energy budget parameterisation
#'
#' @param mdl ZooMSS results object containing model outputs
#' @param groups Character string specifying which groups to plot.
#'   Options: "fish", "zooplankton", "all" (default: "all")
#' @param n_years The number of years (from the end) over which to average (default: 10)
#'
#' @return ggplot object showing growth rate vs body size by group
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#'
#' # Compare growth rates
#' all_growth <- plotGrowthComparison(results, groups = "all")
#' fish_growth <- plotGrowthComparison(results, groups = "fish")
#' zoo_growth <- plotGrowthComparison(results, groups = "zooplankton")
#' }
#'
plotGrowthComparison <- function(mdl, groups = "all", n_years = 10) {

  # Validate inputs
  assertthat::assert_that(
    is.list(mdl),
    msg = "mdl must be a list."
  )

  assertthat::assert_that(
    groups %in% c("fish", "zooplankton", "all"),
    msg = "'groups' must be one of: fish, zooplankton, all"
  )

  # Get group indices
  if (groups == "fish") {
    grp_idx <- mdl$param$fish_grps
  } else if (groups == "zooplankton") {
    grp_idx <- mdl$param$zoo_grps
  } else {
    grp_idx <- seq_len(mdl$param$ngrps)
  }

  grp_names <- mdl$param$Groups$Species[grp_idx]
  grp_colors <- mdl$param$Groups$PlotColour[grp_idx]
  names(grp_colors) <- grp_names

  # Average growth over final n_years
  growth_avg <- averageTimeSeries(mdl, var = "growth", n_years = n_years)

  # Subset to requested groups
  growth_subset <- growth_avg[grp_idx, , drop = FALSE]
  rownames(growth_subset) <- grp_names

  # Create data frame
  species_df <- tibble::as_tibble(t(growth_subset))
  species_df <- species_df %>%
    tibble::add_column("Weight" = mdl$param$w) %>%
    tidyr::pivot_longer(-"Weight", names_to = "Species", values_to = "growth") %>%
    dplyr::filter(.data$growth > 0) %>%
    dplyr::mutate(Species = factor(.data$Species, levels = grp_names))

  # Create plot
  subtitle_text <- switch(groups,
                          "fish" = "Fish Growth Rates by Size",
                          "zooplankton" = "Zooplankton Growth Rates by Size",
                          "all" = "Growth Rates by Size (All Groups)")

  gg <- ggplot2::ggplot(data = species_df,
                        mapping = ggplot2::aes(x = log10(.data$Weight),
                                               y = log10(.data$growth),
                                               colour = .data$Species)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = grp_colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = expression(log[10] ~ "Body Weight (g)"),
                  y = expression(log[10] ~ "Specific Growth Rate (" * yr^-1 * ")"),
                  subtitle = subtitle_text)

  return(gg)
}


#' Plot Defecation Rates by Prey Type
#'
#' @title Visualize prey-dependent defecation rates
#' @description Creates a plot showing how defecation rates vary with prey Carbon
#'   content, illustrating the continuous food quality scaling in the energy budget.
#' @details This function visualizes the prey-dependent defecation rates calculated
#'   using the formula:
#'   \deqn{D_{prey} = D_{high} + (D_{low} - D_{high}) \times (1 - C_{prey}/C_{max})}
#'
#'   The plot shows:
#'   - Each prey group's Carbon content and resulting defecation fraction
#'   - The continuous relationship between Carbon and defecation
#'   - The defecation range (def_high to def_low)
#'
#' @param mdl ZooMSS results object containing model outputs and parameters
#'
#' @return ggplot object showing defecation rate vs prey Carbon content
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' defecation_plot <- plotDefecation(results)
#' }
#'
plotDefecation <- function(mdl) {

 # Validate inputs
  assertthat::assert_that(
    is.list(mdl),
    msg = "mdl must be a list."
  )

  assertthat::assert_that(
    "def_by_prey" %in% names(mdl),
    msg = "Defecation data not available. Ensure model has energy budget parameters."
  )

  # Get parameters
  groups <- mdl$param$Groups
  ngrps <- mdl$param$ngrps
  Carbon <- groups$Carbon
  Carbon_max <- max(Carbon)

  # Calculate defecation for each prey (using first predator as reference)
  def_high <- groups$def_high[1]
  def_low <- groups$def_low[1]

  # Create data frame
  def_data <- data.frame(
    species = groups$Species,
    Carbon = Carbon,
    defecation = def_high + (def_low - def_high) * (1 - Carbon / Carbon_max),
    type = groups$Type
  )

  # Order by Carbon content

  def_data <- def_data[order(def_data$Carbon), ]
  def_data$species <- factor(def_data$species, levels = def_data$species)

  # Get colours
  plot_colours <- groups$PlotColour
  names(plot_colours) <- groups$Species

  # Create plot
  gg <- ggplot2::ggplot(def_data, ggplot2::aes(x = .data$Carbon, y = .data$defecation,
                                                colour = .data$species)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_line(ggplot2::aes(group = 1), colour = "grey50", linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(label = .data$species), hjust = -0.1, vjust = 0.5, size = 3) +
    ggplot2::scale_colour_manual(values = plot_colours) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = expression("Prey Carbon Content (g C g"^-1 ~ "wet weight)"),
                  y = "Defecation Fraction",
                  title = "Defecation Rate by Prey Quality",
                  subtitle = paste0("D ranges from ", def_high, " (high-C prey) to ",
                                   def_low, " (low-C prey)")) +
    ggplot2::xlim(0, max(Carbon) * 1.5) +
    ggplot2::ylim(def_high - 0.05, def_low + 0.05)

  return(gg)
}


#' Plot Energy Budget Fractions
#'
#' @title Visualize energy budget allocation across groups
#' @description Creates a stacked bar plot showing how assimilated energy is
#'   partitioned into metabolism, growth, and reproduction for each functional group.
#' @details This function visualizes the energy budget fractions:
#'   - **f_M**: Metabolic fraction (maintenance costs)
#'   - **K_growth**: Growth fraction (somatic growth)
#'   - **R_frac**: Reproduction fraction (reproductive investment, fish only)
#'
#'   The constraint f_M + K_growth + R_frac = 1 ensures all assimilated energy
#'   is accounted for.
#'
#' @param mdl ZooMSS results object containing model outputs and parameters
#'
#' @return ggplot object showing energy budget allocation by group
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' energy_plot <- plotEnergyBudget(results)
#' }
#'
plotEnergyBudget <- function(mdl) {

  # Validate inputs
  assertthat::assert_that(
    is.list(mdl),
    msg = "mdl must be a list."
  )

  # Get parameters
  groups <- mdl$param$Groups
  ngrps <- mdl$param$ngrps

  # Calculate R_frac
  R_frac <- 1 - groups$f_M - groups$K_growth

  # Create data frame for stacked bar
  energy_data <- data.frame(
    species = rep(groups$Species, 3),
    fraction = c(groups$f_M, groups$K_growth, R_frac),
    component = factor(rep(c("Metabolism (f_M)", "Growth (K)", "Reproduction (R_frac)"),
                          each = ngrps),
                      levels = c("Reproduction (R_frac)", "Growth (K)", "Metabolism (f_M)"))
  )

  # Set species order
  energy_data$species <- factor(energy_data$species, levels = groups$Species)

  # Define colours for components
  component_colours <- c("Metabolism (f_M)" = "#E41A1C",
                        "Growth (K)" = "#377EB8",
                        "Reproduction (R_frac)" = "#4DAF4A")

  # Create plot
  gg <- ggplot2::ggplot(energy_data,
                        ggplot2::aes(x = .data$species, y = .data$fraction,
                                     fill = .data$component)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = component_colours) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "Functional Group",
                  y = "Fraction of Assimilated Energy",
                  title = "Energy Budget Allocation by Group",
                  fill = "Component",
                  subtitle = "Showing f_M + K + R_frac = 1")

  return(gg)
}


#' Plot Assimilation Efficiency Matrix
#'
#' @title Visualize prey-specific assimilation efficiencies
#' @description Creates a heatmap showing the net assimilation efficiency
#'   (accounting for both defecation and growth fraction) for each predator-prey
#'   combination.
#' @details This function visualizes the assim_by_prey matrix which represents:
#'   \deqn{assim_{pred,prey} = (1 - D_{prey}) \times K_{pred}}
#'
#'   This is the fraction of ingested prey biomass that contributes to predator
#'   somatic growth, combining prey digestibility and predator growth efficiency.
#'
#' @param mdl ZooMSS results object containing model outputs and parameters
#'
#' @return ggplot object showing assimilation efficiency heatmap
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' assim_plot <- plotAssimilationMatrix(results)
#' }
#'
plotAssimilationMatrix <- function(mdl) {

  # Validate inputs
  assertthat::assert_that(
    is.list(mdl),
    msg = "mdl must be a list."
  )

  assertthat::assert_that(
    "assim_by_prey" %in% names(mdl),
    msg = "Assimilation matrix not available."
  )

  # Get data
  assim_matrix <- mdl$assim_by_prey
  species <- mdl$param$Groups$Species
  carbon <- mdl$param$Groups$Carbon
  group_type <- mdl$param$Groups$Type
  ngrps <- mdl$param$ngrps

  # Order prey by Carbon content (low to high)
  prey_order <- species[order(carbon)]
  
  # Order predators by type (Zooplankton first, then Fish) and Carbon within type
  pred_order <- species[order(group_type, carbon)]

  # Convert to long format
  assim_long <- expand.grid(
    predator = species,
    prey = species
  )
  assim_long$efficiency <- as.vector(assim_matrix)

  # Set factor levels based on Carbon ordering
  assim_long$predator <- factor(assim_long$predator, levels = pred_order)
  assim_long$prey <- factor(assim_long$prey, levels = prey_order)

  # Create heatmap
  gg <- ggplot2::ggplot(assim_long,
                        ggplot2::aes(x = .data$prey, y = .data$predator,
                                     fill = .data$efficiency)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = round(.data$efficiency, 3)),
                       size = 2.5, colour = "white") +
    ggplot2::scale_fill_viridis_c(option = "plasma", limits = c(0, 0.5)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "Prey Group (ordered by Carbon content)",
                  y = "Predator Group",
                  fill = "Growth\nEfficiency",
                  title = "Prey-Specific Growth Efficiency",
                  subtitle = expression("(1 - D"[prey]*") × K"[predator]))

  return(gg)
}


