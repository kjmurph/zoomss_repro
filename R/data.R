#' ZooMSS Functional Groups Data
#'
#' @title Default functional groups for the ZooMSS model
#' @description A dataset containing the biological parameters for different
#'   functional groups used in the ZooMSS size-structured marine ecosystem model.
#'   These represent various taxa from flagellates to large fish, each defined
#'   by their feeding behavior, size ranges, and physiological parameters.
#'
#' @format A data frame with 12 rows (functional groups) and 24 columns:
#' \describe{
#'   \item{Species}{Character. Name of the functional group/taxa}
#'   \item{Type}{Character. Broad category (Zooplankton or Fish)}
#'   \item{FeedType}{Character. Feeding strategy (Heterotroph, FilterFeeder, Omnivore, Carnivore)}
#'   \item{Prop}{Numeric. Initial proportion of total biomass}
#'   \item{W0}{Numeric. Log10 minimum body weight (g) for the group}
#'   \item{Wmax}{Numeric. Log10 maximum body weight (g) for the group}
#'   \item{Wmat}{Numeric. Log10 maturation body weight (g)}
#'   \item{SearchCoef}{Numeric. Search coefficient for predation interactions}
#'   \item{SearchExp}{Numeric. Search exponent for predation scaling}
#'   \item{PPMRscale}{Numeric. Predator-prey mass ratio scaling parameter}
#'   \item{PPMR}{Numeric. Predator-prey mass ratio (for fish groups)}
#'   \item{FeedWidth}{Numeric. Feeding kernel width parameter}
#'   \item{Carbon}{Numeric. Carbon content (g C / g wet weight) of the group when consumed as prey}
#'   \item{def_high}{Numeric. Defecation fraction for high-quality (high Carbon) prey (typically 0.30)}
#'   \item{def_low}{Numeric. Defecation fraction for low-quality (low Carbon) prey (typically 0.50)}
#'   \item{f_M}{Numeric. Metabolic fraction of assimilated energy (typically 0.50)}
#'   \item{K_growth}{Numeric. Growth fraction of assimilated energy (0.50 for zooplankton, 0.36 for fish)}
#'   \item{repro_eff}{Numeric. Reproductive efficiency - egg-to-recruit survival fraction (fish only)}
#'   \item{repro_on}{Integer. Flag to enable reproduction (0 = off, 1 = on; fish only)}
#'   \item{mat_ogive_slope}{Numeric. Steepness of maturity ogive function (typically 10)}
#'   \item{Fmort}{Numeric. Fishing mortality rate}
#'   \item{Fmort_W0}{Numeric. Log10 minimum weight for fishing mortality}
#'   \item{Fmort_Wmax}{Numeric. Log10 maximum weight for fishing mortality}
#'   \item{PlotColour}{Character. Color code for plotting the functional group}
#' }
#'
#' @details The GroupInputs dataset defines 12 functional groups spanning from
#'   small microzooplankton (flagellates, ciliates) through various mesozooplankton
#'   groups (copepods, euphausiids, chaetognaths) to gelatinous zooplankton (salps, jellyfish)
#'   and three fish size classes (small, medium, large). Each group is characterized by:
#'
#'   - **Size ranges**: W0 to Wmax define the body size spectrum
#'   - **Feeding behavior**: Different strategies for resource acquisition
#'   - **Interaction parameters**: Search rates and predator-prey relationships
#'   - **Energy budget parameters**: Defecation, metabolism, and growth fractions following DBPM theory
#'   - **Reproduction parameters**: For fish groups, explicit reproduction with maturity ogive
#'
#'   The energy budget follows a two-stage partitioning:
#'   1. Ingestion -> Defecation (D) + Assimilation (A = 1 - D)
#'   2. Assimilation -> Metabolism (f_M) + Growth (K_growth) + Reproduction (R_frac = 1 - f_M - K_growth)
#'
#'   Defecation is calculated using continuous scaling based on prey Carbon content:
#'   def_effective = def_high + (def_low - def_high) * (1 - Carbon_prey / Carbon_max)
#'
#'   These parameters are based on marine ecological literature, DBPM theory,
#'   and represent typical values for temperate marine ecosystems.
#'
#' @source Marine ecological literature, DBPM (Blanchard et al.), and ZooMSS model development
#' @family ZooMSS-data
#' @examples
#' data(GroupInputs)
#' head(GroupInputs)
#'
#' # View size ranges across groups
#' plot(GroupInputs$W0, GroupInputs$Wmax,
#'      col = GroupInputs$PlotColour,
#'      xlab = "Log10 Min Weight", ylab = "Log10 Max Weight")
#' text(GroupInputs$W0, GroupInputs$Wmax, GroupInputs$Species, pos = 3, cex = 0.7)
"GroupInputs"
