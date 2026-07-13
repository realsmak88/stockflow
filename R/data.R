###############################################################
#
# stockflow - Documentation of datasets
#
# FICTITIOUS DATASETS. This data is entirely SIMULATED
# (see data-raw/simulate_datasets.R): it reproduces the
# STRUCTURE of fisheries monitoring data (columns, types,
# factor levels, plausible ranges) but contains
# NO real observations. Sites and locations are fictitious. It
# is used solely to illustrate and test the package's functions.
#
# Three illustrative taxa: Cymbium spp. (gastropod),
# Octopus vulgaris (octopus), Penaeus notialis (pink shrimp).
#
###############################################################


#' Individual biological data - Cymbium spp.
#'
#' Individual measurements (size, weight), sex and maturity stage
#' simulated for landings of \emph{Cymbium} spp. (fictitious data).
#'
#' @format A data.frame with 3130 rows and 10 variables:
#' \describe{
#'   \item{landing_site}{Landing site (factor).}
#'   \item{Date}{Sampling date (Date).}
#'   \item{fishing_zone}{Declared fishing ground (factor).}
#'   \item{gear}{Fishing gear (factor).}
#'   \item{season}{Hydrological season (factor).}
#'   \item{species}{\emph{Cymbium} species (factor, 4 levels).}
#'   \item{LCQ}{Shell length (cm).}
#'   \item{weight}{Individual weight (g).}
#'   \item{sex}{Sex: M / F (factor).}
#'   \item{maturity}{Maturity stage: IM (immature), M (mature),
#'     P (spawning) (factor).}
#' }
#' @source Fictitious simulated data (see \code{data-raw/simulate_datasets.R}).
"cymbium_bio"


#' Individual biological data - Octopus vulgaris
#'
#' Individual measurements, sex and maturity of the common octopus
#' \emph{Octopus vulgaris} (fictitious data).
#'
#' @format A data.frame with 1682 rows and 12 variables:
#' \describe{
#'   \item{landing_site}{Landing site (factor).}
#'   \item{Date}{Sampling date (Date).}
#'   \item{fishing_zone}{Fishing ground (factor).}
#'   \item{gear}{Fishing gear (factor).}
#'   \item{season}{Hydrological season (factor).}
#'   \item{species}{Species (factor).}
#'   \item{LM}{Mantle length (cm).}
#'   \item{weight}{Individual weight (g).}
#'   \item{sex}{Sex: M / F (factor).}
#'   \item{maturity}{Maturity stage (factor, scale 1-4).}
#'   \item{total_weight_kg}{Total weight of the sampled catch (kg).}
#'   \item{sample_weight_kg}{Weight of the measured sample (kg).}
#' }
#' @source Fictitious simulated data (see \code{data-raw/simulate_datasets.R}).
"octopus_bio"


#' Individual biological data - Penaeus notialis
#'
#' Individual measurements, sex, maturity and gonad weight of the
#' pink shrimp \emph{Penaeus notialis} (fictitious data).
#'
#' @format A data.frame with 6286 rows and 13 variables:
#' \describe{
#'   \item{landing_site}{Landing site (factor).}
#'   \item{Date}{Sampling date (Date).}
#'   \item{fishing_zone}{Fishing ground (factor).}
#'   \item{gear}{Fishing gear (factor).}
#'   \item{season}{Hydrological season (factor).}
#'   \item{species}{Species (factor).}
#'   \item{LCT}{Cephalothoracic length (mm).}
#'   \item{weight}{Individual weight (g).}
#'   \item{sex}{Sex: M / F (factor).}
#'   \item{maturity}{Maturity stage (factor, scale 1-4).}
#'   \item{gonad_weight}{Gonad weight (g).}
#'   \item{total_weight_kg}{Total weight of the sampled catch (kg).}
#'   \item{sample_weight_kg}{Weight of the measured sample (kg).}
#' }
#' @source Fictitious simulated data (see \code{data-raw/simulate_datasets.R}).
"penaeus_bio"


#' Size frequencies - Cymbium spp.
#'
#' Size measurements from the sampling of catches of
#' \emph{Cymbium} spp. (input structure for LFQ analyses).
#'
#' @format A data.frame with 8243 rows and 9 variables (same descriptive
#'   columns as \code{\link{cymbium_bio}}, with \code{LCQ} for the shell
#'   length and \code{sex}).
#' @seealso \code{\link{prepare_tropfish}} to build an lfq object.
#' @source Fictitious simulated data (see \code{data-raw/simulate_datasets.R}).
"cymbium_freq"


#' Size frequencies - Octopus vulgaris
#'
#' Size measurements from the sampling of catches of the octopus
#' \emph{Octopus vulgaris}.
#'
#' @format A data.frame with 1689 rows and 9 variables (\code{LM} for the
#'   mantle length).
#' @source Fictitious simulated data (see \code{data-raw/simulate_datasets.R}).
"octopus_freq"


#' Size frequencies - Penaeus notialis
#'
#' Size measurements from the sampling of catches of the
#' pink shrimp \emph{Penaeus notialis}.
#'
#' @format A data.frame with 31043 rows and 8 variables (\code{LCT} for the
#'   cephalothoracic length).
#' @source Fictitious simulated data (see \code{data-raw/simulate_datasets.R}).
"penaeus_freq"


#' Annual total catches by group
#'
#' Annual catch series (tonnes) for the three target groups.
#'
#' @format A data.frame with 55 rows and 4 variables:
#' \describe{
#'   \item{year}{year.}
#'   \item{Cymbium}{Catches of \emph{Cymbium} spp. (t).}
#'   \item{Octopus}{Catches of \emph{Octopus vulgaris} (t).}
#'   \item{Penaeus}{Catches of \emph{Penaeus notialis} (t).}
#' }
#' @source Fictitious simulated data (see \code{data-raw/simulate_datasets.R}).
"total_catches"


#' Annual CPUE indices by group
#'
#' Annual catch per unit effort (abundance index) for the
#' three target groups.
#'
#' @format A data.frame with 47 rows and 4 variables:
#' \describe{
#'   \item{year}{year.}
#'   \item{Cymbium}{CPUE \emph{Cymbium} spp.}
#'   \item{Octopus}{CPUE \emph{Octopus vulgaris}.}
#'   \item{Penaeus}{CPUE \emph{Penaeus notialis}.}
#' }
#' @source Fictitious simulated data (see \code{data-raw/simulate_datasets.R}).
"annual_cpue"
