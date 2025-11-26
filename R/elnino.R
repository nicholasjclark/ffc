#' El Niño Sea Surface Temperature Data
#'
#' Monthly sea surface temperature measurements from El Niño regions 1 and 2
#' (OISST) spanning 1982-2018. The data captures the characteristic seasonal
#' patterns and inter-annual variability associated with El Niño Southern
#' Oscillation (ENSO) events.
#'
#' @format A tsibble with 444 rows and 3 columns:
#' \describe{
#'   \item{year}{Integer year from 1982 to 2018}
#'   \item{month}{Integer month from 1 to 12}
#'   \item{temperature}{Sea surface temperature in degrees Celsius}
#' }
#'
#' @details
#' The dataset is structured as a tsibble with year as the index and month as
#' the key variable. This format is ideal for functional forecasting where
#' seasonal patterns (monthly temperature curves) evolve over years.
#'
#' El Niño regions 1 and 2 are located in the eastern tropical Pacific Ocean
#' and are critical for monitoring ENSO events. Temperature anomalies in these
#' regions are strong predictors of global climate patterns.
#'
#' @source
#' Data originally from \code{rainbow::ElNino_OISST_region_1and2}.
#' OISST (Optimally Interpolated Sea Surface Temperature) data from NOAA.
#'
#' @references
#' Hyndman, R.J. and Shang, H.L. (2010) Rainbow plots, bagplots, and boxplots
#' for functional data. Journal of Computational and Graphical Statistics,
#' 19(1), 29-45.
#'
#' @examples
#' data(elnino_sst)
#' head(elnino_sst)
#'
#' # View seasonal patterns
#' library(ggplot2)
#' ggplot(elnino_sst, aes(x = month, y = temperature, group = year, color = year)) +
#'   geom_line(alpha = 0.5) +
#'   scale_color_viridis_c()
#'
"elnino_sst"
