#' Queensland mortality data
#'
#' A dataset containing number of deaths, population at risk and
#' ages for residents of Queensland, Australia from 1979 to 2020
#'
#' The age group 100 also includes people who died aged older than 100.
#' The data come from the Australian Human Mortality Database (\url{https://aushd.org}).
#'
#' @format A `data.frame` containing the following fields:
#' \describe{
#' \item{year}{integer, year of records}
#' \item{age}{integer, age at death, in years}
#' \item{sex}{factor differentiating the sexes}
#' \item{deaths}{integer, number of deaths recorded}
#' \item{population}{numeric population recorded at 30 June each year}
#' }
#' @source
#' Australian Human Mortality Database
"qld_mortality"
