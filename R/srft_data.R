#'
#'srft data
#'
#' A data-set that contains for estimated glomerular filtration rate and
#' a number of explanatory variables for 22,910 primary care patients
#'
#' @format A data frame with 392870 rows and 6 variables:
#' \describe{
#' \item{id}{identification number}
#' \item{egfr}{estimated glomerular filtration rate,
#' in mL/min per 1.73\eqn{m^2} of body surface area)}
#' \item{sex}{sex of the patient; 0: male, 1: female}
#' \item{bage}{baseline age, in yars}
#' \item{fu}{follow-up time; in years}
#' \item{pwl}{a piece-wise linear term defined as max(age - 56.5, 0),
#' where age is age at measurement}
#' }
#'@source \url{https://doi.org/10.1093/biostatistics/kxu053}
"srft_data"
