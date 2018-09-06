#'
#' Orthodont data
#'
#' A data-set that contains dental data for 27 children at ages 8, 10, 12 and 14
#'
#' @format A data frame with 108 rows and 5 columns:
#' \describe{
#' \item{distance}{outcome variable}
#' \item{age}{age of the child}
#' \item{Subject}{id as character}
#' \item{Sex}{sex of the child}
#' \item{id}{id as numeric}
#' }
#' @examples
#'   \dontrun{
#' data(Orthodont)
#' set.seed(123)
#' NormalMVD <- ngme(                   fixed       = distance ~ age,
#'                                      random      = ~ 1|id,
#'                                      data        = Orthodont,
#'                                      reffects    = "Normal",
#'                                      use.process = F,
#'                                      silent      = F,
#'                                      nIter       = 2,
#'                                      controls    = list(estimate.fisher = FALSE,
#'                                                         subsample.type  = 0))
#' print(NormalMVD$fixed_est)
#' }
"Orthodont"
