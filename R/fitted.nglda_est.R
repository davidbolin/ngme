#' @title Fitted values.
#'
#' @description Calculates fitted values.
#'
#' @param object A fitted object by calling \code{"summary.nglda_est"}.
#' @param type A character string for the type of fitted values;
#'   \code{"marginal"} for fixed effects.
#' @param ... Additional arguments; none used currently.
#'
#' @details Subject-specific fitted values can be obtained using \code{"nglda_predict"}.
#' @seealso \code{\link{nglda_est}}, \code{\link{fitted}}
#'
#' @examples
#'   \dontrun{
#'   fit <- nglda_est(...)
#'   fitted(fit)
#'   }
#'

 fitted.nglda_est <- function(object, type = "marginal", ...){

   if(type == "marginal"){
     marg_fit <- as.numeric(object$x_fixed_f %*% object$fixed_est)
     marg_fit
   }

 }
