#' @title Fitted values.
#'
#' @description Calculates fitted values.
#'
#' @param object A fitted object by calling \code{"ngme"} function.
#' @param type A character string for the type of fitted values;
#'   \code{"marginal"} for fixed effects.
#' @param ... Additional arguments; none used currently.
#'
#' @details Subject-specific fitted values can be obtained using \code{"predict.ngme"}.
#' @seealso \code{\link{ngme}}, \code{\link{fitted}}
#'
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   fitted(fit)
#'   }
#'

 fitted.ngme <- function(object, type = "marginal", ...){

   if(type == "marginal"){
     marg_fit <- as.numeric(object$x_fixed_f %*% object$fixed_est)
     marg_fit
   }

 }
