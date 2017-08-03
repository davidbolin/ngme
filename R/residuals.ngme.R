#' @title Residual values.
#'
#' @description Calculates residuals.
#'
#' @param object A fitted object returned by the \code{"summary.ngme"} function.
#' @param type A character string for the type of fitted values;
#'   \code{"marginal"} for fixed effects.
#' @param ... Additional arguments; none used currently.
#'
#' @details Subject-specific residaul values are not available currently.
#' @seealso \code{\link{ngme}}, \code{\link{fitted.ngme}}
#'
#'
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   fitted(fit)
#'   }
#'
residuals.ngme <- function(object, type = "marginal", ...){

  if(type == "marginal"){
    resid_marg <- as.numeric(unlist(object$Y) - fitted(object, type = "marginal", ...))
    resid_marg
  }

}
