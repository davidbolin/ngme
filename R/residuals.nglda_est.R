#' @title Residual values.
#'
#' @description Calculates residual values.
#'
#' @param object A fitted object by calling \code{"summary.nglda_est"}.
#' @param type A character string for the type of fitted values;
#'   \code{"marginal"} for fixed effects.
#' @param ... Additional arguments; none used currently.
#'
#' @details Subject-specific residaul values are not available currently.
#' @seealso \code{\link{nglda_est}}, \code{\link{fitted.nglda_est}}
#'
#'
#' @examples
#'   \dontrun{
#'   fit <- nglda_est(...)
#'   fitted(fit)
#'   }
#'
residuals.nglda_est <- function(object, type = "marginal", ...){

  if(type == "marginal"){
    resid_marg <- as.numeric(unlist(object$Y) - fitted(object, type = "marginal", ...))
    resid_marg
  }

}
