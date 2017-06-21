#' @title Print function for the function \code{intervals}.
#'
#' @description Prints the output produced by \code{intervals} function concisely.
#'
#' @param object A fitted object by calling \code{"intervals"}.
#' @param ... Additional arguments; none used currently.

#' @seealso \code{\link{intervals}}
#'
#' @examples
#'   \dontrun{
#'   fit <- nglda_est(...)
#'   intervals(fit)
#'   }
#'

 print.intervals <- function(object, ...){

   if(object$type == "fixed"){
     cat(paste0(object$level*100, "% confidence intervals for the fixed effects"))

     cat("\n"); cat("\n")

     print(object$result)
   }

 }
