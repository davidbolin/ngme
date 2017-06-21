
#' @title Print function for \code{"nglda_predict"} objects.
#'
#' @description Prints the output produced by \code{"nglda_predict"} concisely.
#'
#' @param object A fitted object by calling \code{"nglda_predict"}.
#' @param ... Additional arguments; none used currently.
#'
#' @details MAE, RMSE, CRPS, Interval coverage, Interval width
#'
#' @seealso \code{\link{nglda_predict}}, \code{\link{print}}
#' @examples
#'   \dontrun{
#'   pred <- nglda_predict(...)
#'   pred
#'   }

print.nglda_predict <- function(object, ...){

  cat("Call:\n")
  print(object$call)

  cat("\n")

  if(object$type == "Filter"){
    cat("Filtering for subjects with the following ID(s) (", length(object$id), " of them in total)", ":\n", sep = "")
    if(length(object$id) > 10){
      cat(head(object$id, 10), "...\n")
    }else{
      cat(object$id, "\n")
    }
  }else{
    cat("Smoothing for subjects with the following ID(s) (", length(object$id), " of them in total)", ":\n", sep = "")
    if(length(object$id) > 10){
      cat(head(object$id, 10), "...\n")
    }else{
      cat(object$id, "\n")
    }
  }

  cat("\n")

  cat("MAE (SD) = ", object$mean.mae," (", object$std.mae, ")\n", sep = "")
  cat("RMSE (SD) = ",  object$mean.rmse," (", object$std.rmse,")\n", sep = "")
  cat("CRPS (SD) = ", object$mean.crps," (", object$std.crps, ")\n", sep = "")
  cat("Interval coverage (SD) = ", object$coverage.mean, " (", object$coverage.std,")\n", sep = "")
  cat("Interval width (SD) = ", object$mean.int.width, " (", object$std.int.width,")\n", sep = "")


}
