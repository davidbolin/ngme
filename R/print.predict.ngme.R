
#' @title Print function for \code{"predict.ngme"} objects.
#'
#' @description Prints the output contained in the \code{"predict.ngme"} objects concisely.
#'
#' @param object A fitted object by calling \code{"predict.ngme"}.
#' @param ... Additional arguments; none used currently.
#'
#' @details MAE, RMSE, CRPS, Interval coverage, Interval width
#'
#' @seealso \code{\link{predict.ngme}}, \code{\link{print}}
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   pred <- predict.ngme(fit, ...)
#'   pred
#'   }

print.predict.ngme <- function(object, ...){

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

  cat("Mean predictor MAE (SD) = ", object$mean.mae.mean.predictor," (", object$std.mae.mean.predictor, ")\n", sep = "")
  cat("Median predictor MAE (SD) = ", object$mean.mae.median.predictor," (", object$std.mae.median.predictor, ")\n", sep = "")
  cat("Mean predictor RMSE (SD) = ",  object$mean.rmse.mean.predictor," (", object$std.rmse.mean.predictor,")\n", sep = "")
  cat("Median predictor RMSE (SD) = ",  object$mean.rmse.median.predictor," (", object$std.rmse.median.predictor,")\n", sep = "")
  cat("CRPS (SD) = ", object$mean.crps," (", object$std.crps, ")\n", sep = "")
  cat("Interval coverage (SD) = ", object$coverage.mean, " (", object$coverage.std,")\n", sep = "")
  cat("Interval width (SD) = ", object$mean.int.width, " (", object$std.int.width,")\n", sep = "")

}
