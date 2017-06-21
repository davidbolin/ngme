#' @title Prediction plots.
#'
#' @description Plots the predicted values for a specific subject.
#'
#' @param object A fitted object by calling \code{"nglda_est"}.
#' @param id A numerical value or character string for ID of the subject
#'   for whom the plot will be generated.
#' @param control A list of control variables.
#'   \itemize{
#'   \item \code{"xlim_x"} A numerical value to control the range of x-axis.
#'   \item \code{"ylim_c"} A numerical value to control the range of y-axis.
#'   }
#' @param ... Additional arguments; none used currently.
#'
#' @seealso \code{\link{nglda_predict}}
#'
#' @examples
#'   \dontrun{
#'   pred <- nglda_predict(...)
#'   plot(fit, 1)
#'   }
#'

plot.nglda_predict <- function(object, id,
                               control = list(xlim_x = 0.1,
                                              ylim_c = 0.01),
                               ...){

  Y         <- object$Y
  locs      <- object$predictions$locs
  id_list   <- object$id_list
  X.summary <- object$predictions$X.summary

  pInd     <- which(object$id == id)
  pInd_all <- which(id_list == id)

  locs_i <- locs[[pInd_all]]

  if(length(locs_i) > 1){

  mean_i <- X.summary[[pInd]]$Mean
  llim_i <- X.summary[[pInd]]$quantiles[[1]]$field
  ulim_i <- X.summary[[pInd]]$quantiles[[2]]$field
  y_i    <- Y[[pInd_all]]

  x_range <- range(locs_i)
  y_range <- range(c(mean_i, llim_i, ulim_i, y_i))

  plot(locs_i,   mean_i, type = "l",
       xlim = c(x_range[1], x_range[2] + control$xlim_x),
       ylim = c(y_range[1] - control$ylim_c, y_range[2] + control$ylim_c),
       xlab = "Time", ylab = "Outcome", ...)
  lines(locs_i,  llim_i, col = 2, ...)
  lines(locs_i,  ulim_i, col = 2, ...)
  points(locs_i, y_i, ...)
  }else{
    print("The subject has 1 measurement, no plot is produced")
  }

}

