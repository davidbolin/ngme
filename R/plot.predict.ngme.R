
#' @title Prediction plots.
#'
#' @description Plots the predicted values for a specific subject.
#'
#' @param object A fitted object returned by the \code{"predict.ngme"} function.
#' @param id A numerical value or character string for ID of the subject
#'   for whom the plot will be generated.
#' @param col_m A character value for defining the colour of prediction mean 
#'   or median.
#' @param col_c A character value for defining the colour of prediction intervals.
#' @param col_p A character value for defining the colour of observed data.
#' @param ... Additional arguments; none used currently.
#'
#' @seealso \code{\link{predict.ngme}}
#'
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   pred <- predict(fit, ...)
#'   plot(pred, 1)
#'   }
#'

plot.predict.ngme <- function(object, 
                              id,
                              col_m = "black",
                              col_c = "red",
                              col_p = "black",
                              ...){

  Y         <- object$Y
  locs      <- object$predictions$locs
  id_list   <- object$id_list
  X.summary <- object$predictions$X.summary

  pInd     <- which(object$id == id)
  pInd_all <- which(id_list == id)

  locs_i <- locs[[pInd_all]]

  if(length(locs_i) > 1){#if the subject has more than one measurement -- produce a plot

  mean_i <- X.summary[[pInd]]$Mean
  llim_i <- X.summary[[pInd]]$quantiles[[1]]$field
  ulim_i <- X.summary[[pInd]]$quantiles[[2]]$field
  y_i    <- Y[[pInd_all]]

  x_range <- range(locs_i)
  y_range <- range(c(mean_i, llim_i, ulim_i, y_i))
  y_range_inc <- diff(y_range)/100

  Time    <- locs_i
  Outcome <- mean_i
  
    plot(Time, 
         Outcome, 
         type = "l",
         col = col_m,
         ylim = c(y_range[1] - y_range_inc, y_range[2] + y_range_inc),
         ...
         )
    lines(locs_i,  llim_i, col = col_c, ...)
    lines(locs_i,  ulim_i, col = col_c, ...)
    points(locs_i, y_i,    col = col_p, ...)
        
  }else{#the subject has one measurement - do not produce a plot
    print("The subject has 1 measurement, no plot is produced")
  }

}

