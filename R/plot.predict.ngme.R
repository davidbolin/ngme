
#' @title Prediction plots.
#'
#' @description Plots the predicted values for a specific subject.
#'
#' @param object A fitted object returned by the \code{"predict.ngme"} function.
#' @param id A numerical value or character string for ID of the subject
#'   for whom the plot will be generated.
#' @param plot_excursions A logical for plotting excursions.
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
                              id = NULL){
  
  if(is.null(id) == TRUE) id <- names(object$predictions$locs) %>% as.numeric()
  
  Y         <- object$Y
  locs      <- object$predictions$locs
  id_list   <- object$id_list
  X.summary <- object$predictions$X.summary
  
  pInd     <- which(names(locs) %>% as.numeric() %in% id)
  pInd_all <- which(id_list %in% id)
  
  data_plot <- data.frame(id = rep(id, lapply(pInd, function(i) locs[[i]]) %>% lapply(length) %>% unlist()))
  data_plot$locs <- lapply(pInd, function(i) locs[[i]]) %>% unlist()  
  data_plot$mean <- lapply(pInd, function(i) X.summary[[i]]$Mean) %>% unlist()
  data_plot$llim <- lapply(pInd, function(i) X.summary[[i]]$quantiles[[1]]$field) %>% unlist()
  data_plot$ulim <- lapply(pInd, function(i) X.summary[[i]]$quantiles[[2]]$field) %>% unlist()
  data_plot$y    <- lapply(pInd_all, function(i) Y[[i]]) %>% unlist()
  
  if(is.null(object$predictions$Xderivative) == TRUE){
    g <- ggplot(data_plot, aes(x = locs, y = y, group = id))
    g + geom_point() + facet_wrap(~ id, scales = "free", ncol = floor(sqrt(length(id)))) + 
      geom_line(aes(y = mean)) + 
      geom_line(aes(y = ulim), linetype = "dashed") + 
      geom_line(aes(y = llim), linetype = "dashed") + 
      labs(x = "Time", y = "Outcome")
  }else{
    data$probability <- lapply(pInd, function(i) object$predict$Xderivative.summary[[i]]$excursions$P) %>% unlist()
    
    g1 <- ggplot(data_plot, aes(x = locs, y = y, group = id)) + 
      geom_point() + facet_wrap(~ id, scales = "free", ncol = 1) + 
      geom_line(aes(y = mean)) + 
      geom_line(aes(y = ulim), linetype = "dashed") + 
      geom_line(aes(y = llim), linetype = "dashed") + 
      labs(x = "Time", y = "Outcome")
    g2 <- ggplot(data_plot, aes(x = locs, y = probability, group = id)) + 
      geom_point() + facet_wrap(~ id, scales = "free", ncol = 1) + 
      labs(x = "Time", y = "Probability")
    grid.arrange(g1, g2, ncol = 2)
  }
  
}

