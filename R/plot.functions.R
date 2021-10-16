#' @title Trace plots.
#'
#' @description Plots the parameter estimates at successive iterations
#'    of the stochastic algorithm.
#'
#' @param object A fitted object returned by the \code{"ngme"} function.
#' @param param A character string for which parameter set the traces
#'   are to be plotted;
#'   \code{"fixed"} for fixed effects,
#'   \code{"random"} for time-invariant random effects,
#'   \code{"process"} for stochastic process,
#'   \code{"error"} for error term.
#' @param n.exclude A numeric value for excluding a portion of the
#'    trajectory (from left) for plotting.
#' @param index When plotting fixed effects, the indices of the effects to plot.

#' @param ... Additional arguments; none used currently.
#'
#' @seealso \code{\link{ngme}}
#'
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   fitted(fit)
#'   }
#' @export
#' @method plot ngme

plot.ngme <- function(object, param = "fixed", n.exclude = 0, index = NULL, ...){

  # if process was not specified in the model
  if(param == "process" && object$use_process == FALSE){
    stop("No process was included in the model fit")
  }

  if(param == "fixed"){

    if(n.exclude == 0){
      if(is.null(index) == TRUE){
        plot.ts(object$fixed_est_vec, main = "Fixed effects", xlab = "Iteration")
      }else{
        plot.ts(object$fixed_est_vec[, index], main = "Fixed effects", xlab = "Iteration")
        }
    }else{
      if(is.null(index) == TRUE){
        plot.ts(as.matrix(object$fixed_est_vec)[-c(1:n.exclude), ], main = "Fixed effects", xlab = "Iteration")
      }else{
        plot.ts(as.matrix(object$fixed_est_vec)[-c(1:n.exclude), index], main = "Fixed effects", xlab = "Iteration")
      }
    }

  }else if(param == "random"){

    # preparing the ranef_Sigma_vec: excluding duplicated columns and adding column names
    ncol_Sigma        <- ncol(summary(object)$random_results$Sigma_est)
    if(ncol_Sigma > 1){
      cols_omit         <- unlist(lapply(1:(ncol_Sigma-1), function(i) (ncol_Sigma * i + 1):(ncol_Sigma * i + i)))
      ranef_results_vec <- object$ranef_Sigma_vec
      ranef_results_vec <- ranef_results_vec[, -cols_omit]

      colnames_Sigma_est          <- colnames(summary(object)$random_results$Sigma_est)
      #colnames_Sigma_est          <- substr(colnames_Sigma_est, 1, 3)
      colnames_Sigma_est_ext      <- paste(rep(colnames_Sigma_est,each = ncol_Sigma), colnames_Sigma_est, sep = "_")
      colnames_Sigma_est_ext      <- colnames_Sigma_est_ext[-cols_omit]
      colnames(ranef_results_vec) <- colnames_Sigma_est_ext
    }else{
      ranef_results_vec <- object$ranef_Sigma_vec
      colnames(ranef_results_vec) <- colnames(summary(object)$random_results$Sigma_est)
    }

    if(summary(object)$random_distr %in% ("Normal")){

      if(n.exclude == 0){
        plot.ts(ranef_results_vec, main = "Random effects", xlab = "Iteration")
      }else{
        plot.ts(as.matrix(ranef_results_vec)[-c(1:n.exclude), ],
                main = "Random effects", xlab = "Iteration", ylab = colnames(ranef_results_vec))
      }

    }else if(summary(object)$random_distr %in% c("NIG", "tdist")){

      ranef_results_vec <- cbind(ranef_results_vec, object$ranef_mu_vec, object$ranef_nu_vec)
      colnames(ranef_results_vec)[(ncol(ranef_results_vec)-ncol_Sigma ) : (ncol(ranef_results_vec))] <-
        c(paste("mu", colnames(summary(object)$random_results$Sigma_est), sep = "_"), "nu")

      if(n.exclude == 0){
        plot.ts(ranef_results_vec, main = "Random effects", xlab = "Iteration")
      }else{
        plot.ts(as.matrix(ranef_results_vec)[-c(1:n.exclude), ], main = "Random effects", xlab = "Iteration")
      }

    }

  }else if(param == "process"){

    if(object$process_distr %in% c("Normal", "CH")){

      if(object$operator_type == "fd2"){
        process_results_vec <- matrix(object$operator_tau_vec, ncol = 1)
        colnames(process_results_vec) <- "tau"
      }else if(object$operator_type == "matern" || object$operator_type == "exponential"){
        process_results_vec <- cbind(object$operator_tau_vec,
                                     object$operator_kappa_vec)
        colnames(process_results_vec) <- c("tau", "kappa")
      }

      if(n.exclude == 0){
        plot.ts(process_results_vec, main = "Process", xlab = "Iteration")
      }else{
        plot.ts(process_results_vec[-c(1:n.exclude), ], main = "Process", xlab = "Iteration")
      }

    }else if(object$process_distr %in% c("NIG", "GAL")){

      if(object$operator_type == "fd2"){
        process_results_vec <- cbind(object$operator_tau_vec,
                                     object$process_mu_vec,
                                     object$process_nu_vec)
        colnames(process_results_vec) <- c("tau", "mu", "nu")
      }else if(object$operator_type == "matern"){
        process_results_vec <- cbind(object$operator_tau_vec,
                                     object$operator_kappa_vec,
                                     object$process_mu_vec,
                                     object$process_nu_vec)
        colnames(process_results_vec) <- c("tau", "kappa","mu", "nu")
      }

      if(n.exclude == 0){
        plot.ts(process_results_vec, main = "Process", xlab = "Iteration")
      }else{
        plot.ts(process_results_vec[-c(1:n.exclude), ], main = "Process", xlab = "Iteration")
      }
    }

  }else if(param == "error"){

    if(object$error_distr == "Normal"){

      if(n.exclude == 0){
        plot(object$meas_error_sigma_vec, main = "Measurement error",
             xlab = "Iteration", ylab = "sigma", type = "l")
      }else{
        plot(as.matrix(object$meas_error_sigma_vec)[-c(1:n.exclude), ], main = "Measurement error",
             xlab = "Iteration", ylab = "sigma", type = "l")
      }


    }else if(object$error_distr %in% c("NIG", "tdist")){
      meas_error_result_vec <- cbind(object$meas_error_sigma_vec,
                                   object$meas_error_nu_vec)
      colnames(meas_error_result_vec) <- c("sigma", "nu")

      if(n.exclude == 0){
        plot.ts(meas_error_result_vec, main = "Measurement error", xlab = "Iteration")
      }else{
        plot.ts(as.matrix(meas_error_result_vec)[-c(1:n.exclude), ], main = "Measurement error", xlab = "Iteration")
      }

    }
  }
}



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
#' @export
#' @method plot predict.ngme

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

