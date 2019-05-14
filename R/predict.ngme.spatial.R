
#' @title Prediction.
#'
#' @description Obtains predicted values based on filtering and smoothing distributions.
#'
#' @param object A fitted object obtained by calling \code{"ngme"}.
#' @param id A numeric vector containing the ID's of the replicates for whom
#'   predictions are to be obtained. Default is set to \code{"NULL"}
#'   indicating perform predictions for all replicates.
#' @param type A character string for the type of prediction: \code{"Smoothing"} gives spatial prediction
#' based on all available data, \code{"LOOCV"} gives leave-one-out crossvalidation where also a number of 
#' accuracy measures are calculated. 
#' @param quantiles A two-element vector that contains the quantiles
#'   of the predictions to be calculated.
#' @param controls A list of control variables.
#'  \itemize{
#'  \item \code{"return.samples"} A logical variable for returning the
#'    Monte Carlo samples used to compute the predictions; \code{"TRUE"} indicates
#'    return, \code{"FALSE"} do not return.
#'  \item \code{"excursions"} A list of excursion probabilities to compute.
#'     Each list should contain:
#'    \itemize{
#'    \item \code{"type"} Type of excursion for indicating calculation of the probabilities, 
#'     greater than or less than a threshol. User should specify "greater than" by \code{">"} 
#'     "less than" by \code{"<"}.
#'    \item \code{"level"} A numeric value for the threshold compute excursion probability for.
#'    \item \code{"process"} A character string that for which component of the model 
#'     the excursion probabilities to be calculated.
#'    \itemize{
#'    \item \eqn{x\alpha + dU + W + Z} with \eqn{x \alpha} being fixed effects,
#'    \item \eqn{dU} random effects and \eqn{Z} noise, to compute the probability for.
#'    \item \code{'X'} for \eqn{x\alpha + dU + W},
#'    \item \code{'W'} for \eqn{W},
#'    \item \code{'Y'} for \eqn{x\alpha + dU + W + Z},
#'    }
#'    }
#'  \item \code{"crps"} A logical variable for calculating
#'    continuous ranked probability score (CRPS); \code{"TRUE"} indicates
#'    calculate, \code{"FALSE"} do not calculate.
#'  \item \code{"nSim"} A numerical value for the number of samples for the Gibbs sampler
#'    to obtain the predictions.
#'  \item \code{"nBurnin"} A numeric value for the number of samples that are
#'    discarded as burn-in whilst calculating the predictions.
#'  \item \code{silent} A logical value for printing the details;
#'      \code{"TRUE"} indicates do not print, \code{"FALSE"} indicates print.
#'  }
#' @return A list of output.
#'
#' @details This function is a wrapper function that calls
#'    \code{"predictLong"} internally.
#'
#' @seealso \code{\link{ngme}}
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   predict(fit, ...)
#'   }


predict.ngme.spatial <- function(object,
                                 id = NULL,
                                 type = "Smoothing",
                                 data = NULL,
                                 quantiles = c(0.025, 0.975),
                                 controls = list(
                                   return.samples = FALSE,
                                   excursions = NULL,
                                   crps = TRUE,
                                   nSim = 1000,
                                   nBurnin = 100,
                                   silent = TRUE,
                                   seed = NULL
                                 )
)
{
  
  
  if(length(controls) < 6){
    controls_full <- list(
      return.samples = FALSE,
      excursions = NULL,
      crps = TRUE,
      nSim = 1000,
      nBurnin = 100,
      silent = TRUE,
      seed = ceiling(10^8 * runif(1))
    )
    for(i in 1:length(controls)){
      controls_full[names(controls)[i]] <- controls[i]
    }
    controls <- controls_full
  }
  
  id_list <- as.numeric(names(object$Y))
  if(length(id_list)==0){
    id_list = 1
  }
  if(length(id)==0){
    id <- id_list
  }
  pInd <- which(id_list %in% id)
  
  if(type=="LOOCV"){
    preds <- predictLong(
      Y                    = object$Y,
      locs                 = object$locs,
      pInd                 = pInd,
      return.samples       = controls$return.samples,
      type                 = type,
      quantiles            = quantiles,
      excursions           = controls$excursions,
      crps                 = controls$crps,
      crps.skip            = 1,
      mixedEffect_list     = object$mixedEffect_list,
      measurment_list      = object$measurementError_list,
      processes_list       = object$processes_list,
      operator_list        = object$operator_list,
      nSim                 = controls$nSim,
      nBurnin              = controls$nBurnin,
      silent               = controls$silent,
      seed                 = controls$seed)
    
    if(controls$silent == FALSE){
      cat("Calculating accuracy measures", "\n")
    }
    
    Y_for_pred <- lapply(1:length(pInd), function(i) object$Y[[pInd[i]]])
    
    if(dim(object$Y[[1]])[2] == 2){
      n.obs <- lapply(1:length(pInd), function(i) dim(object$Y[[pInd[i]]])[1])
      Y1_for_pred <- lapply(1:length(pInd), function(i) object$Y[[pInd[i]]][,1])
      Y2_for_pred <- lapply(1:length(pInd), function(i) object$Y[[pInd[i]]][,2])
      pred_data1 <- data.frame(id       = rep(id, unlist(lapply(Y1_for_pred, length))),
                              observed = unlist(Y1_for_pred),
                              mean     = unlist(lapply(1:length(id), function(i) preds$X.summary[[i]]$Mean[1:n.obs[[i]]])),
                              median   = unlist(lapply(1:length(id), function(i) preds$X.summary[[i]]$Median[1:n.obs[[i]]])),
                              lower    = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$quantiles[[1]]$field[1:n.obs[[i]]])),
                              upper    = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$quantiles[[2]]$field[1:n.obs[[i]]])),
                              crps     = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$crps[1:n.obs[[i]]])))
      ind1 <- !is.na(pred_data1$observed)
      pred_data1$id       = pred_data1$id[ind1]
      pred_data1$observed = pred_data1$observed[ind1]
      pred_data1$mean     = pred_data1$mean[ind1]
      pred_data1$median   = pred_data1$median[ind1]
      pred_data1$lower    = pred_data1$lower[ind1]
      pred_data1$upper    = pred_data1$upper[ind1]
      pred_data1$crps     = pred_data1$crps[ind1]
                               
      pred_data2 <- data.frame(id       = rep(id, unlist(lapply(Y2_for_pred, length))),
                               observed = unlist(Y2_for_pred),
                               mean     = unlist(lapply(1:length(id), function(i) preds$X.summary[[i]]$Mean[(n.obs[[i]]+1):(2*n.obs[[i]])])),
                               median   = unlist(lapply(1:length(id), function(i) preds$X.summary[[i]]$Median[(n.obs[[i]]+1):(2*n.obs[[i]])])),
                               lower    = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$quantiles[[1]]$field[(n.obs[[i]]+1):(2*n.obs[[i]])])),
                               upper    = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$quantiles[[2]]$field[(n.obs[[i]]+1):(2*n.obs[[i]])])),
                               crps     = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$crps[(n.obs[[i]]+1):(2*n.obs[[i]])])))
      ind2 <- !is.na(pred_data2$observed)
      pred_data2$id       = pred_data2$id[ind1]
      pred_data2$observed = pred_data2$observed[ind1]
      pred_data2$mean     = pred_data2$mean[ind1]
      pred_data2$median   = pred_data2$median[ind1]
      pred_data2$lower    = pred_data2$lower[ind1]
      pred_data2$upper    = pred_data2$upper[ind1]
      pred_data2$crps     = pred_data2$crps[ind1]
      
      abs_diff_mean1   <- with(pred_data1, abs(observed - mean))
      abs_diff_median1 <- with(pred_data1, abs(observed - median))
      sq_diff_mean1   <- with(pred_data1, (observed - mean)^2)
      sq_diff_median1 <- with(pred_data1, (observed - median)^2)
      covered1   <- with(pred_data1, lower < observed & upper > observed)
      int.width1 <- with(pred_data1, upper - lower)
      n_obs1 <- nrow(pred_data1)
      
      abs_diff_mean2   <- with(pred_data2, abs(observed - mean))
      abs_diff_median2 <- with(pred_data2, abs(observed - median))
      sq_diff_mean2   <- with(pred_data2, (observed - mean)^2)
      sq_diff_median2 <- with(pred_data2, (observed - median)^2)
      covered2   <- with(pred_data2, lower < observed & upper > observed)
      int.width2 <- with(pred_data2, upper - lower)
      n_obs2 <- nrow(pred_data2)
      mean.rmse.mean.predictor <- c(sqrt(mean(sq_diff_mean1)),sqrt(mean(sq_diff_mean2)))
      mean.rmse.median.predictor <- c(sqrt(mean(sq_diff_median1)),sqrt(mean(sq_diff_median2)))
      
      out <- list(predictions = preds,
                  id = id,
                  type = type,
                  call = match.call(),
                  pred.data1 = pred_data1,
                  pred.data2 = pred_data2,
                  mean.mae.mean.predictor = c(mean(abs_diff_mean1),mean(abs_diff_mean2)),
                  mean.mae.median.predictor = c(mean(abs_diff_median1),mean(abs_diff_median2)),
                  median.mae.mean.predictor = c(median(abs_diff_mean1),median(abs_diff_mean2)),
                  median.mae.median.predictor = c(median(abs_diff_median1),median(abs_diff_median2)),
                  std.mae.mean.predictor = c(sqrt(var(abs_diff_mean1)/n_obs1),sqrt(var(abs_diff_mean2)/n_obs2)),
                  std.mae.median.predictor = c(sqrt(var(abs_diff_median1)/n_obs1),sqrt(var(abs_diff_median2)/n_obs2)),
                  mean.rmse.mean.predictor = mean.rmse.mean.predictor,
                  mean.rmse.median.predictor = mean.rmse.median.predictor,
                  std.rmse.mean.predictor = c(sqrt( (0.5 * (1/mean.rmse.mean.predictor[1]) * var(sq_diff_mean1)) / n_obs1),
                                              sqrt( (0.5 * (1/mean.rmse.mean.predictor[2]) * var(sq_diff_mean2)) / n_obs2)),
                  std.rmse.median.predictor = c(sqrt( (0.5 * (1/mean.rmse.median.predictor[1]) * var(sq_diff_median1)) / n_obs1),
                                                sqrt( (0.5 * (1/mean.rmse.median.predictor[2]) * var(sq_diff_median2)) / n_obs2)),
                  coverage.mean = 100*c(mean(covered1),mean(covered2)),
                  coverage.std = 100*c(sqrt(var(covered1)/n_obs1),sqrt(var(covered2)/n_obs2)),
                  mean.crps = c(mean(pred_data1$crps),mean(pred_data2$crps)),
                  median.crps = c(median(pred_data1$crps),median(pred_data2$crps)),
                  std.crps = c(sqrt(var(pred_data1$crps)/n_obs1),sqrt(var(pred_data2$crps)/n_obs2)),
                  mean.int.width = c(mean(int.width1),mean(int.width2)),
                  std.int.width = c(sqrt(var(int.width1)/n_obs1),sqrt(var(int.width2)/n_obs2)),
                  Y = object$Y,
                  locs = object$locs,
                  id_list = id_list
      )  
      
    } else {
      pred_data <- data.frame(id       = rep(id, unlist(lapply(Y_for_pred, length))),
                              observed = unlist(Y_for_pred),
                              mean     = unlist(lapply(1:length(id), function(i) preds$X.summary[[i]]$Mean)),
                              median   = unlist(lapply(1:length(id), function(i) preds$X.summary[[i]]$Median)),
                              lower    = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$quantiles[[1]]$field)),
                              upper    = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$quantiles[[2]]$field)),
                              crps     = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$crps))
      )
      
      abs_diff_mean   <- with(pred_data, abs(observed - mean))
      abs_diff_median <- with(pred_data, abs(observed - median))
      sq_diff_mean   <- with(pred_data, (observed - mean)^2)
      sq_diff_median <- with(pred_data, (observed - median)^2)
      covered   <- with(pred_data, lower < observed & upper > observed)
      int.width <- with(pred_data, upper - lower)
      n_obs <- nrow(pred_data)
      mean.rmse.mean.predictor <- sqrt(mean(sq_diff_mean))
      mean.rmse.median.predictor <<- sqrt(mean(sq_diff_median))
      
      out <- list(predictions = preds,
                  id = id,
                  type = type,
                  call = match.call(),
                  pred.data = pred_data,
                  mean.mae.mean.predictor = mean(abs_diff_mean),
                  mean.mae.median.predictor = mean(abs_diff_median),
                  median.mae.mean.predictor = median(abs_diff_mean),
                  median.mae.median.predictor = median(abs_diff_median),
                  std.mae.mean.predictor = sqrt(var(abs_diff_mean)/n_obs),
                  std.mae.median.predictor = sqrt(var(abs_diff_median)/n_obs),
                  mean.rmse.mean.predictor = mean.rmse.mean.predictor,
                  mean.rmse.median.predictor = mean.rmse.median.predictor,
                  std.rmse.mean.predictor = sqrt( (0.5 * (1/mean.rmse.mean.predictor) * var(sq_diff_mean)) / n_obs),
                  std.rmse.median.predictor = sqrt( (0.5 * (1/mean.rmse.median.predictor) * var(sq_diff_median)) / n_obs),
                  coverage.mean = 100 * mean(covered),
                  coverage.std = 100 * sqrt(var(covered)/n_obs),
                  mean.crps = mean(pred_data$crps),
                  median.crps = median(pred_data$crps),
                  std.crps = sqrt(var(pred_data$crps)/n_obs),
                  mean.int.width = mean(int.width),
                  std.int.width = sqrt(var(int.width)/n_obs),
                  Y = object$Y,
                  locs = object$locs,
                  id_list = id_list
      )  
    }
    
    
  }  else {
    #kriging prediction
    fixed <-object$call$fixed
    random <- object$call$random
    group.id <- object$call$group.id
    location.names <- as.character(object$call$location.names)[2:3]
    if(is.null(random)){
      use.random = FALSE
    } else {
      use.random = TRUE
    }
    if(use.random){
      idname <- rev(unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE)))[1]
      id <- data[, idname]  
    } else if(!is.null(group.id)){
      idname = group.id
      id <- data[, idname]  
    } else {
      idname = NULL
    }
    response.name <- all.vars(fixed)[1]
    data[response.name] = 0
    effects <- extract.effects(data = data, fixed = fixed,random=random, idname = idname)
    B_fixed = effects$B_fixed
    B_random = effects$B_random
    
    if(is.null(idname)){
      locs <- list(as.matrix(data[, location.names]))  
    } else {
      locs <- tapply(as.matrix(data[, location.names]), id, function(x) x)  
    }
    
    if(dim(object$Y[[1]])[2] == 2){
      fixed2 <-object$call$fixed2
      random2 <- object$call$random2
      response.name2 <- all.vars(fixed2)[1]
      data[response.name2] = 0
      effects2 <- extract.effects(data = data, fixed = fixed2,random=random2, idname = idname)
      Be.pred <- list()
      for(i in 1:length(B_fixed)){
        B_fixed[[i]] <- as.matrix(bdiag(B_fixed[[i]],effects2$B_fixed[[i]]))
        if(use.random){
          B_random[[i]] <- as.matrix(bdiag(B_random[[i]],effects2$B_random[[i]]))
        }
        Be.pred[[i]] <- kronecker(diag(2),matrix(rep(1, dim(locs[[i]])[1])))
      }
      object$measurementError_list$Bpred <- Be.pred
    }
    
    
    
    if(use.random){
      preds <- predictLong(
        Y                    = object$Y,
        locs                 = object$locs,
        locs.pred            = locs,
        pInd                 = pInd,
        Brandom.pred         = B_random,
        Bfixed.pred          = B_fixed,
        return.samples       = controls$return.samples,
        type                 = type,
        quantiles            = quantiles,
        excursions           = controls$excursions,
        mixedEffect_list     = object$mixedEffect_list,
        measurment_list      = object$measurementError_list,
        processes_list       = object$processes_list,
        operator_list        = object$operator_list,
        nSim                 = controls$nSim,
        nBurnin              = controls$nBurnin,
        silent               = controls$silent,
        seed                 = controls$seed)  
    } else {
      preds <- predictLong(
        Y                    = object$Y,
        locs                 = object$locs,
        locs.pred            = locs,
        pInd                 = pInd,
        Bfixed.pred          = B_fixed,
        return.samples       = controls$return.samples,
        type                 = type,
        quantiles            = quantiles,
        excursions           = controls$excursions,
        mixedEffect_list     = object$mixedEffect_list,
        measurment_list      = object$measurementError_list,
        processes_list       = object$processes_list,
        operator_list        = object$operator_list,
        nSim                 = controls$nSim,
        nBurnin              = controls$nBurnin,
        silent               = controls$silent,
        seed                 = controls$seed)
    }
    
    out <- list(predictions = preds,
                id = id,
                type = type,
                call = match.call(),
                Y = object$Y,
                locs = object$locs,
                id_list = id_list
    )
  }
  
  
  if(dim(object$Y[[1]])[2] == 2){
    preds <- out$predictions
    for(i in 1:length(object$Y)){
      n.p = length(preds$X.summary[[i]]$Mean)/2
      preds$Y.summary[[i]]$Mean <- matrix(preds$Y.summary[[i]]$Mean, n.p,2)
      preds$Y.summary[[i]]$Var <- matrix(preds$Y.summary[[i]]$Var, n.p,2)
      preds$Y.summary[[i]]$Median <- matrix(preds$Y.summary[[i]]$Median, n.p,2)
      for(j in 1:length(preds$Y.summary[[i]]$quantiles)){
        preds$Y.summary[[i]]$quantiles[[j]]$field <- matrix(preds$Y.summary[[i]]$quantiles[[j]]$field, n.p,2)
      }
      
      preds$X.summary[[i]]$Mean <- matrix(preds$X.summary[[i]]$Mean, n.p,2)
      preds$X.summary[[i]]$Var <- matrix(preds$X.summary[[i]]$Var, n.p,2)
      preds$X.summary[[i]]$Median <- matrix(preds$X.summary[[i]]$Median, n.p,2)
      for(j in 1:length(preds$X.summary[[i]]$quantiles)){
        preds$X.summary[[i]]$quantiles[[j]]$field <- matrix(preds$X.summary[[i]]$quantiles[[j]]$field, n.p,2)
      }
      
      preds$W.summary[[i]]$Mean <- matrix(preds$W.summary[[i]]$Mean, n.p,2)
      preds$W.summary[[i]]$Var <- matrix(preds$W.summary[[i]]$Var, n.p,2)
      preds$W.summary[[i]]$Median <- matrix(preds$W.summary[[i]]$Median, n.p,2)
      for(j in 1:length(preds$W.summary[[i]]$quantiles)){
        preds$W.summary[[i]]$quantiles[[j]]$field <- matrix(preds$W.summary[[i]]$quantiles[[j]]$field, n.p,2)
      }
      
      preds$V.summary[[i]]$Mean <- matrix(preds$V.summary[[i]]$Mean, n.p,2)
      preds$V.summary[[i]]$Var <- matrix(preds$V.summary[[i]]$Var, n.p,2)
      preds$V.summary[[i]]$Median <- matrix(preds$V.summary[[i]]$Median, n.p,2)
      for(j in 1:length(preds$V.summary[[i]]$quantiles)){
        preds$V.summary[[i]]$quantiles[[j]]$field <- matrix(preds$V.summary[[i]]$quantiles[[j]]$field, n.p,2)
      }
    }
    
    out$predictions <- preds
  }
  
  
  class(out) <- "predict.ngme.spatial"
  out
  
}

