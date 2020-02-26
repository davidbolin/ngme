
#' @title Prediction.
#'
#' @description Obtains predicted values based on filtering and smoothing distributions.
#'
#' @param object A fitted object obtained by calling \code{"ngme"}.
#' @param id A numeric vector containing the ID's of the subjects for whom
#'   predictions are to be obtained. Default is set to \code{"NULL"}
#'   indicating perform predictions for all the subjects.
#' @param type A character string for the type of prediction: \code{"Nowcast"} 
#'   for nowcasting, \code{"Filter"} for
#'   one step-ahead forecasting, \code{"Smoothing"} for smoothing.
#' @param quantiles A two-element vector that contains the quantiles
#'   of the predictions to be calculated.
#' @param controls A list of control variables.
#'  \itemize{
#'  \item \code{"return.samples"} A logical variable for returning the
#'    Monte Carlo samples used to compute the predictions; \code{"TRUE"} indicates
#'    return, \code{"FALSE"} do not return.
#'  \item \code{"predict.derivatives"} A list for calculating excursion probabilities
#'    \itemize{
#'    \item \code{Bfixed} A list of fixed effects covariate matrices at \eqn{t_{ij}} + \code{delta}.
#'    \item \code{Brandom} A list of random effects covariate matrices at \eqn{t_{ij}}+ \code{delta}.
#'    \item \code{delta} A numeric value indicating \eqn{\Delta}t for calculating the derivative numerically.
#'    }
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
#'    \item \code{'Xderivative'} for the first derivarive of \eqn{x\alpha + dU + W},
#'    and
#'    \item \code{'Wderivative'} for the first derivariate of \eqn{W}.
#'    }
#'    }
#'  \item \code{"crps"} A logical variable for calculating
#'    continuous ranked probability score (CRPS); \code{"TRUE"} indicates
#'    calculate, \code{"FALSE"} do not calculate.
#'  \item \code{"crps.skip"} A numerical value, say a, that indicates every \emph{a}th
#'    element of the sample to be used to compute the crps score.
#'  \item \code{"nSim"} A numerical value for the number of samples for the Gibbs sampler
#'    to obtain the predictions.
#'  \item \code{"nBurnin"} A numeric value for the number of samples that are
#'    discarded as burn-in whilst calculating the predictions.
#'  \item \code{silent} A logical value for printing the details;
#'      \code{"TRUE"} indicates do not print, \code{"FALSE"} indicates print.
#'  \item \code{n.cores} Number of cores to use for paralell computations, default is set to 1.
#'  \item \code{batch.size} Number of subjects to include in each paralell batch.
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

predict.ngme <- function(object,
                         newdata = NULL,
                         id = NULL,
                         type = "Filter",
                         quantiles = c(0.025, 0.975),
                         controls = list(
                            return.samples = FALSE,
                            predict.derivatives = NULL,
                            excursions = NULL,
                            crps = TRUE,
                            crps.skip = 1,
                            nSim = 1000,
                            nBurnin = 100,
                            silent = TRUE,
                            n.cores = 1,
                            batch.size = 100
                            )
                          )
  {

  controls$seed <- ceiling(10^8 * runif(1))

  if(length(controls) < 11){
    controls_full <- list(
      return.samples = FALSE,
      predict.derivaties = NULL,
      excursions = NULL,
      crps = TRUE,
      crps.skip = 1,
      nSim = 1000,
      nBurnin = 100,
      silent = TRUE,
      n.cores = 1,
      batch.size = 100)
    
    for(i in 1:length(controls)){
        controls_full[names(controls)[i]] <- controls[i]
    }
    controls <- controls_full
  }

  ## when new data were provided
  if(is.null(newdata) == FALSE){
  
    new_ids    <- unique(newdata[, object$idname])
    n_new_subj <- length(new_ids)
    n_ava_subj <- length(object$Y)
    
    if(is.null(id)) id <- new_ids
    
    new_effects <- extract.effects(data = newdata, 
                                   fixed = object$fixed_formula,
                                   random = object$random_formula, 
                                   idname = object$idname)
    
    object$Y[(n_ava_subj + 1): (n_ava_subj + n_new_subj)] <- new_effects$Y
    names(object$Y)[(n_ava_subj + 1): (n_ava_subj + n_new_subj)] <- new_ids
    
    object$mixedEffect_list$B_fixed[(n_ava_subj + 1): (n_ava_subj + n_new_subj)] <- new_effects$B_fixed
    object$mixedEffect_list$B_random[(n_ava_subj + 1): (n_ava_subj + n_new_subj)] <- new_effects$B_random
    
    if(object$use_process == TRUE){
      
      new_locs <- tapply(as.matrix(newdata[, object$call$timeVar]), newdata[, object$idname], function(x) x)
      
      object$locs[(n_ava_subj + 1): (n_ava_subj + n_new_subj)] <- new_locs
      names(object$locs)[(n_ava_subj + 1): (n_ava_subj + n_new_subj)] <- new_ids
      
      new_operator_list <- create_operator(new_locs,
                                           name = object$operator_type,
                                           common.grid = object$mesh$common.grid,
                                           extend  = object$mesh$extend,
                                           max.dist = object$mesh$max.dist,
                                           cutoff = object$mesh$cutoff,
                                           n.cores = object$mesh$n.cores)
      
      for(i in 1:n_new_subj){
        object$processes_list$X[[n_ava_subj + i]] <- rep(0, length(new_operator_list$h[[i]]))
        object$processes_list$W[[n_ava_subj + i]] <- rep(0, length(new_operator_list$h[[i]]))
        object$processes_list$V[[n_ava_subj + i]] <- new_operator_list$h[[i]]
        
        object$operator_list$Q[[n_ava_subj + i]]   <- new_operator_list$Q[[i]]
        object$operator_list$h[[n_ava_subj + i]]   <- new_operator_list$h[[i]]
        object$operator_list$loc[[n_ava_subj + i]] <- new_operator_list$loc[[i]]
      } 

    }
  }

  id_list <- as.numeric(names(object$Y))
  if(is.null(id)){
    id <- id_list
  }
  pInd <- which(id_list %in% id)

  #do the prediction in batches
  # pInd.list <- list()
  # pInd.tmp <- pInd
  #
  # if(length(pInd.tmp)>=controls$batch.size){
  #   k = 1
  #   while(length(pInd.tmp)>=controls$batch.size){
  #     if(length(pInd.tmp)>=controls$batch.size){
  #       pInd.list[[k]] <- pInd.tmp[1:controls$batch.size]
  #       if(controls$batch.size<length(pInd.tmp)){
  #         pInd.tmp <- pInd.tmp[(controls$batch.size+1):length(pInd.tmp)]
  #       } else {
  #         pInd.tmp <- NULL
  #       }
  #       k = k+1
  #     } else {
  #       pInd.list[[k]] <- pInd.tmp
  #       pInd.tmp <- NULL
  #     }
  #   }
  # } else {
  #   pInd.list[[1]] <- pInd
  # }

  batch.size <- controls$batch.size
  iterations <- ceiling(length(pInd)/batch.size)

  pInd.list <- lapply(1:iterations, function(i) na.omit(pInd[(batch.size*(i - 1) + 1) : (batch.size * i)]))

  # iterations <- length(pInd.list)

  if(controls$n.cores == 1){

    preds.list <- list()

    for(i in 1:iterations){

      if(controls$silent == FALSE){
        cat("\n")
        cat("Iteration", i, "out of", iterations, "\n")
        cat("\n")
      }

      #cat(object.size(preds.list,units = "MB",standard = "SI"),"\n")
      if(object$use_process == TRUE){
        preds.list[[i]] <- predictLong(
                           Y                    = object$Y,
                           locs                 = object$locs,
                           pInd                 = pInd.list[[i]],
                           locs.pred            = object$locs,
                           Brandom.pred         = object$mixedEffect_list$B_random,
                           Bfixed.pred          = object$mixedEffect_list$B_fixed,
                           return.samples       = controls$return.samples,
                           type                 = type,
                           quantiles            = quantiles,
                           predict.derivatives  = controls$predict.derivatives,
                           excursions           = controls$excursions,
                           crps                 = controls$crps,
                           crps.skip            = controls$crps.skip,
                           mixedEffect_list     = object$mixedEffect_list,
                           measurment_list      = object$measurementError_list,
                           processes_list       = object$processes_list,
                           operator_list        = object$operator_list,
                           nSim                 = controls$nSim,
                           nBurnin              = controls$nBurnin,
                           silent               = controls$silent,
                           seed                 = controls$seed)

      }else{
        preds.list[[i]] <- predictLong(
                          Y                      = object$Y,
                          locs                   = object$locs,
                          pInd                   = pInd.list[[i]],
                          locs.pred              = object$locs,
                          Brandom.pred           = object$mixedEffect_list$B_random,
                          Bfixed.pred            = object$mixedEffect_list$B_fixed,
                          return.samples         = controls$return.samples,
                          type                   = type,
                          quantiles              = quantiles,
                          predict.derivatives    = controls$predict.derivatives,
                          excursions             = controls$excursions,
                          crps                   = controls$crps,
                          crps.skip              = controls$crps.skip,
                          mixedEffect_list       = object$mixedEffect_list,
                          measurment_list        = object$measurementError_list,
                          operator_list          = object$operator_list,
                          nSim                   = controls$nSim,
                          nBurnin                = controls$nBurnin,
                          silent                 = controls$silent,
                          seed                   = controls$seed)
      }
    }
  } else {

    cl <- makeCluster(controls$n.cores)

    registerDoSNOW(cl)

    pb <- txtProgressBar(max = length(pInd.list), style = 3)

    progress <- function(n) setTxtProgressBar(pb, n)

    opts <- list(progress = progress)

    parallel::clusterExport(cl, varlist = c('object','pInd.list', 'controls','type','quantiles'), envir = environment())

    preds.list <- foreach(i = 1:iterations, .options.snow = opts) %dopar%
    {
      if(object$use_process == TRUE){
        pl <- predictLong( Y                    = object$Y,
                           locs                 = object$locs,
                           pInd                 = pInd.list[[i]],
                           locs.pred            = object$locs,
                           Brandom.pred         = object$mixedEffect_list$B_random,
                           Bfixed.pred          = object$mixedEffect_list$B_fixed,
                           return.samples       = controls$return.samples,
                           type                 = type,
                           quantiles            = quantiles,
                           predict.derivatives  = controls$predict.derivatives,
                           excursions           = controls$excursions,
                           crps                 = controls$crps,
                           crps.skip            = controls$crps.skip,
                           mixedEffect_list     = object$mixedEffect_list,
                           measurment_list      = object$measurementError_list,
                           processes_list       = object$processes_list,
                           operator_list        = object$operator_list,
                           nSim                 = controls$nSim,
                           nBurnin              = controls$nBurnin,
                           silent               = controls$silent,
                           seed                 = controls$seed)
      }else{
        pl <- predictLong(Y                     = object$Y,
                          locs                   = object$locs,
                          pInd                   = pInd.list[[i]],
                          locs.pred              = object$locs,
                          Brandom.pred           = object$mixedEffect_list$B_random,
                          Bfixed.pred            = object$mixedEffect_list$B_fixed,
                          return.samples         = controls$return.samples,
                          type                   = type,
                          quantiles              = quantiles,
                          predict.derivatives    = controls$predict.derivatives,
                          excursions             = controls$excursions,
                          crps                   = controls$crps,
                          crps.skip              = controls$crps.skip,
                          mixedEffect_list       = object$mixedEffect_list,
                          measurment_list        = object$measurementError_list,
                          operator_list          = object$operator_list,
                          nSim                   = controls$nSim,
                          nBurnin                = controls$nBurnin,
                          silent                 = controls$silent,
                          seed                   = controls$seed)
      }
      return(pl)
    }
    close(pb)
    stopCluster(cl)
  }


  preds <- merge.pred.lists(preds.list, pInd)


  if(controls$silent == FALSE){
    cat("Calculating accuracy measures", "\n")
  }

  Y_for_pred <- lapply(1:length(pInd), function(i) object$Y[[pInd[i]]])

  pred_data <- data.frame(id       = rep(id, unlist(lapply(Y_for_pred, length))),
                          time     = unlist(lapply(1:length(pInd), function(i) object$locs[[pInd[i]]])),
                          observed = unlist(Y_for_pred),
                          mean     = unlist(lapply(1:length(id), function(i) preds$X.summary[[i]]$Mean)),
                          median   = unlist(lapply(1:length(id), function(i) preds$X.summary[[i]]$Median)),
                          lower    = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$quantiles[[1]]$field)),
                          upper    = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$quantiles[[2]]$field)),
                          crps     = unlist(lapply(1:length(id), function(i) preds$Y.summary[[i]]$crps))
                          )

  # if(type == "Filter"){
  #   pred_data <- pred_data[duplicated(pred_data$id), ]
  # }

  abs_diff_mean   <- with(pred_data, abs(observed - mean))
  abs_diff_median <- with(pred_data, abs(observed - median))

  sq_diff_mean   <- with(pred_data, (observed - mean)^2)
  sq_diff_median <- with(pred_data, (observed - median)^2)

  covered   <- with(pred_data, lower < observed & upper > observed)
  int.width <- with(pred_data, upper - lower)

  n_obs <- nrow(pred_data)

  mean.mae.mean.predictor       <- mean(abs_diff_mean)
  mean.mae.median.predictor     <- mean(abs_diff_median)
  std.mae.mean.predictor        <- sqrt(var(abs_diff_mean)/n_obs)
  std.mae.median.predictor      <- sqrt(var(abs_diff_median)/n_obs)

  mean.rmse.mean.predictor      <- sqrt(mean(sq_diff_mean))
  mean.rmse.median.predictor    <- sqrt(mean(sq_diff_median))
  std.rmse.mean.predictor       <- sqrt( (0.5 * (1/mean.rmse.mean.predictor) * var(sq_diff_mean)) / n_obs)
  std.rmse.median.predictor     <- sqrt( (0.5 * (1/mean.rmse.median.predictor) * var(sq_diff_median)) / n_obs)

  coverage.mean  <- 100 * mean(covered)
  coverage.std   <- 100 * sqrt(var(covered)/n_obs)

  mean.crps      <- mean(pred_data$crps)
  std.crps       <- sqrt(var(pred_data$crps)/n_obs)

  mean.int.width <- mean(int.width)
  std.int.width  <- sqrt(var(int.width)/n_obs)

    out <- list(predictions = preds,
                id = id,
                type = type,
                call = match.call(),
                pred.data = pred_data,
                mean.mae.mean.predictor = mean.mae.mean.predictor,
                mean.mae.median.predictor = mean.mae.median.predictor,
                std.mae.mean.predictor = std.mae.mean.predictor,
                std.mae.median.predictor = std.mae.median.predictor,
                mean.rmse.mean.predictor = mean.rmse.mean.predictor,
                mean.rmse.median.predictor = mean.rmse.median.predictor,
                std.rmse.mean.predictor = std.rmse.mean.predictor,
                std.rmse.median.predictor = std.rmse.median.predictor,
                coverage.mean = coverage.mean,
                coverage.std = coverage.std,
                mean.crps = mean.crps,
                std.crps = std.crps,
                mean.int.width = mean.int.width,
                std.int.width = std.int.width,
                Y = object$Y,
                locs = object$locs,
                id_list = id_list
                )

  class(out) <- "predict.ngme"
  out

}

merge.pred.lists <- function(preds.list, pInd){
  #merge lists
  preds <- preds.list[[1]]
  if(length(preds.list)>1){
    for(i in 2:length(preds.list)){
      for(j in 1:length(names(preds))){
        preds[[names(preds)[j]]] <- append(preds[[names(preds)[j]]],preds.list[[i]][[names(preds)[j]]])
      }
    }
  }
  #sort in correct order based on pInd
  ind = match(pInd,preds$pInd)
  for(j in 1:length(names(preds))){
    preds[[names(preds)[j]]] <- preds[[names(preds)[j]]][ind]
  }
  return(preds)
}




