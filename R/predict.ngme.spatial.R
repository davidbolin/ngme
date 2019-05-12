
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
                                   silent = TRUE
                                 )
)
{
  
  controls$seed <- ceiling(10^8 * runif(1))
  
  if(length(controls) < 11){
    controls_full <- list(
      return.samples = FALSE,
      excursions = NULL,
      crps = TRUE,
      crps.skip = 1,
      nSim = 1000,
      nBurnin = 100,
      silent = TRUE,
      n.cores = 1,
      batch.size = 100
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
    
    pred_data <- data.frame(id       = rep(id, unlist(lapply(Y_for_pred, length))),
                            time     = unlist(lapply(1:length(pInd), function(i) object$locs[[pInd[i]]])),
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
    
    mean.mae.mean.predictor       <- mean(abs_diff_mean)
    mean.mae.median.predictor     <- mean(abs_diff_median)
    median.mae.mean.predictor       <- median(abs_diff_mean)
    median.mae.median.predictor     <- median(abs_diff_median)
    std.mae.mean.predictor        <- sqrt(var(abs_diff_mean)/n_obs)
    std.mae.median.predictor      <- sqrt(var(abs_diff_median)/n_obs)
    
    mean.rmse.mean.predictor      <- sqrt(mean(sq_diff_mean))
    mean.rmse.median.predictor    <- sqrt(mean(sq_diff_median))
    std.rmse.mean.predictor       <- sqrt( (0.5 * (1/mean.rmse.mean.predictor) * var(sq_diff_mean)) / n_obs)
    std.rmse.median.predictor     <- sqrt( (0.5 * (1/mean.rmse.median.predictor) * var(sq_diff_median)) / n_obs)
    
    coverage.mean  <- 100 * mean(covered)
    coverage.std   <- 100 * sqrt(var(covered)/n_obs)
    
    mean.crps      <- mean(pred_data$crps)
    median.crps      <- median(pred_data$crps)
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
                median.mae.mean.predictor = median.mae.mean.predictor,
                median.mae.median.predictor = median.mae.median.predictor,
                std.mae.mean.predictor = std.mae.mean.predictor,
                std.mae.median.predictor = std.mae.median.predictor,
                mean.rmse.mean.predictor = mean.rmse.mean.predictor,
                mean.rmse.median.predictor = mean.rmse.median.predictor,
                std.rmse.mean.predictor = std.rmse.mean.predictor,
                std.rmse.median.predictor = std.rmse.median.predictor,
                coverage.mean = coverage.mean,
                coverage.std = coverage.std,
                mean.crps = mean.crps,
                median.crps = median.crps,
                std.crps = std.crps,
                mean.int.width = mean.int.width,
                std.int.width = std.int.width,
                Y = object$Y,
                locs = object$locs,
                id_list = id_list
    )
    
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
    mf_fixed <- model.frame(formula = fixed, data = data)
    x_fixed_f  <- as.matrix(model.matrix(attr(mf_fixed, "terms"), data = mf_fixed))
    colnames(x_fixed_f)[1] <- gsub("[[:punct:]]", "", colnames(x_fixed_f)[1])
    
    # excluding the intercept and the covariates that are specified in random
    if(use.random){
      cov_list_fixed  <- attr(terms(fixed), "term.labels")
      cov_list_random <- unlist(strsplit(attr(terms(random), "term.labels"), " | ", fixed = TRUE))
      cov_list_random <- c(strsplit(cov_list_random[1], " + ", fixed=TRUE)[[1]], cov_list_random[2])
      cov_list_random <- cov_list_random[-length(cov_list_random)]  
      to_del_x_fixed <- c("Intercept", cov_list_fixed[(cov_list_fixed %in% cov_list_random)])
      x_fixed <- x_fixed_f[, !(colnames(x_fixed_f) %in% to_del_x_fixed), drop = FALSE]
      #random effects design matrix
      random_names             <- unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE))
      random_names_id_excluded <- random_names[!(random_names %in% idname)]
      random_formula           <- as.formula(paste("~", paste(random_names_id_excluded, collapse = "+")))
      
      mf_random <- model.frame(formula = random_formula, data = data)
      x_random  <- as.matrix(model.matrix(attr(mf_random, "terms"), data = mf_random))
      colnames(x_random)[1] <- gsub("[[:punct:]]", "", colnames(x_random)[1])
      
      idlist <- unique(id)
      data_random <- data.frame(cbind(id, x_random))
      B_random    <- split(data_random[, -1], data_random[,1])
      B_random    <- lapply(B_random, function(x) as.matrix(x))
      
      data_fixed <- data.frame(cbind(id, x_fixed))
      B_fixed    <- split(data_fixed[, -1], data_fixed[,1])
      B_fixed    <- lapply(B_fixed, function(x) as.matrix(x))  
    } else {
      x_fixed <- x_fixed_f
      x_random <- NA
      if(is.null(idname)){ # no replicates
        B_fixed = list(x_fixed)
      } else {
        data_fixed <- data.frame(cbind(id, x_fixed))
        B_fixed    <- split(data_fixed[, -1], data_fixed[,1])
        B_fixed    <- lapply(B_fixed, function(x) as.matrix(x))
      }
    }
    
    if(is.null(idname)){
      locs <- list(as.matrix(data[, location.names]))  
    } else {
      locs <- tapply(as.matrix(data[, location.names]), id, function(x) x)  
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
  
  
  
  
  
  class(out) <- "predict.ngme.spatial"
  out
  
}

