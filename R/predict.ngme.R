
#' @title Prediction.
#'
#' @description Obtains predicted values based on filtering and smoothing distributions.
#'
#' @param object A fitted object obtained by calling \code{"ngme"}.
#' @param id A numeric vector containing the ID's of the subjects for whom
#'   predictions are to be obtained. Default is set to \code{"NULL"}
#'   indicating perform predictions for all the subjects.
#' @param type A character string for the type of prediction: \code{"Filter"} for
#'   filtering, \code{"Smoothing"} for smoothing.
#' @param quantiles A two-elemnent vector that contains the quantiles
#'   of the predictions to be calculated.
#' @param controls A list of control variables.
#'  \itemize{
#'  \item \code{"return.samples"} A logical variable for returning the
#'    Monte Carlo samples used to compute the predictions; \code{"TRUE"} indicates
#'    return, \code{"FALSE"} do not return.
#'  \item \code{"predict.derivatives"} A logical variable (????) to obtain
#'    predictions for the derivatives of the process
#'  \item \code{"excursions"} A list of excursion probabilities to compute.
#'    Each list should contain:
#'    \itemize{
#'    \item \code{"type"} - type of excursion '>' or '<',
#'    \item \code{"level"} - level to compute excursion probability for,
#'    \item \code{"process"} - which expression for the model,
#'    \eqn{x\alpha + dU + W + Z} with \eqn{x \alpha} being fixed effects,
#'    \eqn{dU} random effects and \eqn{Z} noise, to compute the probability for.
#'    \code{'X'} for \eqn{x\alpha + dU + W},
#'    \code{'W'} for \eqn{W},
#'    \code{'Y'} for \eqn{x\alpha + dU + W + Z},
#'    \code{'Xderivative'} for the first derivarive of \eqn{x\alpha + dU + W},
#'    and
#'    \code{'Wderivative'} for the first derivariate of \eqn{W}.
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
#'  \item \code{n.cores} Number of cores to use for parallell computations (default =1);
#'  \item \code{batch.size} Number of patients to include in each parallell batch;
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
                            return.preds = TRUE,
                            n.cores = 1,
                            batch.size = 10
                            )
                          )
  {

  controls$seed <- ceiling(10^8 * runif(1))

  if(length(controls) < 12){
    controls_full <- list(
      return.samples = FALSE,
      predict.derivaties = NULL,
      excursions = NULL,
      crps = TRUE,
      crps.skip = 1,
      nSim = 1000,
      nBurnin = 100,
      silent = TRUE,
      return.preds = TRUE,
      n.cores = 1,
      batch.size = 10
      )
    for(i in 1:length(controls)){
        controls_full[names(controls)[i]] <- controls[i]
    }
    controls <- controls_full
  }


  id_list <- as.numeric(names(object$Y))
  if(is.null(id)){
    id <- id_list
  }
  pInd <- which(id_list %in% id)


  #do the prediction in batches of 100 patients
  pInd.list <- list()
  pInd.tmp <- pInd

  if(length(pInd.tmp)>=controls$batch.size){
    k = 1
    while(length(pInd.tmp)>=controls$batch.size){
      if(length(pInd.tmp)>=controls$batch.size){
        pInd.list[[k]] <- pInd.tmp[1:controls$batch.size]
        if(controls$batch.size<length(pInd.tmp)){
          pInd.tmp <- pInd.tmp[(controls$batch.size+1):length(pInd.tmp)]
        } else {
          pInd.tmp <- NULL
        }
        k = k+1
      } else {
        pInd.list[[k]] <- pInd.tmp
        pInd.tmp <- NULL
      }
    }
  } else {
    pInd.list[[1]] <- pInd
  }

  iterations <- length(pInd.list)
  if(controls$n.cores==1){
    preds.list <- list()
    for(i in 1:iterations)
    {

      #cat(object.size(preds.list,units = "MB",standard = "SI"),"\n")
      if(object$use_process == TRUE){
        preds.list[[i]] <- predictLong( Y                    = object$Y,
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
                           silent               = TRUE,
                           seed                 = controls$seed)

      }else{
        preds.list[[i]] <- predictLong(Y                     = object$Y,
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
                          silent                 = TRUE,
                          seed                   = controls$seed)
      }
    }
  } else {
    cl <- makeCluster(controls$n.cores)
    registerDoSNOW(cl)

    pb <- txtProgressBar(max = length(pInd.list), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    clusterExport(cl, list = c('object','pInd.list', 'controls','type','quantiles'),envir=environment())

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
                           silent               = TRUE,
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
                          silent                 = TRUE,
                          seed                   = controls$seed)
      }
      return(pl)
    }
    close(pb)
    stopCluster(cl)
  }


  preds <- merge.pred.lists(preds.list)
  mae <- covered <- int.width <- crps <- rmse <- NULL

  n.obs <- rep(0, length(pInd))

  if(controls$silent == FALSE){
    cat("Calculating accuracy measures", "\n")
  }

  for(i in 1:length(pInd)){
    mae <- c(mae, abs(preds$X.summary[[i]]$Median - object$Y[[pInd[i]]]))
    n.obs[i] <- length(object$Y[[pInd[i]]])
    covered <- c(covered,(preds$Y.summary[[i]]$quantiles[[1]]$field < object$Y[[pInd[i]]]) & (preds$Y.summary[[i]]$quantiles[[2]]$field > object$Y[[pInd[i]]]))
    int.width <- c(int.width, preds$Y.summary[[i]]$quantiles[[2]]$field - preds$Y.summary[[i]]$quantiles[[1]]$field)
    crps <- c(crps, preds$Y.summary[[i]]$crps)
    rmse <- c(rmse, (preds$X.summary[[i]]$Mean - Y[[pInd[i]]])^2)
  }

  sum.n.obs <- sum(n.obs)

  mean.mae       <- mean(mae)
  std.mae        <- sqrt(var(mae)/sum.n.obs)
  coverage.mean  <- 100 * mean(covered)
  coverage.std   <- 100 * sqrt(var(covered)/sum.n.obs)
  mean.rmse      <- sqrt(mean(rmse))
  std.rmse       <- sqrt(var(rmse)/sum.n.obs)
  mean.crps      <- mean(crps)
  std.crps       <- sqrt(var(crps)/sum.n.obs)
  mean.int.width <- mean(int.width)
  std.int.width  <- sqrt(var(int.width)/sum.n.obs)

  if(controls$return.preds == TRUE){
    out <- list(predictions = preds,
                id = id,
                type = type,
                call = match.call(),
                mae = mae,
                n.obs = n.obs,
                covered = covered,
                int.width = int.width,
                crps = crps,
                rmse = rmse,
                mean.mae = mean.mae,
                std.mae = std.mae,
                coverage.mean = coverage.mean,
                coverage.std = coverage.std,
                mean.rmse = mean.rmse,
                std.rmse = std.rmse,
                mean.crps = mean.crps,
                std.crps = std.crps,
                mean.int.width = mean.int.width,
                std.int.width = std.int.width,
                Y = object$Y,
                locs = object$locs,
                id_list = id_list)
  }else{
    out <- list(predictions = preds,
                id = id,
                type = type,
                call = match.call(),
                mae = mae,
                n.obs = n.obs,
                covered = covered,
                int.width = int.width,
                crps = crps,
                rmse = rmse,
                mean.mae = mean.mae,
                std.mae = std.mae,
                coverage.mean = coverage.mean,
                coverage.std = coverage.std,
                mean.rmse = mean.rmse,
                std.rmse = std.rmse,
                mean.crps = mean.crps,
                std.crps = std.crps,
                mean.int.width = mean.int.width,
                std.int.width = std.int.width,
                Y = object$Y,
                locs = object$locs,
                id_list = id_list)
  }

  class(out) <- "predict.ngme"
  out

}

merge.pred.lists <- function(preds.list){
  preds <- preds.list[[1]]
  if(length(preds.list)>1){
    for(i in 2:length(preds.list)){
      for(j in 1:length(names(preds))){
        preds[[names(preds)[j]]] <- append(preds[[names(preds)[j]]],preds.list[[i]][[names(preds)[j]]])
      }
    }
  }
  return(preds)
}
