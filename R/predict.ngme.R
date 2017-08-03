
#' @title Prediction.
#'
#' @description Obtains predicted values based on filtering and smoothing distributions.
#'
#' @param object A fitted object obtained by calling \code{"ngme"}.
#' @param id A numeric vector containing the ID's of the subjects for whom
#'   predictions are to be obtained.
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
                          id,
                          type = "Filter",
                          quantiles = c(0.025, 0.975),
                          controls = list(
                            return.samples = TRUE,
                            predict.derivaties = NULL,
                            excursions = NULL,
                            crps = TRUE,
                            crps.skip = 1,
                            nSim = 10,
                            nBurnin = 10,
                            silent = TRUE
                            )
                          )
  {

  if(length(controls) < 8){
    controls_full <- list(
      return.samples = TRUE,
      predict.derivaties = NULL,
      excursions = NULL,
      crps = TRUE,
      crps.skip = 1,
      nSim = 10,
      nBurnin = 10,
      silent = TRUE
    )
    for(i in 1:length(controls)){
        controls_full[names(controls)[i]] <- controls[i]
    }
  }
  controls <- controls_full
  
  Y                <- object$Y
  locs             <- object$locs
  mixedEffect_list <- object$mixedEffect_list
  measurementError_list <- object$measurementError_list
  processes_list   <- object$processes_list
  operator_list    <- object$operator_list

  id_list <- as.numeric(names(Y))
  pInd <- which(id_list %in% id)

  B_fixed <- mixedEffect_list$B_fixed
  B_random <- mixedEffect_list$B_random

  if(object$use_process == TRUE){
    predict <- predictLong(Y = Y,
                           locs = locs,
                           pInd = pInd,
                           locs.pred = locs,
                           Brandom.pred = B_random,
                           Bfixed.pred = B_fixed,
                           return.samples = controls$return.samples,
                           type = type,
                           quantiles = quantiles,
                           predict.derivatives = controls$predict.derivatives,
                           excursions = controls$excursions,
                           crps = controls$crps,
                           crps.skip = controls$crps.skip,
                           mixedEffect_list = mixedEffect_list,
                           measurment_list = measurementError_list,
                           processes_list = processes_list,
                           operator_list = operator_list,
                           nSim  = controls$nSim,
                           nBurnin = controls$nBurnin,
                           silent  = controls$silent
                           )
  }else{
    predict <- predictLong(Y = Y,
                           locs = locs,
                           pInd = pInd,
                           locs.pred = locs,
                           Brandom.pred = B_random,
                           Bfixed.pred = B_fixed,
                           return.samples = controls$return.samples,
                           type = type,
                           quantiles = quantiles,
                           predict.derivatives = controls$predict.derivatives,
                           excursions = controls$excursions,
                           crps = controls$crps,
                           crps.skip = controls$crps.skip,
                           mixedEffect_list = mixedEffect_list,
                           measurment_list = measurementError_list,
                           operator_list = operator_list,
                           nSim  = controls$nSim,
                           nBurnin = controls$nBurnin,
                           silent  = controls$silent
                           )
  }
  
  mae <- covered <- int.width <- crps <- rmse <- NULL

  n.obs <- rep(0, length(pInd))

  for(i in 1:length(pInd)){
    mae <- c(mae, abs(predict$X.summary[[i]]$Median - Y[[pInd[i]]]))
    n.obs[i] <- length(Y[[pInd[i]]])
    covered <- c(covered,(predict$Y.summary[[i]]$quantiles[[1]]$field < Y[[pInd[i]]]) & (predict$Y.summary[[i]]$quantiles[[2]]$field > Y[[pInd[i]]]))
    int.width <- c(int.width, predict$Y.summary[[i]]$quantiles[[2]]$field - predict$Y.summary[[i]]$quantiles[[1]]$field)
    crps <- c(crps, predict$Y.summary[[i]]$crps)
    rmse <- c(rmse, (predict$X.summary[[i]]$Mean - Y[[pInd[i]]])^2)
  }

  mean.mae       <- mean(mae)
  std.mae        <- sqrt(var(mae)/sum(n.obs))
  coverage.mean  <- 100 * mean(covered)
  coverage.std   <- 100 * sqrt(var(covered)/sum(n.obs))
  mean.rmse      <- sqrt(mean(rmse))
  std.rmse       <- sqrt(var(rmse)/sum(n.obs))
  mean.crps      <- mean(crps)
  std.crps       <- sqrt(var(crps))
  mean.int.width <- mean(int.width)
  std.int.width  <- sqrt(var(int.width))

  out <- list(predictions = predict,
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
              Y = Y,
              locs = locs,
              id_list = id_list
              )
  class(out) <- "predict.ngme"
  out

}
