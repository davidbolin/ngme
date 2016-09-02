#' @param   Y           - list with the observations
#' @param   locs        - list with positions of the observations (Y)
#' @param   pInd        - indices of longitudinal samples to do prediction for
#' @param   locs.pred   - list with positions to predict
#' @param   Brandom.pred - random effect covaraites at prediction locations
#' @param   Bfixed.pred  - fixed effect covaraites at prediction locations
#' @param   quantiles   - list of posterior quantiles to compute
#' @param   return.samples - return samples used for prediction?
#' @param   type        - Type of prediction: Filter or Smoothing
#' @param mixedEffect_list -
#' @param   meas_noise  - the aviable noise classes: Normal or NIG
#' @param   noise       - the distribution of the mixed effect
#' @param   B_random    - list for the random effect covariates (needs to be matrix, can be NULL)
#' @param   B_fixed     - list for the fixed  effect covariates (needs to be matrix, can be NULL)
#' @param   beta_random - initial parameters of the random effect (mean parameter) (if not specified set to zero)
#' @param   beta_fixed  - initial parameters of the fixed  effect (if not specified set to zero)
#' @param   Sigma       - initial parameters of the covariance of random effect (if not specified set to I )
#' @param   nu          - shape parameter for noise (NIG only)
#' @param   mu          - shift parameter for noise (NIG only)
#' @param   U           - (list) inital guess of the random effect
#' @param   V           - (list) inital guess of the variance effect
#'
#' @param measurment_list   - list for measurement error:
#' @param     sigma       - measurement noise variance
#' @param     nu          - shape parameter for noise (NIG only)
#' @param     Vs          - (list) inital guess for the noise of the measurement
#' @param operator_list   - list created using create_operator function!
#' @param processes_list  - for the stochastic noise driving the
#' @param noise           - either Normal, NIG or GAL (change name to type rather then noise)
#' @param nu              - shape parameter for NIG or GAL
#' @param mu              - assymetric parameter for NIG or GAL
#'
#' @param nSim            - number of samples of the gibbs sampler
#' @param nBurnin         - number of samples to discard before starting prediction
#' @parma silent          - print iteration info
predictLong <- function( Y,
                         locs,
                         pInd,
                         locs.pred,
                         Brandom.pred,
                         Bfixed.pred,
                         return.samples = FALSE,
                         type = "Filter",
                         quantiles = NULL,
                         mixedEffect_list,
                         measurment_list,
                         processes_list,
                         operator_list,
                         nSim  = 1,
                         nBurnin = 10,   # steps before starting prediction
                         silent  = FALSE # print iteration info
)
{
  if(type=='Filter'){
    pred_type = 1
  } else if(type == 'Smoothing'){
    pred_type = 0
  } else {
    stop('Type needs to be either Filter or Smoothing.')
  }

  if(missing(locs.pred)){
    locs.pred <- locs
  }
  obs_list <- list()

  if(!missing(pInd) && !is.null(pInd)){
    Y                   <- Y[pInd]
    locs                <- locs[pInd]
    locs.pred           <- locs.pred[pInd]
    Brandom.pred        <- Brandom.pred[pInd]
    Bfixed.pred         <- Bfixed.pred[pInd]
    measurment_list$Vs  <- measurment_list$Vs[pInd]
    mixedEffect_list$B_fixed <- mixedEffect_list$B_fixed[pInd]
    mixedEffect_list$B_random <- mixedEffect_list$B_random[pInd]
    mixedEffect_list$U  <- mixedEffect_list$U
    processes_list$X    <- processes_list$X[pInd]
    processes_list$V    <- processes_list$V[pInd]
    n.patient = length(pInd)
  } else {
    if(is.list(Y)){
      n.patient = length(Y)
    } else {
      n.patient = 1
    }
  }

  for(i in 1:n.patient){
    if(is.list(locs)){
      n.pred.i = length(locs[[i]])
    } else {
      n.pred.i = length(locs)
    }

    if(type == "Filter"){
        pred.ind <- matrix(nrow = n.pred.i,ncol = 2)
        obs.ind  <- matrix(nrow = n.pred.i,ncol = 2)
        for(j in 1:(n.pred.i-1)){
          # pred.ind shows which values to save for the j:th prediction
          ind <- (1:length(locs.pred[[i]]))[(locs.pred[[i]] >= locs[[i]][j]) & (locs.pred[[i]] < locs[[i]][j+1])]
          pred.ind[j,] <- c(ind[1]-1,length(ind)) #first index and number of indices.
          # obs.ind shows which data to use for the j:th prediction
          obs.ind[j,] <- c(0,j)
        }
        ind <- (1:length(locs.pred[[i]]))[locs.pred[[i]] >= locs[[i]][n.pred.i]]
        pred.ind[n.pred.i,] <- c(ind[1]-1,length(ind))
        obs.ind[n.pred.i,] <- c(0,n.pred.i)
      } else {
        pred.ind <- matrix(c(0,length(locs.pred[[i]])),nrow = 1,ncol = 2)
        obs.ind  <- matrix(c(0,n.pred.i),nrow = 1,ncol = 2)
      }

      obs_list[[i]] <- list(A = spde.A(locs[[i]], operator_list$loc,
                                       right.boundary = operator_list$right.boundary,
                                       left.boundary = operator_list$left.boundary),
                            Apred = spde.A(locs.pred[[i]],operator_list$loc,
                                       right.boundary = operator_list$right.boundary,
                                       left.boundary = operator_list$left.boundary),
                            Y=Y[[i]],
                            pred_ind = pred.ind,
                            obs_ind = obs.ind,
                            locs = locs[[i]],
                            Brandom_pred = Brandom.pred[[i]],
                            Bfixed_pred = Bfixed.pred[[i]])
  }
  input <- list( obs_list         = obs_list,
                 operator_list    = operator_list,
                 measurementError_list  = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 processes_list   = processes_list,
                 nSim             = nSim,
                 nBurnin          = nBurnin,   # steps before starting gradient estimation
                 silent           = silent, # print iteration info)
                 pred_type        = pred_type
  )

  #output <- predictLong_cpp(input)
  output <- list()
  out_list <- list()

  if(return.samples){
    out_list$X.samples <- output$XVec
    out_list$W.samples <- output$WVec
  }

  out_list$locs <- locs.pred
  out_list$X.summary <- list()
  out_list$W.summary <- list()

  for(i in 1:length(locs)){

    out_list$X.summary[[i]] <- list()
    out_list$W.summary[[i]] <- list()
    out_list$X.summary[[i]]$Mean <- apply(output$XVec[[i]],1,mean)
    out_list$W.summary[[i]]$Mean <- apply(output$WVec[[i]],1,mean)
    out_list$X.summary[[i]]$Var  <- apply(output$XVec[[i]],1,var)
    out_list$W.summary[[i]]$Var  <- apply(output$WVec[[i]],1,var)
    out_list$X.summary[[i]]$Median <- apply(output$XVec[[i]],1,median)
    out_list$W.summary[[i]]$Median <- apply(output$WVec[[i]],1,median)

    if(!is.null(quantiles)){
      x.list <- list()
      w.list <- list()
      for(c in 1:length(quantiles)){
        c.i <- list()
        c.i$level = quantiles[c]
        c.i$field <- apply(output$XVec[[i]],1,quantile,probs=c(quantiles[c]))
        x.list[[c]] = c.i
        c.i$field <- apply(output$WVec[[i]],1,quantile,probs=c(quantiles[c]))
        w.list[[c]] = c.i
      }
      out_list$X.summary[[i]]$quantiles <- x.list
      out_list$W.summary[[i]]$quantiles <- w.list
    }
  }
  return(out_list)
}