#' @param   pInd        - indices of longitudinal samples to do prediction for
#' @param   locs.pred   - list with positions to predict
#' @param   Brandom.pred - random effect covaraites at prediction locations
#' @param   Bfixed.pred  - fixed effect covaraites at prediction locations
#' @param   quantiles   - list of posterior quantiles to compute
#' @param   excursions   - list of excursion probabilities to compute. Each list should contain:
#'                type  - type of excursion '>' or '<'.
#'                level - level to compute excursion probability for
#'                process - which process to compute the probability for, 'X', 'W', 'Y','Xderivative' or 'Wderivative'
#'
#' @param   return.samples - return samples used for prediction?
#' @param   type        - Type of prediction: Filter or Smoothing
# All other parameters explained in help text for estimateGH.R
predictLong <- function( Y,
                         locs,
                         pInd,
                         locs.pred,
                         Brandom.pred,
                         Bfixed.pred,
                         return.samples = FALSE,
                         type = "Filter",
                         quantiles = NULL,
                         excursions = NULL,
                         crps = FALSE,
                         mixedEffect_list,
                         measurment_list,
                         processes_list,
                         operator_list,
                         nSim  = 1,
                         nBurnin = 10,   # steps before starting prediction
                         silent  = FALSE, # print iteration info
                         max.num.threads = 2
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
  common.grid = FALSE
  if(length(operator_list$loc)==1){
    common.grid = TRUE
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
    if(!common.grid){
      operator_list$Q <- operator_list$Q[pInd]
      operator_list$loc <- operator_list$loc[pInd]
      operator_list$h <- operator_list$h[pInd]
      if(operator_list$type == "Matern"){
        operator_list$C <- operator_list$C[pInd]
        operator_list$G <- operator_list$G[pInd]
      }
    }
    n.patient = length(pInd)
  } else {
    if(is.list(Y)){
      n.patient = length(Y)
    } else {
      n.patient = 1
    }
  }
  ind1 <- seq(1,nSim,by=2)
  ind2 <- seq(2,nSim,by=2)

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
      if(length(Y[[i]]) != length(locs[[i]])){
        stop("Length of Y and locs differ.")
      }
      if(common.grid){
        obs_list[[i]] <- list(A = spde.A(locs[[i]], operator_list$loc[[1]],
                                         right.boundary = operator_list$right.boundary,
                                         left.boundary = operator_list$left.boundary),
                              Apred = spde.A(locs.pred[[i]],operator_list$loc[[1]],
                                             right.boundary = operator_list$right.boundary,
                                             left.boundary = operator_list$left.boundary),
                              Y=Y[[i]],
                              pred_ind = pred.ind,
                              obs_ind = obs.ind,
                              locs = locs[[i]],
                              Brandom_pred = Brandom.pred[[i]],
                              Bfixed_pred = Bfixed.pred[[i]])
      } else {
        obs_list[[i]] <- list(A = spde.A(locs[[i]], operator_list$loc[[i]],
                                         right.boundary = operator_list$right.boundary,
                                         left.boundary = operator_list$left.boundary),
                              Apred = spde.A(locs.pred[[i]],operator_list$loc[[i]],
                                             right.boundary = operator_list$right.boundary,
                                             left.boundary = operator_list$left.boundary),
                              Y=Y[[i]],
                              pred_ind = pred.ind,
                              obs_ind = obs.ind,
                              locs = locs[[i]],
                              Brandom_pred = Brandom.pred[[i]],
                              Bfixed_pred = Bfixed.pred[[i]])
      }

  }
  input <- list( obs_list         = obs_list,
                 operator_list    = operator_list,
                 measurementError_list  = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 processes_list   = processes_list,
                 nSim             = nSim,
                 nBurnin          = nBurnin,   # steps before starting gradient estimation
                 silent           = silent, # print iteration info)
                 pred_type        = pred_type,
                 n_threads        = max.num.threads
  )

  output <- predictLong_cpp(input)
  out_list <- list()

  if(return.samples){
    out_list$Y.samples <- output$YVec
    out_list$X.samples <- output$XVec
    out_list$W.samples <- output$WVec
  }

  out_list$locs <- locs.pred
  out_list$Y.summary <- list()
  out_list$X.summary <- list()
  out_list$W.summary <- list()

  for(i in 1:length(locs)){
    out_list$Y.summary[[i]] <- list()
    out_list$X.summary[[i]] <- list()
    out_list$W.summary[[i]] <- list()

    out_list$Y.summary[[i]]$Mean <- apply(output$YVec[[i]],1,mean)
    out_list$X.summary[[i]]$Mean <- apply(output$XVec[[i]],1,mean)
    out_list$W.summary[[i]]$Mean <- apply(output$WVec[[i]],1,mean)

    out_list$Y.summary[[i]]$Var  <- apply(output$YVec[[i]],1,var)
    out_list$X.summary[[i]]$Var  <- apply(output$XVec[[i]],1,var)
    out_list$W.summary[[i]]$Var  <- apply(output$WVec[[i]],1,var)

    out_list$Y.summary[[i]]$Median <- apply(output$YVec[[i]],1,median)
    out_list$X.summary[[i]]$Median <- apply(output$XVec[[i]],1,median)
    out_list$W.summary[[i]]$Median <- apply(output$WVec[[i]],1,median)

    if(!is.null(quantiles)){
      y.list <- list()
      x.list <- list()
      w.list <- list()
      for(c in 1:length(quantiles)){
        c.i <- list()
        c.i$level = quantiles[c]
        c.i$field <- apply(output$YVec[[i]],1,quantile,probs=c(quantiles[c]))
        y.list[[c]] = c.i
        c.i$field <- apply(output$XVec[[i]],1,quantile,probs=c(quantiles[c]))
        x.list[[c]] = c.i
        c.i$field <- apply(output$WVec[[i]],1,quantile,probs=c(quantiles[c]))
        w.list[[c]] = c.i
      }
      out_list$Y.summary[[i]]$quantiles <- y.list
      out_list$X.summary[[i]]$quantiles <- x.list
      out_list$W.summary[[i]]$quantiles <- w.list
    }
    if(!is.null(excursions)){
      for(c in 1:length(excursions)){
        ex.i <- list(type = excursions[[c]]$type, level = excursions[[c]]$level)

        if(excursions[[c]]$process == 'X'){
          proc <- output$XVec[[i]]
        } else if(excursions[[c]]$process == 'Y'){
          proc <- output$YVec[[i]]
        } else if(excursions[[c]]$process == 'W'){
          proc <- output$WVec[[i]]
        } else if(excursions[[c]]$process == 'Xderivative'){
          proc <-t(t(apply(output$XVec[[i]],2,diff))/diff(locs.pred[[i]]))
        } else if(excursions[[c]]$process == 'Wderivative'){
          proc <-t(t(apply(output$WVec[[i]],2,diff))/diff(locs.pred[[i]]))
        }
        if(ex.i$type == '>'){
          if(length(proc)>1){
            ex.i$P <- apply(proc>ex.i$level,1,mean)
          } else {
            ex.i$P <- mean(proc>ex.i$level)
          }
        } else {
          if(length(proc)>1){
            ex.i$P <- apply(proc<ex.i$level,1,mean)
          } else {
            ex.i$P <- mean(proc<ex.i$level)
          }
        }
        if(excursions[[c]]$process == 'X' || excursions[[c]]$process == 'Xderivative'){
          out_list$X.summary[[i]]$excursions <- ex.i
        } else if(excursions[[c]]$process == 'Y'){
          out_list$Y.summary[[i]]$excursions <- ex.i
        } else if(excursions[[c]]$process == 'W' || excursions[[c]]$process == 'Wderivative'){
          out_list$W.summary[[i]]$excursions <- ex.i
        }
      }
    }
    if(crps){
      out_list$Y.summary[[i]]$crps <- apply(abs(matrix(rep(Y[[i]],each=length(ind1)),ncol=length(ind1),byrow=TRUE)-output$YVec[[i]][,ind1]),1,mean) - 0.5*apply(abs(output$YVec[[i]][,ind1]-output$YVec[[i]][,ind2]),1,mean)
    }
  }
  return(out_list)
}
