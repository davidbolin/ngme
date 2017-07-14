#' @title STUFF
#' 
#' @description STUFF
#' 
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
                         Brandom.pred = NULL,
                         Bfixed.pred,
                         return.samples = FALSE,
                         type = "Filter",
                         quantiles = NULL,
                         predict.derivatives = NULL,
                         excursions = NULL,
                         crps = FALSE,
                         crps.skip = 10,
                         mixedEffect_list,
                         measurment_list,
                         processes_list,
                         operator_list = NULL,
                         nSim  = 1,
                         nBurnin = 10,   # steps before starting prediction
                         silent  = FALSE, # print iteration info
                         max.num.threads = 2,
                         repeat.mix = 10
)
{

  if(type=='Filter'){
    pred_type = 1
  } else if(type == 'Smoothing'){
    pred_type = 0
  } else {
    stop('Type needs to be either Filter or Smoothing.')
  }
  use.random.effect = TRUE
  if(is.null(mixedEffect_list$B_random)){
    use.random.effect = FALSE
  }
  if(!use.random.effect && !missing(Brandom.pred)){
    stop("Model not specified using random effects")
  }

  use.process = TRUE
  if(missing(processes_list) || is.null(processes_list)){
    use.process = FALSE
  }

  if(!missing(processes_list) && is.null(processes_list)){
    stop("Operator list missing")
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
    Bfixed.pred         <- Bfixed.pred[pInd]
    measurment_list$Vs  <- measurment_list$Vs[pInd]
    mixedEffect_list$B_fixed <- mixedEffect_list$B_fixed[pInd]
    if(use.random.effect){
      Brandom.pred        <- Brandom.pred[pInd]
      mixedEffect_list$B_random <- mixedEffect_list$B_random[pInd]
      if(!is.null(mixedEffect_list$U)){
        mixedEffect_list$U  <- mixedEffect_list$U
      }
    }
    if(!is.null(predict.derivatives)){
      predict.derivatives$B_fixed <- predict.derivatives$B_fixed[pInd]
      if(use.random.effect){
        predict.derivatives$B_random <- predict.derivatives$B_random[pInd]
      }
    }
    if(use.process){
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
    }

    n.patient = length(pInd)
  } else {
    if(is.list(Y)){
      n.patient = length(Y)
    } else {
      n.patient = 1
    }
  }

  ind1 <- 1:nSim
  ind2 <- 1+(nSim/2+ind1-1)%%nSim

  ind3 <- seq(1,nSim,by=2)
  ind4 <- seq(2,nSim,by=2)

  if(sum(abs(unlist(lapply(locs.pred,length))-unlist(lapply(lapply(locs.pred,unique),length))))>0){
    stop("Prediction locations should be unique")
  }

  for(i in 1:n.patient){
    if(is.list(locs)){
      li = locs[[i]]
    } else {
      li = locs
    }
    if(is.matrix(li)){
      n.pred.i = dim(li)[1]
    } else {
      n.pred.i = length(li)
    }

    if(type == "Filter"){
        pred.ind <- obs.ind <- NULL
        ind <- (1:length(locs.pred[[i]]))[locs.pred[[i]] <= li[1]] #all prediction locations before first obs
        if(length(ind)>0){
          pred.ind <- Matrix::rBind(pred.ind,c(0,length(ind)))
          obs.ind <- Matrix::rBind(obs.ind,c(0,0))
        }
        # pred.ind shows which values to save for the j:th prediction
        # obs.ind shows which data to use for the j:th prediction
        if(n.pred.i>1){
          for(j in 2:n.pred.i){
            ind <- (1:length(locs.pred[[i]]))[(locs.pred[[i]] > li[j-1]) & (locs.pred[[i]] <= li[j])]
            if(length(ind)>0){
              pred.ind <- Matrix::rBind(pred.ind,c(ind[1]-1,length(ind))) #first index and number of indices.
              obs.ind <- Matrix::rBind(obs.ind,c(0,j-1))
            }
          }
        }
        #obs.ind[n.pred.i,] <- c(0,n.pred.i-1)
        if(max(locs.pred[[i]])>max(li)){
          ind <- (1:length(locs.pred[[i]]))[locs.pred[[i]] >= max(li)]
          pred.ind <- Matrix::rBind(pred.ind,c(0,length(ind)))
          obs.ind <- Matrix::rBind(obs.ind,c(0,length(locs)-1))
        }
        n.pred.i = dim(pred.ind)[1]
      } else {
        if(is.matrix(locs.pred[[i]])){
          pred.ind <- matrix(c(0,dim(locs.pred[[i]])[1]),nrow = 1,ncol = 2)
        } else {
          pred.ind <- matrix(c(0,length(locs.pred[[i]])),nrow = 1,ncol = 2)
        }
        obs.ind  <- matrix(c(0,n.pred.i),nrow = 1,ncol = 2)
      }

      obs_list[[i]] <- list(Y=Y[[i]],
                            pred_ind = pred.ind,
                            obs_ind = obs.ind,
                            locs = locs[[i]],
                            Bfixed_pred = Bfixed.pred[[i]])
      if(use.process){
        obs_list[[i]]$A = build.A.matrix(operator_list,locs,i)
        obs_list[[i]]$Apred = build.A.matrix(operator_list,locs.pred,i)
      }

      if(use.random.effect){
        obs_list[[i]]$Brandom_pred = Brandom.pred[[i]]
      }

      if(!is.null(predict.derivatives)){
        if(use.process){
          obs_list[[i]]$Apred = build.A.matrix(operator_list,locs.pred+predict.derivatives$delta,i)
        }
        obs_list[[i]]$Bfixed_pred1 = predict.derivatives$Bfixed[[i]]
        if(use.random.effect){
          obs_list[[i]]$Brandom_pred1 = predict.derivatives$Brandom[[i]]
        }
    }
  }
  delta = 1
  predict_derivative = 0
  if(!is.null(predict.derivatives)){
    delta = predict.derivatives$delta
    predict_derivative = 1
  }
  input <- list( obs_list         = obs_list,
                 measurementError_list  = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 nSim             = nSim,
                 nBurnin          = nBurnin,   # steps before starting gradient estimation
                 silent           = silent, # print iteration info)
                 pred_type        = pred_type,
                 n_threads        = max.num.threads,
                 mix_samp = repeat.mix,
                 use_random_effect = use.random.effect,
                 derivative_scaling = delta,
                 predict_derivative = predict_derivative)
  if(use.process){
    input$processes_list   = processes_list
    input$operator_list    = operator_list
  }

  output <- predictLong_cpp(input)
  out_list <- list()

  if(return.samples){
    out_list$Y.samples <- output$YVec
    out_list$X.samples <- output$XVec
    out_list$W.samples <- output$WVec
    out_list$V.samples <- output$VVec
    if(!is.null(predict.derivatives)){
      out_list$Xderivative.samples <- output$XVec_deriv
      out_list$Wderivative.samples <- output$WVec_deriv
    }
  }

  out_list$locs <- locs.pred
  out_list$Y.summary <- list()
  out_list$X.summary <- list()
  out_list$W.summary <- list()
  out_list$V.summary <- list()
  if(!is.null(predict.derivatives)){
    out_list$Xderivative.summary <- list()
    out_list$Wderivative.summary <- list()
  }

  for(i in 1:length(locs)){
    out_list$Y.summary[[i]] <- list()
    out_list$X.summary[[i]] <- list()
    out_list$W.summary[[i]] <- list()
    out_list$V.summary[[i]] <- list()

    out_list$Y.summary[[i]]$Mean <- apply(output$YVec[[i]],1,mean)
    out_list$X.summary[[i]]$Mean <- apply(output$XVec[[i]],1,mean)
    out_list$W.summary[[i]]$Mean <- apply(output$WVec[[i]],1,mean)
    out_list$V.summary[[i]]$Mean <- apply(output$VVec[[i]],1,mean)

    out_list$Y.summary[[i]]$Var  <- apply(output$YVec[[i]],1,var)
    out_list$X.summary[[i]]$Var  <- apply(output$XVec[[i]],1,var)
    out_list$W.summary[[i]]$Var  <- apply(output$WVec[[i]],1,var)
    out_list$V.summary[[i]]$Var <- apply(output$VVec[[i]],1,var)

    out_list$Y.summary[[i]]$Median <- apply(output$YVec[[i]],1,median)
    out_list$X.summary[[i]]$Median <- apply(output$XVec[[i]],1,median)
    out_list$W.summary[[i]]$Median <- apply(output$WVec[[i]],1,median)
    out_list$V.summary[[i]]$Median <- apply(output$VVec[[i]],1,median)

    if(!is.null(predict.derivatives)){
      out_list$Xderivative.summary[[i]] <- list()
      out_list$Wderivative.summary[[i]] <- list()
      out_list$Xderivative.summary[[i]]$Mean <- apply(output$XVec_deriv[[i]],1,mean)
      out_list$Wderivative.summary[[i]]$Mean <- apply(output$WVec_deriv[[i]],1,mean)
      out_list$Xderivative.summary[[i]]$Var  <- apply(output$XVec_deriv[[i]],1,var)
      out_list$Wderivative.summary[[i]]$Var  <- apply(output$WVec_deriv[[i]],1,var)
      out_list$Xderivative.summary[[i]]$Median <- apply(output$XVec_deriv[[i]],1,median)
      out_list$Wderivative.summary[[i]]$Median <- apply(output$WVec_deriv[[i]],1,median)
    }

    if(!is.null(quantiles)){
      y.list <- list()
      x.list <- list()
      w.list <- list()
      v.list <- list()
      if(!is.null(predict.derivatives)){
        xd.list <- list()
        wd.list <- list()
      }
      for(c in 1:length(quantiles)){
        c.i <- list()
        c.i$level = quantiles[c]
        c.i$field <- apply(output$YVec[[i]],1,quantile,probs=c(quantiles[c]))
        y.list[[c]] = c.i
        c.i$field <- apply(output$XVec[[i]],1,quantile,probs=c(quantiles[c]))
        x.list[[c]] = c.i
        c.i$field <- apply(output$WVec[[i]],1,quantile,probs=c(quantiles[c]))
        w.list[[c]] = c.i
        c.i$field <- apply(output$VVec[[i]],1,quantile,probs=c(quantiles[c]))
        v.list[[c]] = c.i
        if(!is.null(predict.derivatives)){
          c.i$field <- apply(output$XVec_deriv[[i]],1,quantile,probs=c(quantiles[c]))
          xd.list[[c]] = c.i
          c.i$field <- apply(output$WVec_deriv[[i]],1,quantile,probs=c(quantiles[c]))
          wd.list[[c]] = c.i
        }
      }
      out_list$Y.summary[[i]]$quantiles <- y.list
      out_list$X.summary[[i]]$quantiles <- x.list
      out_list$W.summary[[i]]$quantiles <- w.list
      out_list$V.summary[[i]]$quantiles <- v.list
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
          proc <- output$XVec_deriv[[i]]
        } else if(excursions[[c]]$process == 'Wderivative'){
          proc <- output$WVec_deriv[[i]]
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
        if(excursions[[c]]$process == 'X'){
          out_list$X.summary[[i]]$excursions <- ex.i
        } else if(excursions[[c]]$process == 'Y'){
          out_list$Y.summary[[i]]$excursions <- ex.i
        } else if(excursions[[c]]$process == 'W'){
          out_list$W.summary[[i]]$excursions <- ex.i
        } else if(excursions[[c]]$process == 'Xderivative'){
          out_list$Xderivative.summary[[i]]$excursions <- ex.i
        } else if(excursions[[c]]$process == 'Wderivative'){
          out_list$Wderivative.summary[[i]]$excursions <- ex.i
        }
      }
    }
    if(crps){
      if(dim(output$YVec[[i]])[1]>1){
        out_list$Y.summary[[i]]$crps <- apply(abs(matrix(rep(Y[[i]],each=length(ind1)),ncol=length(ind1),byrow=TRUE)-output$YVec[[i]][,ind1]),1,mean) - 0.5*apply(abs(output$YVec[[i]][,ind1]-output$YVec[[i]][,ind2]),1,mean)
      } else {
        out_list$Y.summary[[i]]$crps <- mean(abs(matrix(rep(Y[[i]],each=length(ind1)),ncol=length(ind1),byrow=TRUE)-output$YVec[[i]][,ind1])) - 0.5*mean(abs(output$YVec[[i]][,ind1]-output$YVec[[i]][,ind2]))
      }
    }
  }
  return(out_list)
}

#' @title STUFF
#' 
#' @description STUFF
#' 
updateLists <- function(mixedEffect_list,
                        processes_list,
                        operator_list,
                        measurement_list,
                        Bfixed,
                        Brandom,
                        locs)
{
  n.pred <- length(Bfixed)

  mixedEffect_list$B_fixed <- Bfixed
  if(!missing(Brandom)){
    mixedEffect_list$B_random <- Brandom
  }

  X <- V <- Vin <- list()
  for(i in 1:n.pred){
    Vin[[i]] <- rep(1,length(locs[[i]]))
    if(length(operator_list$h) == 1){
      X[[i]] <- rep(0, length(operator_list$h[[1]]))
      V[[i]] <- operator_list$h[[1]]
    } else {
      error("Not yet implemented")
    }
  }

  processes_list$X <- X
  if(!is.null(processes_list$V))
    processes_list$V <- V

  if(measurement_list$noise != "Normal")
    measurement_list$Vs <- Vin

  return(list(processes_list = processes_list,
              measurement_list = measurement_list,
              mixedEffect_list = mixedEffect_list))
}
