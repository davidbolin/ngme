#'
#' @title Obtain predictions
#'
#' @description A function to obtain predictions based on either filtering
#'    or smoothing distributions.
#'
#' @param pInd A numeric vector that contains the indices of longitudinal
#'    subjects for whom the predictions are to be obtained.
#' @param locs.pred A numeric list that contains the timings of the repeated
#'    measurements.
#' @param Brandom.pred A numeric list that contains random effects covaraite
#'    matrices.
#' @param Bfixed.pred  A numeric list that contains fixed effects covaraite
#'    matrices.
#' @return.samples A logical variable for returning the
#'    Monte Carlo samples used to compute the predictions; \code{"TRUE"} indicates
#'    return, \code{"FALSE"} do not return.
#' @param type A character string for the type of prediction: \code{"Filter"} for
#'   filtering, \code{"Smoothing"} for smoothing.
#' @param quantiles A two-elemnent vector that contains the quantiles
#'   of the predictions to be calculated.
#' @param predict.derivatives STUFF
#' @param Y.val Observations to use when calculating CRPS
#' @param excursions A list of excursion probabilities to compute.
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
#' @param crps A logical variable for calculating
#'    continuous ranked probability score (CRPS); \code{"TRUE"} indicates
#'    calculate, \code{"FALSE"} do not calculate.
#' @param crps.skip A numerical value, say a, that indicates every \emph{a}th
#'    element of the sample to be used to compute the crps score.
#' @inheritParams estimateLong
#' @param max.num.threads STUFF
#' @param repeat.mix STUFF
#'
#' @return A list of output.
#'
#' @details This function calls \code{"predictLong_cpp"} internally.
#'    It is wrapped by \code{"predict.ngme"}, and not advised to be used.
#'
#' @seealso \code{\link{predict.ngme}}
#' @examples
#'   \dontrun{
#'   predictLong(...)
#'   }

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
                         Y.val,
                         mixedEffect_list,
                         measurment_list,
                         processes_list,
                         operator_list = NULL,
                         nSim  = 1,
                         nBurnin = 10,   # steps before starting prediction
                         silent  = FALSE, # print iteration info
                         max.num.threads = 2,
                         repeat.mix = 10,
                         seed    = NULL
)
{
  ind.general = 0
  if(type == "LOOCV"){
    pred_type = 3
    ind.general = 1
  }else if(type == "Nowcast"){
    pred_type = 2
  }else if(type=='Filter'){
    pred_type = 1
  } else if(type == 'Smoothing'){
    pred_type = 0
  } else {
    stop('Type needs to be either LOOCV, Filter, Smoothing, or Nowcast.')
  }
  use.random.effect = TRUE
  if(is.null(mixedEffect_list$B_random)){
    use.random.effect = FALSE
  }
  if(!use.random.effect && !missing(Brandom.pred)){
    stop("Model not specified using random effects")
  }
  if(type=="Smoothing" && crps && missing(Y.val)){
    warning("CRPS is calculated for smoothing prediction without specifying Y.val. Are you sure this is correct?")
  }
  use.process = TRUE
  if(missing(processes_list) || is.null(processes_list)){
    use.process = FALSE
  }

  if(!missing(processes_list) && is.null(processes_list)){
    stop("Operator list missing")
  }
  
  bivariate = FALSE
  if(use.process && operator_list$manifold == "R2"){
    if(type == "Nowcast" || type == "Filter"){
      stop("Nowcasting and filtering can only be done for temporal models.")
    }
    if(!is.null(predict.derivatives)){
      stop("Derivatives can only be predicted for temporal models")
    }
    if(operator_list$type == "matern bivariate"){
      bivariate = TRUE
      ind.general = 1
    }
  }
  if(missing(locs.pred)){
    locs.pred <- locs
  } else {
    if(pred_type == 3){
      warning("locs.pred supplied for LOOCV.")
    }
  }
  
  if(missing(Bfixed.pred)){
    if(pred_type == 3){
      Bfixed.pred = mixedEffect_list$B_fixed
    } else {
        stop("Must supply Bfixed.pred")
    }
  } else {
    if(pred_type == 3){
      warning("Bfixed.pred supplied for LOOCV.")
    }
  }

  if(use.random.effect){
    if(missing(Brandom.pred)){
      if(pred_type == 3){
        Brandom.pred = mixedEffect_list$B_random
      } else {
        stop("Must supply Brandom.pred")
      }
    } else {
      if(pred_type == 3){
        warning("Brandom.pred supplied for LOOCV.")
      }
    }  
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

      operator_list$Q <- operator_list$Q[pInd]
      operator_list$loc <- operator_list$loc[pInd]
      operator_list$h <- operator_list$h[pInd]
      if(tolower(operator_list$type) %in% c("matern","exponential","matern bivariate","matern.asym")){
        operator_list$C <- operator_list$C[pInd]
        operator_list$Ci <- operator_list$Ci[pInd]
        operator_list$G <- operator_list$G[pInd]
      }
    }

    n.patient = length(pInd)
  } else {
    
    if(is.list(Y)){
      n.patient = length(Y)
      pInd = 1:length(Y)
    } else {
      n.patient = 1
      pInd = 1
    }
  }

  if(sum(abs(unlist(lapply(locs.pred,length))-unlist(lapply(lapply(locs.pred,unique),length))))>0){
    stop("Prediction locations should be unique")
  }
  if(!is.null(predict.derivatives)){
    locs.pred.delta <- locs.pred
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
    if(type == "LOOCV"){
      #for CV, pred.ind and obs.ind have to contain the actual incides
      if(is.matrix(locs.pred[[i]])){
        no <- dim(locs.pred[[i]])[1]
      } else {
        no <- length(locs.pred[[i]])
      }
      if(bivariate){
        pred.ind <- cbind(diag(no,diag(no)))
        obs.ind <- cbind(1 - diag(no),1 - diag(no))
        obs.ind <- obs.ind[,!is.nan(c(Y[[i]]))] #remove missing observations
      } else {
        pred.ind <- diag(no)
        obs.ind <- 1 - diag(no)  
      }
      
    } else if(type== "Nowcast"){
      pred.ind <- obs.ind <- NULL
      ind <- (1:length(locs.pred[[i]]))[locs.pred[[i]] < li[1]] #all prediction locations before first obs
      if(length(ind)>0){
        pred.ind <- rbind(pred.ind,c(0,length(ind)))
        obs.ind <- rbind(obs.ind,c(0,0))
      }
      # pred.ind shows which values to save for the j:th prediction
      # obs.ind shows which data to use for the j:th prediction
      if(n.pred.i>1){
        for(j in 2:n.pred.i){
          #indices for prediction locatiocs in  [s_(j-1), s_j), should use data up to and including s_(j-1)
          ind <- (1:length(locs.pred[[i]]))[(locs.pred[[i]] >= li[j-1]) & (locs.pred[[i]] < li[j])]
          if(length(ind)>0){
            pred.ind <- rbind(pred.ind,c(ind[1]-1,length(ind))) #first index and number of indices.
            obs.ind <- rbind(obs.ind,c(0,j-1))
          }
        }
      }
      #do the prediction for all locations in [s_N,max(loc.pred)]
      if(max(locs.pred[[i]])>=max(li)){
        ind <- (1:length(locs.pred[[i]]))[locs.pred[[i]] >= max(li)]
        pred.ind <- rbind(pred.ind,c(ind[1]-1,length(ind)))
        obs.ind <- rbind(obs.ind,c(0,length(li)))
      }
      n.pred.i = dim(pred.ind)[1]
    } else if(type == "Filter"){
        pred.ind <- obs.ind <- NULL
        ind <- (1:length(locs.pred[[i]]))[locs.pred[[i]] <= li[1]] #all prediction locations up to the first obs
        if(length(ind)>0){ #If there are values to predict up to the first obs
          pred.ind <- rbind(pred.ind,c(0,length(ind)))
          obs.ind <- rbind(obs.ind,c(0,0))
        }
        # pred.ind shows which values to save for the j:th prediction
        # obs.ind shows which data to use for the j:th prediction
        if(n.pred.i>1){
          for(j in 2:n.pred.i){
            ind <- (1:length(locs.pred[[i]]))[(locs.pred[[i]] > li[j-1]) & (locs.pred[[i]] <= li[j])]
            if(length(ind)>0){
              pred.ind <- rbind(pred.ind,c(ind[1]-1,length(ind))) #first index and number of indices.
              obs.ind <- rbind(obs.ind,c(0,j-1))
            }
          }
        }
        #obs.ind[n.pred.i,] <- c(0,n.pred.i-1)
        if(max(locs.pred[[i]])>max(li)){
          ind <- (1:length(locs.pred[[i]]))[locs.pred[[i]] > max(li)]
          pred.ind <- rbind(pred.ind,c(ind[1]-1,length(ind)))
          obs.ind <- rbind(obs.ind,c(0,length(li)))
        }
        n.pred.i = dim(pred.ind)[1]
      } else { #Smoothing
        if(is.matrix(locs.pred[[i]])){
          pred.ind <- matrix(c(0,dim(locs.pred[[i]])[1]),nrow = 1,ncol = 2)
        } else {
          pred.ind <- matrix(c(0,length(locs.pred[[i]])),nrow = 1,ncol = 2)
        }
        obs.ind  <- matrix(c(0,n.pred.i),nrow = 1,ncol = 2)  
        if(bivariate){
          pred.ind <- matrix(rep(1,2*dim(locs.pred[[i]])[1]),nrow=1)
          obs.ind <- matrix(rep(1,2*n.pred.i),nrow=1)
          obs.ind <- obs.ind[,!is.nan(c(Y[[i]])),drop=FALSE] #remove missing observations
        }
      }
  
      obs_list[[i]] <- list(Y=c(Y[[i]]),
                            pred_ind = pred.ind,
                            obs_ind = obs.ind,
                            locs = locs[[i]],
                            Bfixed_pred = Bfixed.pred[[i]])
      if(use.process){
        A <- build.A.matrix(operator_list,locs,i)
        Ap <- build.A.matrix(operator_list,locs.pred,i)
        if(bivariate){
          A1 <- A[!is.nan(Y[[i]][,1]),]
          A2 <- A[!is.nan(Y[[i]][,2]),]
          obs_list[[i]]$A = bdiag(A1,A2)
          obs_list[[i]]$Apred = bdiag(Ap,Ap)
          Yi <- c(Y[[i]])
          obs_list[[i]]$Y = Yi[!is.nan(Yi)]
          obs_list[[i]]$locs = cbind(locs[[i]],locs[[i]])
        } else {
          obs_list[[i]]$A = A
          obs_list[[i]]$Apred = Ap
        }
      }
      
      if(use.random.effect){
        obs_list[[i]]$Brandom_pred = Brandom.pred[[i]]
      }

      if(!is.null(predict.derivatives)){
        if(use.process){
          locs.pred.delta[[i]] <- locs.pred.delta[[i]]+predict.derivatives$delta
          obs_list[[i]]$Apred1 = build.A.matrix(operator_list,locs.pred.delta,i)
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
                 predict_derivative = predict_derivative,
                 ind_general = ind.general)
  if(use.process){
    input$processes_list   = processes_list
    input$operator_list    = operator_list
  }
  if(is.null(seed) == FALSE)
    input <- setseed_ME(input, seed)

  output <- predictLong_cpp(input)
  out_list <- list()


  if(return.samples){
    out_list$Y.samples <- output$YVec
    out_list$X.samples <- output$XVec
    out_list$W.samples <- output$WVec
    out_list$V.samples <- output$VVec
    if(type == "Smoothing"){
      out_list$U.samples <- output$UVec
      if(use.process){
        out_list$Wnoise.samples <- output$WnoiseVec
        out_list$Wnoise.summary <- list()
      }
    }
    if(!is.null(predict.derivatives)){
      out_list$Xderivative.samples <- output$XVec_deriv
      out_list$Wderivative.samples <- output$WVec_deriv
    }
  }
  out_list$pInd <- pInd
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
    if(type == "Smoothing"){
      out_list$U.summary[[i]] <- list()
      out_list$U.summary[[i]]$Mean <- apply(output$UVec[[i]],1,mean)
      out_list$U.summary[[i]]$Var  <- apply(output$UVec[[i]],1,var)
      out_list$U.summary[[i]]$Median <- apply(output$UVec[[i]],1,median)
      if(use.process){
        out_list$Wnoise.summary[[i]] <- list()
        out_list$Wnoise.summary[[i]]$Mean <- apply(output$WnoiseVec[[i]],1,mean)
        out_list$Wnoise.summary[[i]]$Var  <- apply(output$WnoiseVec[[i]],1,var)
        out_list$Wnoise.summary[[i]]$Median <- apply(output$WnoiseVec[[i]],1,median)
      }
    }

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
      wnoise.list <- list()
      v.list <- list()
      u.list <- list()
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
        if(type=="Smoothing"){
          c.i$field <- apply(output$UVec[[i]],1,quantile,probs=c(quantiles[c]))
          u.list[[c]] = c.i
          if(use.process){
            c.i$field <- apply(output$WnoiseVec[[i]],1,quantile,probs=c(quantiles[c]))
            wnoise.list[[c]] = c.i
          }
        }
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
      if(type=="Smoothing"){
        out_list$U.summary[[i]]$quantiles <- u.list
        if(use.process){
          out_list$Wnoise.summary[[i]]$quantiles <- wnoise.list
        }
      }
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
      ind1 = 1:round(nSim/2)
      ind2 <- 1+(nSim/2+ind1-1)%%nSim
      if(dim(output$YVec[[i]])[1]>1){
        out_list$Y.summary[[i]]$crps <- apply(abs(matrix(rep(Y.val[[i]],each=length(ind1)),ncol=length(ind1),byrow=TRUE)-output$YVec[[i]][,ind1]),1,mean) - 0.5*apply(abs(output$YVec[[i]][,ind1]-output$YVec[[i]][,ind2]),1,mean)
      } else {
        out_list$Y.summary[[i]]$crps <- mean(abs(matrix(rep(Y.val[[i]],each=length(ind1)),ncol=length(ind1),byrow=TRUE)-output$YVec[[i]][,ind1])) - 0.5*mean(abs(output$YVec[[i]][,ind1]-output$YVec[[i]][,ind2]))
      }
    }
  }
  return(out_list)
}

#' @title STUFF
#'
#' @description STUFF
#'
#' @param mixedEffect_list A list of inputs for the random effects.
#' @param processes_list A list of inputs for the process.
#' @param operator_list A list of inputs for the operator.
#' @param measurement_list A list of inputs for the measurement error.
#' @param locs A numeric list that contains the timings at which the outcomes
#'    are collected.
#' @param Brandom A numeric list of random effects covariate matrices.
#' @param Bfixed A numeric list of fixed effects covariate matrices.
#'
#' @return A list of output.
#'
#' @details STUFF
#'
#' @examples
#'   \dontrun{
#'   updateLists(...)
#'   }

updateLists <- function(mixedEffect_list,
                        processes_list,
                        operator_list,
                        measurement_list,
                        Bfixed,
                        Brandom,
                        locs,
                        pred.ind = NULL)
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
      X[[i]] <- rep(0, length(operator_list$h[[pred.ind[i]]]))
      V[[i]] <- operator_list$h[[pred.ind[i]]]
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
