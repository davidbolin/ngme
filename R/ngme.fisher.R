
#' @title Estimation of Fisher information matrix
#'
#' @description Estimates Fisher information matrix for ngme result.
#'
#' @param fit An ngme object. 
#' @param std.threshold A threshold for the MC standard deviation of the estimates. The estimation is run until
#' all diagonal elements F.i of the inverse Fisher information matrix satisfy std(F.i)*std.threshold < F.i
#' @param only.effects If TRUE, the criteria for std.threshold is only applied to the fixed effect estimates.
#' @param observed If TRUE, the observed Fisher information matrix is estimated. Otherwise the ordinary Fisher-Information matrix is estimated.
#' @param nIter A numeric value for the number of iterations to be used to obtain the Fisher-Information matrix. If the observed Fisher information matrix
#' is estimate, nIter should be 1 since no new data has to be simulated. 
#' @param silent If TRUE, some diagnostic information is shown during estimation.
#' @param nSim A numeric value for the number of samples of the Gibbs sampler that is used internally.
#' @param n.rep A numeric value for the numer of MC estimates the initial estimate should be based on. If std.threshold is used
#' further estimates are computed until the criteria are satisfied.
#' @param n.cores The number of cores to use for the estimation. n.cores estimates are compued in parallel. 
#' @param nBurnin The number of samples to discard as burnin in the Gibbs sampler.
#' 
#' @return The ngme result object with the Fisher information added.
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'   ngme(...)
#'   }

ngme.fisher <- function(fit,
                        std.threshold = NULL,
                        observed = TRUE,
                        silent = FALSE,
                        n.cores = 1,
                        nSim = 1000,
                        n.rep = 10,
                        nIter = NULL,
                        only.effects = TRUE,
                        nBurnin = 100)
{
  if(missing(fit) || is.null(fit)){
    stop('No result object supplied.')
  }
  if(class(fit) != "ngme"){
    stop('The result object is not of class ngme.')
  }
  
  if(observed){
    type = 2
  } else {
    type = 1
  }
  
  
  
  if(observed && is.null(nIter)){
    nIter = 1
  } else if(observed == FALSE && is.null(nIter)){
    nIter = nSim
  }
  
  if(only.effects){
    fn <- names(fit$fixed_est)
  }
  
  if(silent == FALSE){
    cat("Compute initial estimate.\n")
  }
  cl <- makeCluster(n.cores, outfile= "")
  registerDoSNOW(cl)
  f.list <- foreach(i = 1:n.rep) %dopar%
  {
    out <- list()
    if(is.null(fit$processes_list)){
      fit.f <- estimateLong(fit$Y,
                            fit$locs,
                            fit$mixedEffect_list,
                            fit$measurementError_list,
                            nIter = nIter,
                            silent = 1,
                            nBurnin = 0,
                            nSim = nSim,
                            nBurnin_base = nBurnin,
                            estimate_fisher = type)
    } else {
      fit.f <- estimateLong(fit$Y,
                            fit$locs,
                            fit$mixedEffect_list,
                            fit$measurementError_list,
                            fit$processes_list,
                            fit$operator_list,
                            nIter = nIter,
                            silent = 1,
                            nBurnin = 0,
                            nSim = nSim,
                            nBurnin_base = nBurnin,
                            estimate_fisher = type)  
      out$X <- fit.f$processes_list$X
      out$V <- fit.f$processes_list$V
    }
    out$U <- fit.f$mixedEffect_list$U
    out$Fmat <- fit.f$FisherMatrix
    return(out)
  }
  stopCluster(cl)
  
  #collect initial estimates
  for(i in 1:n.cores){
    if(i==1){
      #sigma2.elements <- f.list[[i]]$sigma2.elements
      F.estimate <- f.list[[i]]$Fmat
      if(only.effects){
        ind.effects <- union(grep("beta_",c(rownames(F.estimate))),
                             grep("mu_",c(rownames(F.estimate)))
                             )
      } else {
        ind.effects <- 1:dim(F.estimate)[1]
      }
      F.estimate <- F.estimate[ind.effects,ind.effects]
      sigma2.elements <- matrix(diag(solve(F.estimate)),dim(F.estimate)[1],1)
    } else {
      Fi <- f.list[[i]]$Fmat[ind.effects,ind.effects]
      s2i <- matrix(diag(solve(Fi)),dim(Fi)[1],1)
      sigma2.elements <- cbind(sigma2.elements,s2i)
      F.estimate <- F.estimate + Fi
    }
  }
  sigma2.est <- rowMeans(sigma2.elements)
  sigma2.est.var <- rowMeans((sigma2.elements-sigma2.est)^2)/dim(sigma2.elements)[2]
  
  
  #update process samples in fit
  if(!is.null(fit$processes_list)){
    fit$processes_list$X <- f.list[[1]]$X
    fit$processes_list$V <- f.list[[1]]$V
  }
  fit$mixedEffect_list$U <- f.list[[1]]$U
  
  cat("\n")
  if(!is.null(std.threshold)){
    while(min(sigma2.est/sqrt(sigma2.est.var)) < std.threshold){
      ratios <- sigma2.est/sqrt(sigma2.est.var)
      ind <- which(ratios<std.threshold)
      if(silent == FALSE){
        rn <- rownames(F.estimate)[ind.effects]
        cat("Non-converged parameters:\n")
        df <- data.frame(est = sigma2.est[ind],std = sqrt(sigma2.est.var[ind]),ratio = ratios[ind],row.names = rn[ind])
        print(t(df))
        
      }
      cl <- makeCluster(n.cores)
      registerDoSNOW(cl)
      f.list <- foreach(i = 1:n.cores) %dopar%
      {
        out <- list()
        if(is.null(fit$processes_list)){
          fit.f <- estimateLong(fit$Y,
                                fit$locs,
                                fit$mixedEffect_list,
                                fit$measurementError_list,
                                fit$operator_list,
                                nIter = nIter,
                                silent = 1,
                                nBurnin = 0,
                                nSim = nSim,
                                nBurnin_base = nBurnin,
                                estimate_fisher = type)
        } else {
          fit.f <- estimateLong(fit$Y,
                                fit$locs,
                                fit$mixedEffect_list,
                                fit$measurementError_list,
                                fit$processes_list,
                                fit$operator_list,
                                nIter = nIter,
                                silent = 1,
                                nBurnin = 0,
                                nSim = nSim,
                                nBurnin_base = nBurnin,
                                estimate_fisher = type)  
          out$X <- fit.f$processes_list$X
          out$V <- fit.f$processes_list$V
        }
        out$U <- fit.f$mixedEffect_list$U
        out$Fmat <- fit.f$FisherMatrix
        return(out)
      }
      stopCluster(cl)
      #collect initial estimates
      for(i in 1:n.cores){
        Fi <- f.list[[i]]$Fmat[ind.effects,ind.effects]
        s2i <- matrix(diag(solve(Fi)),dim(Fi)[1],1)
        sigma2.elements <- cbind(sigma2.elements,s2i)
        F.estimate <- F.estimate + Fi
      }
      sigma2.est <- rowMeans(sigma2.elements)
      sigma2.est.var <- rowMeans((sigma2.elements-sigma2.est)^2)/dim(sigma2.elements)[2]
      
      #update process samples in fit
      if(!is.null(fit$processes_list)){
        fit$processes_list$X <- f.list[[1]]$X
        fit$processes_list$V <- f.list[[1]]$V
      }
      fit$mixedEffect_list$U <- f.list[[1]]$U
    }
  }
  
  
  
  if(silent == FALSE){
    cat("All estimates have converged. Saving results. \n")
  }
  
  # names for fixed effects in Fisher Matrix
  fisher_est <- F.estimate/dim(sigma2.elements)[2]
  
  #names <- c(names(fit$fixed_est), colnames(fit$ranef_Sigma))
  #colnames(fisher_est)[1:length(names)] <- rownames(fisher_est)[1:length(names)] <- names
  
  fit$fisher_est <- fisher_est
  fit$estimate_fisher = type
  fit$Fisher.inv.diag.var = sigma2.est.var
  names(fit$Fisher.inv.diag.var) <- colnames(fit$fisher_est)
  
  return(fit)
}
