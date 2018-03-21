
#' @title Estimation of Fisher information matrix
#'
#' @description Estimates Fisher information matrix for ngme result.
#'
#' @param nIter A numeric value for the number of iterations to be used to
#'       obtain the Fisher-Information matrix.
#' @param observed If TRUE, the observed Fisher information matrix is estimated. Otherwise the ordinary Fisher-Information matrix is estimated.
#' @param controls A list of control variables.
#'  \itemize{
#'     \item \code{nSim} A numeric value for the number of samples of the Gibbs sampler
#'       to estimate the gradient.
#'     \item \code{pSubsample} A numeric value for the portion of data to be used in each
#'       gradient iteration. \code{pSubsample = 1} indicates use of all subjects' data.
#'     \item \code{nSim.fisher} A numeric value for the number of samples of the Gibbs sampler
#'       to obtain the Fisher-Information matrix.
#'     \item \code{nBurnin.base} A numerical value for burn-in simulations that are performed
#'       for a subject that is sampled for the first time in the estimation method.
#'     \item \code{subsample.type} A numeric value for the type of subsampling;
#'       1: uniform sampling,
#'       2: sample size weighted,
#'       3: weighted sampling by gradient size,
#'       4: grouped sub-sampler.
#'     \item \code{pSubsample2} A numeric value for the portion of the data
#'       to be used in each gradient subsampling weighted by gradient.

#'  }
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

  cl <- makeCluster(n.cores)
  registerDoSNOW(cl)

  if(silent == FALSE){
    cat("Compute initial estimate.\n")
  }
  cl <- makeCluster(n.cores)
  registerDoSNOW(cl)
  f.list <- foreach(i = 1:n.rep) %dopar%
  {
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
    #update process samples in fit
    out <- list()
    out$X <- fit.f$processes_list$X
    out$V <- fit.f$processes_list$V
    out$U <- fit.f$mixedEffect_list$U

    out$Fmat <- fit.f$FisherMatrix
    out$sigma2.elements <- matrix(diag(solve(out$Fmat)),dim(out$Fmat)[1],1)
    return(out)
  }
  stopCluster(cl)

  #collect initial estimates
  for(i in 1:n.cores){
    if(i==1){
      sigma2.elements <- f.list[[i]]$sigma2.elements
      F.estimate <- f.list[[i]]$Fmat
    } else {
      sigma2.elements <- cbind(sigma2.elements,f.list[[i]]$sigma2.elements)
      F.estimate <- F.estimate + f.list[[i]]$Fmat
    }
  }
  sigma2.est <- rowMeans(sigma2.elements)
  sigma2.est.var <- rowMeans((sigma2.elements-sigma2.est)^2)/dim(sigma2.elements)[2]

  if(only.effects){
    ind.effects <- grep("beta_",c(rownames(F.estimate)))
  } else {
    ind.effects <- 1:dim(F.estimate)[1]
  }
  #update process samples in fit
  fit$processes_list$X <- f.list[[1]]$X
  fit$processes_list$V <- f.list[[1]]$V
  fit$mixedEffect_list$U <- f.list[[1]]$U

  cat("\n")
  if(!is.null(std.threshold)){
    while(min(sigma2.est[ind.effects]/sigma2.est.var[ind.effects]) < 1/std.threshold^2){
      ratios <- sigma2.est[ind.effects]/sigma2.est.var[ind.effects]
      ind <- which(ratios<1/std.threshold^2)
      if(silent == FALSE){
        rn <- rownames(F.estimate)[ind.effects]
        cat("Non-converged parameters:\n")
        df <- data.frame(est = sigma2.est[ind],var = sigma2.est.var[ind],ratio = ratios[ind],row.names = rn[ind])
        print(t(df))

      }
      cl <- makeCluster(n.cores)
      registerDoSNOW(cl)
      f.list <- foreach(i = 1:n.cores) %dopar%
      {
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
        #update process samples in fit
        out <- list()
        out$X <- fit.f$processes_list$X
        out$V <- fit.f$processes_list$V
        out$U <- fit.f$mixedEffect_list$U

        out$Fmat <- fit.f$FisherMatrix
        out$sigma2.elements <- matrix(diag(solve(out$Fmat)),dim(out$Fmat)[1],1)
        return(out)
      }
      stopCluster(cl)
      #collect initial estimates
      for(i in 1:n.cores){
        sigma2.elements <- cbind(sigma2.elements,f.list[[i]]$sigma2.elements)
        F.estimate <- F.estimate + f.list[[i]]$Fmat
      }
      sigma2.est <- rowMeans(sigma2.elements)
      sigma2.est.var <- rowMeans((sigma2.elements-sigma2.est)^2)/dim(sigma2.elements)[2]

      #update process samples in fit
      fit$processes_list$X <- f.list[[1]]$X
      fit$processes_list$V <- f.list[[1]]$V
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
