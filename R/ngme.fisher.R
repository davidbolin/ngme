
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
                        std.threshold = 0.1,
                        observed = TRUE,
                        silent = FALSE,
                        n.cores = 1,
                        nSim = 1000,
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

  if(n.cores == 1){
    n.start = 4
    if(silent == FALSE){
      cat("initial estimate: ")
    }

    for(i in 1:n.start){
      if(silent == FALSE){
        cat(i," ")
      }
      if(i==1){
        fit.f <- estimateLong(fit$Y,
                              fit$locs,
                              fit$mixedEffect_list,
                              fit$measurementError_list,
                              fit$processes_list,
                              fit$operator_list,
                              nIter = 1,
                              silent = 1,
                              nBurnin = 0,
                              nSim = nSim,
                              nBurnin_learningrate = 0,
                              nBurnin_base = nBurnin,
                              pSubsample2 = 1,
                              estimate_fisher = type)
        #update process samples in fit
        fit$processes_list$X <- fit.f$processes_list$X
        fit$processes_list$V <- fit.f$processes_list$V
        fit$mixedEffect_list$U <- fit.f$mixedEffect_list$U
      } else {
        fit.f <- estimateLong(fit$Y,
                              fit$locs,
                              fit$mixedEffect_list,
                              fit$measurementError_list,
                              fit$processes_list,
                              fit$operator_list,
                              nIter = 1,
                              silent = 1,
                              nBurnin = 0,
                              nSim = nSim,
                              nBurnin_learningrate = 0,
                              nBurnin_base = nBurnin,
                              pSubsample2 = 1,
                              estimate_fisher = type)
      }

      Fmat <- fit.f$FisherMatrix
      Fi <- solve(Fmat)
      #cat("est: ", diag(Fi),"\n")
      if(i==1){
        F.estimate <- Fmat
        sigma2.elements <- matrix(diag(Fi),dim(Fi)[1],1)
      } else {
        F.estimate <- F.estimate + Fmat
        sigma2.elements <- cbind(sigma2.elements,matrix(diag(Fi),dim(Fi)[1],1))
      }

    }
    sigma2.est <- rowMeans(sigma2.elements)
    sigma2.est.var <- rowMeans((sigma2.elements-sigma2.est)^2)/dim(sigma2.elements)[2]

    i = n.start+1
    cat("\n")
    while(min(sigma2.est/sqrt(sigma2.est.var)) < 1/std.threshold){
      ratios <- sigma2.est/sqrt(sigma2.est.var)
      ind <- which(ratios<1/std.threshold)
      if(silent == FALSE){
        cat("Non-converged parameters:")
        cat("Estimates:\n")
        cat(sigma2.est[ind],"\n")
        cat("Stds:\n")
        cat(sqrt(sigma2.est.var[ind]),"\n")
        cat("Ratio:\n")
        cat(ratios[ind],"\n")
      }

      fit.f <- estimateLong(fit$Y,
                            fit$locs,
                            fit$mixedEffect_list,
                            fit$measurementError_list,
                            fit$processes_list,
                            fit$operator_list,
                            nIter = 1,
                            silent = 1,
                            nBurnin = 0,
                            nSim = nSim,
                            nBurnin_learningrate = 0,
                            nBurnin_base = nBurnin,
                            pSubsample2 = 1,
                            estimate_fisher = type)
      #update process samples in fit
      fit$processes_list$X <- fit.f$processes_list$X
      fit$processes_list$V <- fit.f$processes_list$V
      fit$mixedEffect_list$U <- fit.f$mixedEffect_list$U
      Fmat = fit.f$FisherMatrix
      Fi <- solve(Fmat)
      #cat("est: ", diag(Fi),"\n")
      F.estimate <- F.estimate + Fmat
      sigma2.elements <- cbind(sigma2.elements,matrix(diag(Fi),dim(Fi)[1],1))
      sigma2.est <- rowMeans(sigma2.elements)
      sigma2.est.var <- rowMeans((sigma2.elements-sigma2.est)^2)/dim(sigma2.elements)[2]
      i = i+1
    }

  } else { #Parallell computations
    cl <- makeCluster(n.cores)
    registerDoSNOW(cl)
    if(silent == FALSE){
      cat("Compute initial estimate.\n")
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
                            nIter = 1,
                            silent = 1,
                            nBurnin = 0,
                            nSim = nSim,
                            nBurnin_base = 10*nBurnin,# Large burnin first time
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

    #update process samples in fit
    fit$processes_list$X <- f.list[[1]]$X
    fit$processes_list$V <- f.list[[1]]$V
    fit$mixedEffect_list$U <- f.list[[1]]$U

    cat("\n")
    while(min(sigma2.est/sqrt(sigma2.est.var)) < 1/std.threshold){
      ratios <- sigma2.est/sqrt(sigma2.est.var)
      ind <- which(ratios<1/std.threshold)
      if(silent == FALSE){
        cat("Non-converged parameters:")
        cat("Estimates:\n")
        cat(sigma2.est[ind],"\n")
        cat("Stds:\n")
        cat(sqrt(sigma2.est.var[ind]),"\n")
        cat("Ratio:\n")
        cat(ratios[ind],"\n")
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
                              nIter = 1,
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


  # names for fixed effects in Fisher Matrix
  fisher_est <- F.estimate/dim(sigma2.elements)[2]

  names <- c(names(fit$fixed_est), colnames(fit$ranef_Sigma))
  colnames(fisher_est)[1:length(names)] <- rownames(fisher_est)[1:length(names)] <- names

  fit$fisher_est <- fisher_est
  fit$estimate_fisher = type

  return(fit)
}
