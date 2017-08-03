
#' @title Parameter estimation.
#'
#' @description Estimates model parameters using maximum likelihood that is 
#'   implemented by a computationally efficient stochastic gradient algorithm.
#'
#' @param fixed A two-sided formula to specify the fixed effects design matrix.
#' @param random A one-sided formula to specify the random effects design matrix.
#' @param use.process A logical variable for inclusion of the stochastic process in
#'   the mixed model: \code{"TRUE"} indicates inclusion, \code{"FALSE"} exclusion.
#' @param reffects A character string that indicates the distribution of the
#'   random effects. Available options are:  
#'   \code{"Normal"} for Normal, and 
#'   \code{"NIG"} for Normal-inverve Gaussian distributions.
#' @param process A character vector with two elements to specify 
#'   the process. Whilst the first element is for the covariance structure, the
#'   second element is for the process distribution. Available options for the
#'   first are: 
#'   \code{"fd2"} for integrated Random-Walk (integrated Brownian Motion), 
#'   \code{"matern"} for Matern family; 
#'   for the second element are:
#'   \code{"Normal"} for Normal, 
#'   \code{"NIG"} for Normal-inverse Gaussian,
#'   \code{"GAL"} for generalised-asymmetric Laplace, and 
#'   \code{"CH"} for Cauchy
#'   distributions. 
#'   \code{prcess} is ignored when \code{use.process} is set to FALSE.
#' @param error A character string to specify the distribution of the error term.
#'   Available options are: 
#'   \code{"Normal"} for Normal, 
#'   \code{"NIG"} for Normal-inverse Gaussian, 
#'   \code{"tdist"} for t-distribution.
#' @param data A data-frame from which the response and covariates to be extracted.
#' @param timevar A character string that indicates the name of the time variable. 
#'   It does not need to be specified when \code{use.process} is set to FALSE.
#' @param silent A logical value for printing the details of the iterations;
#'   \code{"TRUE"} indicates do not print, \code{"FALSE"} indicates print.  
#' @param nIter A numeric value for the number of iteration that will be
#'   used by the stochastic gradient.
#' @param controls A list of control variables.
#'  \itemize{
#'     \item \code{learning.rate} A numeric value for the parameter of stochastic gradient.
#'     \item \code{polyak.rate} A numeric value for moving average of parameters;
#'       -1: inactive, 0: pure mean.
#'     \item \code{nBurnin} A numeric value for the number of steps before starting
#'       gradient estimation.
#'     \item \code{nSim} A numeric value for the number of samples of the Gibbs sampler
#'       to estimate the gradient.
#'     \item \code{nIter.gauss} A numerical value for the number of iterations to be used
#'       to obtain the initial values for models with non-Gaussian processes.
#'     \item \code{pSubsample} A numeric value for the portion of data to be used in each
#'       gradient iteration. \code{pSubsample = 1} indicates use of all subjects' data.  
#'     \item \code{nPar.burnin} A numeric value; "M-step" updates will be used until this
#'       iteration.
#'     \item \code{nIter.fisher} A numeric value for the number of iterations to be used to
#'       obtain the Fisher-Information matrix.
#'     \item \code{nSim.fisher} A numeric value for the number of samples of the Gibbs sampler
#'       to obtain the Fisher-Information matrix.
#'     \item \code{step0} A numeric value for stepsize for the optimizer; step0 / i^alpha.
#'     \item \code{alpha} A numeric value for stepsize for the optimizer; step0 / i^alpha.
#'     \item \code{nBurnin.learningrate} A numeric value until which the learning will
#'       not be started.
#'     \item \code{nBurnin.base} A numerical value for burn-in simulations that are performed  
#'       for a subject that is sampled for the first time in the estimation method. 
#'     \item \code{subsample.type} A numeric value for the type of subsampling;
#'       1: uniform sampling, 
#'       2: sample size weighted,
#'       3: weighted sampling by gradient size, 
#'       4: grouped sub-sampler.
#'     \item \code{pSubsample2} A numeric value for the portion of the data
#'       to be used in each gradient subsampling weighted by gradient.
#'     \item \code{seed} A numerical value for starting the Gibbs samplers from fixed seed.
#'     \item \code{standardize.mixedEffects} A logical variable for standardising the covariates;
#'       \code{"FALSE"} indicates no standardisation, \code{"TRUE"} standardisation.
#'     \item \code{estimate.fisher} A logical variable for whether Fisher-Information matrix
#'       to be obtained; \code{"FALSE"} indicates do not obtain, \code{"TRUE"} obtain.
#'     \item \code{individual.sigma} A logical variable for specifying patient-specific mixture
#'       random variable for the error-term; \code{"FALSE"} indicates do not obtain,
#'       \code{"TRUE"} obtain.
#'     \item \code{n.process} A numerical value for the number of basis functions to
#'       approximate the stochastic process.
#'  }
#' @details This function is a user-friendly wrapper that calls the \code{estimateLong} function. 
#'     Generic functions \code{summary}, \code{print} and \code{plot} are available for the 
#'     output returned by the function \code{ngme}. For Matern covariance function, 
#'     currently the shape parameter is set to 0.5 which corresponds to exponential correlation 
#'     function. 
#' @return A list of outputs. 
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'   ngme(...)
#'   }

ngme <- function(fixed,
                 random,
                 use.process = TRUE,
                 reffects = "Normal",
                 process = c("Normal", "fd2"),
                 error = "Normal",
                 data,
                 timeVar = NULL,
                 silent = TRUE,
                 nIter,
                 controls = list(learning.rate = 0,
                                 polyak.rate = 0.1,
                                 nBurnin = 100,
                                 nSim = 2,
                                 nIter.gauss = 1000,
                                 pSubsample = 0.1,
                                 nPar.burnin = 0,
                                 nIter.fisher = 1000,
                                 nSim.fisher = 1000,
                                 step0 = 0.3,
                                 alpha = 0.3,
                                 nBurnin.learningrate = NULL,
                                 nBurnin.base = 0,
                                 subsample.type = 4,
                                 pSubsample2 = 0.3,
                                 seed = NULL,
                                 standardize.mixedEffects = FALSE,
                                 estimate.fisher = TRUE,
                                 individual.sigma = FALSE,
                                 n.process = NULL)
                 )
{
  
  # being sure that controls includes everything
  if(length(controls) < 20){
    controls.full <- list(learning.rate = 0,
                          polyak.rate = 0.1,
                          nBurnin = 100,
                          nSim = 2,
                          nIter.gauss = 1000,
                          pSubsample = 0.1,
                          nPar.burnin = 0,
                          nIter.fisher = 1000,
                          nSim.fisher = 1000, 
                          step0 = 0.3,
                          alpha = 0.3,
                          nBurnin.learningrate = NULL,
                          nBurnin.base = 0,
                          subsample.type = 4,
                          pSubsample2 = 0.3,
                          seed = NULL,
                          standardize.mixedEffects = FALSE,
                          estimate.fisher = TRUE,
                          individual.sigma = FALSE,
                          n.process = NULL)
    for(i in 1:length(controls.full)){
      if(!(names(controls.full)[i] %in% names(controls))){
        controls[names(controls.full)[i]] <- controls.full[i]
      }
    }
    
  }
  
  # correct input for distributions
  if(!(process[1] %in% c("NIG", "Normal", "GAL", "CH"))){
    stop("Process distribution should be one of the following: NIG, Normal, GAL, CH")
  }
  if(!(reffects %in% c("NIG", "Normal"))){
    stop("Random-effects distribution should be one of the following: NIG, Normal")
  }
  if(!(error %in% c("NIG", "Normal", "tdist"))){
    stop("Measurement error distribution should be one of the following: NIG, Normal, tdist")
  }
  
  # correct input for timeVar
  if(use.process == TRUE & is.null(timeVar) == TRUE){
    stop("timeVar should be specified, since the model consists of process")
  }
    
  # extract id variable
  idname <- rev(unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE)))[1]
  id <- data[, idname]
  
  # response matrix and fixed effects design matrix
  mf_fixed <- model.frame(formula = fixed, data = data)
  y        <- as.matrix(model.extract(mf_fixed, "response"))
  x_fixed_f  <- as.matrix(model.matrix(attr(mf_fixed, "terms"), data = mf_fixed))
  colnames(x_fixed_f)[1] <- gsub("[[:punct:]]", "", colnames(x_fixed_f)[1])
  
  # excluding the intercept and the covariates that are specified in random
  cov_list_fixed  <- attr(terms(fixed), "term.labels")
  cov_list_random <- unlist(strsplit(attr(terms(random), "term.labels"), " | ", fixed = TRUE))
  cov_list_random <- c(strsplit(cov_list_random[1], " + ", fixed=TRUE)[[1]], cov_list_random[2])
  cov_list_random <- cov_list_random[-length(cov_list_random)]
  
  to_del_x_fixed <- c("Intercept", cov_list_fixed[(cov_list_fixed %in% cov_list_random)])
  x_fixed <- x_fixed_f[, !(colnames(x_fixed_f) %in% to_del_x_fixed)]
  
  #random effects design matrix
  random_names             <- unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE))
  random_names_id_excluded <- random_names[!(random_names %in% idname)]
  random_formula           <- as.formula(paste("~", paste(random_names_id_excluded, collapse = "+")))
  
  mf_random <- model.frame(formula = random_formula, data = data)
  x_random  <- as.matrix(model.matrix(attr(mf_random, "terms"), data = mf_random))
  colnames(x_random)[1] <- gsub("[[:punct:]]", "", colnames(x_random)[1])
  
  idlist <- unique(id)
  
  # converting the followings to lists:
  # fixed effects design matrix, random effects design matrix, response matrix, time variable
  data_fixed <- data.frame(cbind(id, x_fixed))
  B_fixed    <- split(data_fixed[, -1], data_fixed[,1])
  B_fixed    <- lapply(B_fixed, function(x) as.matrix(x))
  
  data_random <- data.frame(cbind(id, x_random))
  B_random    <- split(data_random[, -1], data_random[,1])
  B_random    <- lapply(B_random, function(x) as.matrix(x))
  
  Y    <- tapply(y, id, function(x) x)
  locs <- tapply(data[, timeVar], id, function(x) x)
  
  ## Obtain starting values

  # setup the lists that are going to be passed into the 
  # functions that will obtain the starting value
  if(!silent){
    cat("Setup lists\n")
  }

  Vin <- lapply(Y, function(x) rep(1, length(x)))
  measurement_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1)
  
  mixedEffect_list <- list(B_random = B_random,
                           B_fixed  = B_fixed,
                           noise = "Normal",
                           Sigma_epsilon = 1
                           )
  if(use.process){
    if(is.null(controls$n.process)){
      n.process <- max(round(mean(unlist(lapply(locs, length)))), 1)
      operator_list <- create_operator(locs, n.process, name = process[2])
    }else{
      operator_list <- create_operator(locs, controls$n.process, name = process[2])
    }
    
    process_list = list(noise = "Normal",
                        nu  = 1,
                        mu  = 0
                        )
    process_list$V <- list()
    process_list$X <- list()
    
    for(i in 1:length(locs))
    {
      process_list$X[[i]] <- rep(0, length(operator_list$h[[1]]))
      process_list$V[[i]] <- operator_list$h[[1]]
    }
  }
  
  # starting values for measurement error and mixed effects using OLS:
  if(!silent){
    cat("Calculate starting values\n")    
  }

  mixedEffect_list       <- ME.startvalues(Y, mixedEffect_list)
  measurement_list$sigma <- mixedEffect_list$sigma
  
  # the case where the model formulation consists of W(t)
  if(use.process){
    
    #starting values for process
    operator_list$type  <- process[2]
    operator_list <- operator.startvalues(Y, 
                                          locs,
                                          mixedEffect_list, 
                                          operator_list, 
                                          measurement_list
                                          )
    # at least one of random effects, process or measurement error non-Gaussian   
    if(reffects != "Normal" || process[1] != "Normal" || error != "Normal"){
      
      # first fit the Gaussian model to obtain initials
      if(!silent){
        cat("Estimate Gaussian")
      }

      fit <- estimateLong(Y, 
                          locs,
                          mixedEffect_list,
                          measurement_list,
                          process_list,
                          operator_list,
                          nIter = nIter,
                          silent = silent,
                          learning_rate = controls$learning.rate,
                          polyak_rate = controls$polyak.rate,
                          nBurnin = controls$nBurnin,
                          nSim = controls$nSim,
                          pSubsample = controls$pSubsample,
                          nPar_burnin = controls$nPar.burnin,
                          step0 = controls$step0,
                          alpha = controls$alpha,
                          nBurnin_learningrate = controls$nBurnin.learningrate,
                          nBurnin_base = controls$nBurnin.base, 
                          subsample.type = controls$subsample.type,
                          pSubsample2 = controls$pSubsample2,
                          seed = controls$seed, 
                          standardize.mixedEffects = controls$standardize.mixedEffects,
                          estimate_fisher = FALSE
                          )
      
      # then fit the non-Gaussin model
      if(!silent){
        cat("Estimate non-Gaussian")
      }
      
      fit$mixedEffect_list$noise <- reffects
      fit$mixedEffect_list$nu    <- as.matrix(10)
      fit$mixedEffect_list$mu    <- matrix(0, dim(B_random[[1]])[2], 1)
      
      fit$processes_list$noise <- process[1]
      fit$processes_list$mu    <- 0
      fit$processes_list$nu    <- 10
      
      fit$measurementError_list$noise    <- error
      fit$measurementError_list$nu       <-  10
      fit$measurementError_list$common_V <- controls$individual.sigma
      fit$measurementError_list$Vs       <- Vin
      
      # Obtain parameter estimates
      fit <- estimateLong(Y, 
                          locs,
                          fit$mixedEffect_list,
                          fit$measurementError_list,
                          fit$processes_list,
                          fit$operator_list,
                          nIter = nIter,
                          silent = silent,
                          learning_rate = controls$learning.rate,
                          polyak_rate = controls$polyak.rate,
                          nBurnin = controls$nBurnin,
                          nSim = controls$nSim,
                          pSubsample = controls$pSubsample,
                          nPar_burnin = controls$nPar.burnin,
                          step0 = controls$step0,
                          alpha = controls$alpha,
                          nBurnin_learningrate = controls$nBurnin.learningrate,
                          nBurnin_base = controls$nBurnin.base, 
                          subsample.type = controls$subsample.type,
                          pSubsample2 = controls$pSubsample2,
                          seed = controls$seed, 
                          standardize.mixedEffects = controls$standardize.mixedEffects,
                          estimate_fisher = FALSE
                          )
      # Obtain Fisher matrix
      if(controls$estimate.fisher){
        fit.f <- estimateLong(Y, 
                              locs,
                              fit$mixedEffect_list,
                              fit$measurementError_list,
                              fit$processes_list,
                              fit$operator_list,
                              nIter = controls$nIter.fisher,
                              silent = silent,
                              learning_rate = controls$learning.rate,
                              polyak_rate = -1,
                              nBurnin = controls$nBurnin,
                              nSim = controls$nSim.fisher,
                              pSubsample = controls$pSubsample,
                              nPar_burnin = controls$nPar.burnin,
                              step0 = controls$step0,
                              alpha = controls$alpha,
                              nBurnin_learningrate = controls$nBurnin.learningrate,
                              nBurnin_base = controls$nBurnin.base, 
                              subsample.type = controls$subsample.type,
                              pSubsample2 = controls$pSubsample2,
                              seed = controls$seed, 
                              standardize.mixedEffects = controls$standardize.mixedEffects,
                              estimate_fisher = TRUE
                              )
        fit$FisherMatrix <- fit.f$FisherMatrix
      }
    } else { # random effects, process and measurement error are all Gaussian
      
      if(!silent){
        cat("Estimate Model")
      }
      
      # Obtain parameter estimates
      fit <- estimateLong(Y, 
                          locs,
                          mixedEffect_list,
                          measurement_list,
                          process_list,
                          operator_list,
                          nIter = nIter,
                          silent = silent,
                          learning_rate = controls$learning.rate,
                          polyak_rate = controls$polyak.rate,
                          nBurnin = controls$nBurnin,
                          nSim = controls$nSim,
                          pSubsample = controls$pSubsample,
                          nPar_burnin = controls$nPar.burnin,
                          step0 = controls$step0,
                          alpha = controls$alpha,
                          nBurnin_learningrate = controls$nBurnin.learningrate,
                          nBurnin_base = controls$nBurnin.base, 
                          subsample.type = controls$subsample.type,
                          pSubsample2 = controls$pSubsample2,
                          seed = controls$seed, 
                          standardize.mixedEffects = controls$standardize.mixedEffects,
                          estimate_fisher = FALSE
                          )
      # Obtain Fisher matrix
      if(controls$estimate.fisher){
        fit.f <- estimateLong(Y, 
                              locs,
                              fit$mixedEffect_list,
                              fit$measurementError_list,
                              fit$processes_list,
                              fit$operator_list,
                              nIter = controls$nIter.fisher,
                              silent = silent,
                              learning_rate = controls$learning.rate,
                              polyak_rate = -1,
                              nBurnin = controls$nBurnin,
                              nSim = controls$nSim.fisher,
                              pSubsample = controls$pSubsample,
                              nPar_burnin = controls$nPar.burnin,
                              step0 = controls$step0,
                              alpha = controls$alpha,
                              nBurnin_learningrate = controls$nBurnin.learningrate,
                              nBurnin_base = controls$nBurnin.base, 
                              subsample.type = controls$subsample.type,
                              pSubsample2 = controls$pSubsample2,
                              seed = controls$seed, 
                              standardize.mixedEffects = controls$standardize.mixedEffects,
                              estimate_fisher = TRUE
                              )
        fit$FisherMatrix <- fit.f$FisherMatrix
      }
    }
  } else {## the case of the model excludes W(t)
    
    # either random effects or measurement error is non-Gaussian
    if(random.effect.distribution != "Normal" || measurement.distribution != "Normal"){
      
      if(!silent){
        cat("Estimate Gaussian") 
      }
      
      # first fit the Gaussian model to obtain the initials 
      fit <- estimateLong(Y, 
                          locs,
                          mixedEffect_list,
                          measurement_list,
                          nIter = controls$nIter.gauss,
                          silent = silent,
                          learning_rate = controls$learning.rate,
                          polyak_rate = controls$polyak.rate,
                          nBurnin = controls$nBurnin,
                          nSim = controls$nSim,
                          pSubsample = controls$pSubsample,
                          nPar_burnin = controls$nPar.burnin,
                          step0 = controls$step0,
                          alpha = controls$alpha,
                          nBurnin_learningrate = controls$nBurnin.learningrate,
                          nBurnin_base = controls$nBurnin.base, 
                          subsample.type = controls$subsample.type,
                          pSubsample2 = controls$pSubsample2,
                          seed = controls$seed, 
                          standardize.mixedEffects = controls$standardize.mixedEffects,
                          estimate_fisher = FALSE
                          )
      
      # then fit the non-Gaussian model
      if(!silent){
        cat("Estimate non-Gaussian")
      }
      
      fit$mixedEffect_list$noise <- reffects
      fit$mixedEffect_list$nu    <- as.matrix(10)
      fit$mixedEffect_list$mu    <- matrix(0, dim(B_random[[1]])[2], 1)
      
      fit$measurementError_list$noise    <- error
      fit$measurementError_list$nu       <- 10
      fit$measurementError_list$common_V <- controls$individual.sigma
      fit$measurementError_list$Vs       <- Vin
      
      # Obtain parameter estimates
      fit <- estimateLong(Y, 
                          locs,
                          fit$mixedEffect_list,
                          fit$measurementError_list,
                          fit$processes_list,
                          fit$operator_list,
                          nIter = nIter,
                          silent = silent,
                          learning_rate = controls$learning.rate,
                          polyak_rate = controls$polyak.rate,
                          nBurnin = controls$nBurnin,
                          nSim = controls$nSim,
                          pSubsample = controls$pSubsample,
                          nPar_burnin = controls$nPar.burnin,
                          step0 = controls$step0,
                          alpha = controls$alpha,
                          nBurnin_learningrate = controls$nBurnin.learningrate,
                          nBurnin_base = controls$nBurnin.base, 
                          subsample.type = controls$subsample.type,
                          pSubsample2 = controls$pSubsample2,
                          seed = controls$seed, 
                          standardize.mixedEffects = controls$standardize.mixedEffects,
                          estimate_fisher = FALSE
                          )
      # Obtain Fisher matrix
      if(controls$estimate.fisher){
        fit.f <- estimateLong(Y, 
                              locs,
                              fit$mixedEffect_list,
                              fit$measurementError_list,
                              fit$processes_list,
                              fit$operator_list,
                              nIter = nIter,
                              silent = silent,
                              learning_rate = controls$learning.rate,
                              polyak_rate = -1,
                              nBurnin = controls$nBurnin,
                              nSim = controls$nSim,
                              pSubsample = controls$pSubsample,
                              nPar_burnin = controls$nPar.burnin,
                              step0 = controls$step0,
                              alpha = controls$alpha,
                              nBurnin_learningrate = controls$nBurnin.learningrate,
                              nBurnin_base = controls$nBurnin.base, 
                              subsample.type = controls$subsample.type,
                              pSubsample2 = controls$pSubsample2,
                              seed = controls$seed, 
                              standardize.mixedEffects = controls$standardize.mixedEffects,
                              estimate_fisher = TRUE
                              )
        fit$FisherMatrix <- fit.f$FisherMatrix
      }
    } else {# both random effects and measurement error are Gaussian
      
      # Obtain parameter estimates
      fit <- estimateLong(Y, 
                          locs,
                          mixedEffect_list,
                          measurement_list,
                          nIter = nIter,
                          silent = silent,
                          learning_rate = controls$learning.rate,
                          polyak_rate = controls$polyak.rate,
                          nBurnin = controls$nBurnin,
                          nSim = controls$nSim,
                          pSubsample = controls$pSubsample,
                          nPar_burnin = controls$nPar.burnin,
                          step0 = controls$step0,
                          alpha = controls$alpha,
                          nBurnin_learningrate = controls$nBurnin.learningrate,
                          nBurnin_base = controls$nBurnin.base, 
                          subsample.type = controls$subsample.type,
                          pSubsample2 = controls$pSubsample2,
                          seed = controls$seed, 
                          standardize.mixedEffects = controls$standardize.mixedEffects,
                          estimate_fisher = FALSE
                          )
      # Obtain Fisher matrix
      if(controls$estimate.fisher){
        fit.f <- estimateLong(Y, 
                              locs,
                              res$mixedEffect_list,
                              res$measurementError_list,
                              learning_rate = controls$learning.rate,
                              nIter = nIter,
                              silent = silent,
                              learning_rate = controls$learning.rate,
                              polyak_rate = -1,
                              nBurnin = controls$nBurnin,
                              nSim = controls$nSim,
                              pSubsample = controls$pSubsample,
                              nPar_burnin = controls$nPar.burnin,
                              step0 = controls$step0,
                              alpha = controls$alpha,
                              nBurnin_learningrate = controls$nBurnin.learningrate,
                              nBurnin_base = controls$nBurnin.base, 
                              subsample.type = controls$subsample.type,
                              pSubsample2 = controls$pSubsample2,
                              seed = controls$seed, 
                              standardize.mixedEffects = controls$standardize.mixedEffects,
                              estimate_fisher = TRUE
                              )
        fit$FisherMatrix <- fit.f$FisherMatrix
      }
    }
  }
  
  
  #
  # Preparing the output
  #
  
  pSubsample_fit <- fit$pSubsample
  nIter_fit      <- fit$nIter
  nSim_fit       <- fit$nSim
  nBurnin_fit    <- fit$nBurnin
  step0_fit      <- fit$step0
  alpha_fit      <- fit$alpha
  
  # fixed effects estimates - and chains
  fixed_est1 <- as.numeric(fit$mixedEffect_list$beta_fixed)
  fixed_est2 <- as.numeric(fit$mixedEffect_list$beta_random)
  names(fixed_est1) <- colnames(x_fixed)
  names(fixed_est2) <- to_del_x_fixed
  fixed_est <- rep(NA, ncol(x_fixed_f))
  names(fixed_est) <- colnames(x_fixed_f)
  fixed_est[names(fixed_est) %in% names(fixed_est1)] <- fixed_est1
  fixed_est[names(fixed_est) %in% names(fixed_est2)] <- fixed_est2
  
  fixed_est1_vec <- fit$mixedEffect_list$betaf_vec
  fixed_est2_vec <- fit$mixedEffect_list$betar_vec
  colnames(fixed_est1_vec) <- colnames(x_fixed)
  colnames(fixed_est2_vec) <- to_del_x_fixed
  fixed_est_vec <- matrix(NA, ncol = ncol(x_fixed_f), nrow = nIter)
  colnames(fixed_est_vec) <- colnames(x_fixed_f)
  fixed_est_vec[, colnames(fixed_est_vec) %in% names(fixed_est1)] <- fixed_est1_vec
  fixed_est_vec[, colnames(fixed_est_vec) %in% names(fixed_est2)] <- fixed_est2_vec
  
  # random effects
  ranef_Sigma           <- fit$mixedEffect_list$Sigma
  colnames(ranef_Sigma) <- rownames(ranef_Sigma) <- colnames(x_random)
  ranef_Sigma_vec       <- fit$mixedEffect_list$Sigma_vec
  ranef_sigma_epsilon   <- fit$mixedEffect_list$Sigma_epsilon
  
  if(reffects %in% c("NIG")){
    ranef_mu     <- fit$mixedEffect_list$mu
    ranef_mu_vec <- fit$mixedEffect_list$mu_vec
    
    ranef_nu <- fit$mixedEffect_list$nu
    ranef_nu_vec <- fit$mixedEffect_list$nu_vec
  }else{
    ranef_mu <- ranef_mu_vec <- ranef_nu <- ranef_nu_vec <- NA
  }
  
  # process and operator
  if(use.process == TRUE){
    
    #operator
    operator_tau     <- fit$operator_list$tau
    operator_tau_vec <- fit$operator_list$tauVec[-1]
    
    #operator - matern
    if(process[2] %in% c("matern")){
      operator_kappa <- fit$operator_list$kappa
      operator_kappa_vec <- fit$operator_list$kappaVec
    }else{
      operator_kappa <- operator_kappa_vec <- NA
    }
    
    #process
    if(process[1] %in% c("NIG", "GAL")){
      process_nu <- fit$processes_list$nu
      process_nu_vec <- fit$processes_list$nu_vec
      
      process_mu <- fit$processes_list$mu
      process_mu_vec <- fit$processes_list$mu_vec
    }else{
      process_nu <- process_nu_vec <- process_mu <- process_mu_vec <- NA
    }
    
  }else{
    operator_tau <- operator_tau_vec <- operator_kappa <- operator_kappa_vec <-
      process_nu <- process_nu_vec <- process_mu <- process_mu_vec <- NA
  }
  
  # measurement error
  meas_error_sigma <- fit$measurementError_list$sigma
  meas_error_sigma_vec <- fit$measurementError_list$sigma_vec
  
  if(error %in% c("NIG", "tdist")){
    meas_error_nu <- fit$measurementError_list$nu
    meas_error_nu_vec <- fit$measurementError_list$nu_vec
  }else{
    meas_error_nu <- meas_error_nu_vec <- NA
  }
  
  # checking Fisher Matrix is estimated?
  if(controls$estimate.fisher == TRUE){
    
    # names for fixed effects in Fisher Matrix
    fisher_est <- fit$FisherMatrix
    colnames(fisher_est)[1:ncol(x_fixed_f)] <-
      rownames(fisher_est)[1:ncol(x_fixed_f)] <-
      c(colnames(x_fixed), colnames(x_random))
    
  }else{
    fisher_est <- NA
  }
  
  out <- list(
    use_process = use.process,
    estimate_fisher = controls$estimate.fisher,
    x_fixed_f = x_fixed_f,
    x_random = x_random,
    random_distr = reffects,
    process_distr = process[1],
    operator_type = process[2],
    error_distr = error,
    Y = Y,
    locs = locs,
    mixedEffect_list = fit$mixedEffect_list,
    measurementError_list = fit$measurementError_list,
    processes_list = fit$processes_list,
    operator_list = fit$operator_list,
    pSubsample_fit = pSubsample_fit,
    nIter_fit = nIter_fit,
    nSim_fit = nSim_fit,
    nBurnin_fit = nBurnin_fit,
    step0_fit = step0_fit,
    alpha_fit = alpha_fit,
    fixed_est = fixed_est,
    fixed_est_vec = fixed_est_vec,
    ranef_Sigma = ranef_Sigma,
    ranef_Sigma_vec = ranef_Sigma_vec,
    ranef_sigma_epsilon = ranef_sigma_epsilon,
    ranef_mu = ranef_mu,
    ranef_mu_vec = ranef_mu_vec,
    ranef_nu = ranef_nu,
    ranef_nu_vec = ranef_nu_vec,
    operator_tau = operator_tau,
    operator_tau_vec = operator_tau_vec,
    operator_kappa = operator_kappa,
    operator_kappa_vec = operator_kappa_vec,
    process_nu = process_nu,
    process_nu_vec = process_nu_vec,
    process_mu = process_mu,
    process_mu_vec = process_mu_vec,
    meas_error_sigma = meas_error_sigma,
    meas_error_sigma_vec = meas_error_sigma_vec,
    meas_error_nu = meas_error_nu,
    meas_error_nu_vec = meas_error_nu_vec,
    fisher_est = fisher_est,
    call = match.call()
  )
  
  class(out) <- "ngme"
  out
  
}
