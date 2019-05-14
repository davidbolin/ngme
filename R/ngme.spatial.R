
#' @title Parameter estimation of non-Gaussian spatial models.
#'
#' @description Likelihood-based parameter estimation of spatial non-Gaussian models.
#'
#' @param fixed A two-sided formula to specify the fixed effects design matrix.
#' @param random A one-sided formula to specify the random effects design matrix (if any).
#' @param use.process A logical variable for inclusion of the stochastic process in
#'   the mixed model. Default is \code{"TRUE"}.
#' @param reffects A character string that indicates the distribution of the
#'   random effects if present in the model. Available options are:
#'   \code{"Normal"} for Normal distribution, and
#'   \code{"NIG"} for Normal-inverse Gaussian.
#' @param process A character string specifying the distribution of the driving noise of the 
#' spatial process  Available options are
#'   first are:
#'   \code{"Normal"} for Normal distribution,
#'   \code{"NIG"} for Normal-inverse Gaussian,
#'   \code{"GAL"} for generalised-asymmetric Laplace, and
#'   \code{"CH"} for Cauchy.
#' @param error A character string to specify the distribution of the error term.
#'   Available options are:
#'   \code{"Normal"} for Normal distribution,
#'   \code{"NIG"} for Normal-inverse Gaussian,
#'   \code{"tdist"} for t.
#' @param error_assymetric if true the non-Gaussian error is assymetric
#' @param data A data-frame from which the response and covariates to be extracted.
#' @param location.names A character vector with the names of the spatial coordinates.
#' @param silent A logical value for printing the details of the iterations;
#'   \code{"TRUE"} indicates do not print, \code{"FALSE"} print.
#' @param nIter A numeric value for the number of iteration that will be
#'   used by the stochastic gradient.
#' @param mesh A mesh object of class inla.mesh
#' @param controls A list of control variables for parameter estimation. See \code{ngme} for details.
#' @param controls.init A list of control variables to be used to fit the normal model
#'      to get the initial values for fitting a model with at least one of random effects,
#'      process and error being non-Gaussian. See \code{ngme} for details.
#' @param init.fit A fitted \code{ngme.spatial} object to be used as a starting value for estimation.
#' @details The model that is estimated is of the form
#' 
#' \deqn{Y_{ij} = B(s_i)beta + U(s_i)beta_j + X_j(s_i) + e_{ij}}
#' Here i denots the index of the spatial location for the measurement and j denotes the number of the 
#' replicate in case of repeated measurements. The common mean value \eqn{B(s_i)beta} is specified using 
#' fixed effects, and the random effect part \eqn{U(s_i)beta_j} specifies a mean value varying between replicates. 
#' The process \eqn{X_j(s)} is a spatial process (independent between replicates) with a Matern covariance 
#' structure, specified as a stochastic PDE with possibly non-Gaussian driving noise. 
#' The measurement noise \eqn{e_{ij}} is assumed to be independent between observations and replicates. 
#' @return A list of outputs.

ngme.spatial <- function(fixed,
                 random = NULL,
                 fixed2 = NULL,
                 random2 = NULL,
                 group.id = NULL,
                 use.process = TRUE,
                 reffects = "Normal",
                 process = c("Normal", "matern"),
                 error = "Normal",
                 error_assymetric = FALSE,
                 data,
                 location.names = NULL,
                 silent = TRUE,
                 nIter = 1000,
                 mesh = NULL,
                 controls = list(learning.rate = 0,
                                 polyak.rate = 0.1,
                                 nBurnin = 100,
                                 nSim = 2,
                                 pSubsample = NULL,
                                 nPar.burnin = 0,
                                 step0 = 0.3,
                                 alpha = 0.3,
                                 nBurnin.learningrate = NULL,
                                 nBurnin.base = 0,
                                 subsample.type = 1,
                                 pSubsample2 = 0.3,
                                 individual.sigma = FALSE),
                 controls.init = list(learning.rate.init = 0,
                                      polyak.rate.init = 0.1,
                                      nBurnin.init = 100,
                                      nSim.init = 2,
                                      nIter.init = 1000,
                                      pSubsample.init = 0.1,
                                      nPar.burnin.init = 0,
                                      step0.init = 0.3,
                                      alpha.init = 0.3,
                                      nBurnin.learningrate.init = NULL,
                                      nBurnin.base.init = 0,
                                      subsample.type.init = 0,
                                      pSubsample2.init = 0.3,
                                      individual.sigma.init = FALSE),
                 init.fit = NULL,
                 debug = FALSE
)
{
  # generate a seed
  gen_seed <- ceiling(10^8 * runif(1))
  controls$seed <- controls.init$seed.init <- gen_seed
  
  # being sure that controls includes everything
  if(length(controls) < 13){
    controls.full <- list(learning.rate = 0,
                          polyak.rate = 0.1,
                          nBurnin = 100,
                          nSim = 2,
                          pSubsample = NULL,
                          nPar.burnin = 0,
                          step0 = 0.3,
                          alpha = 0.3,
                          nBurnin.learningrate = NULL,
                          nBurnin.base = 0,
                          subsample.type = 4,
                          pSubsample2 = 0.3,
                          individual.sigma = FALSE
    )
    for(i in 1:length(controls.full)){
      if(!(names(controls.full)[i] %in% names(controls))){
        controls[names(controls.full)[i]] <- controls.full[i]
      }
    }
    
  }
  
  # check for controls.init
  if(is.null(init.fit) == TRUE ||
     (reffects == "Normal" & (use.process == TRUE & process[1] == "Normal") & error == "Normal") ||
     (reffects == "Normal" & use.process == FALSE & error == "Normal")
  ){
    
    if(length(controls.init) < 14){
      controls.init.full <- list(learning.rate.init = 0,
                                 polyak.rate.init = 0.1,
                                 nBurnin.init = 100,
                                 nSim.init = 2,
                                 nIter.init = 1000,
                                 pSubsample.init = NULL,
                                 nPar.burnin.init = 0,
                                 step0.init = 0.3,
                                 alpha.init = 0.3,
                                 nBurnin.learningrate.init = NULL,
                                 nBurnin.base.init = 0,
                                 subsample.type.init = 4,
                                 pSubsample2.init = 0.3,
                                 individual.sigma.init = FALSE)
      for(i in 1:length(controls.init.full)){
        if(!(names(controls.init.full)[i] %in% names(controls.init))){
          controls.init[names(controls.init.full)[i]] <- controls.init.full[i]
        }
      }
      
    }
    
  }
  
  ## check mesh
  if(class(mesh) != "inla.mesh" && use.process){
    stop("Provide 'mesh' if process is used")
  }
  
  # correct input for distributions
  if(!(process[1] %in% c("NIG", "Normal", "GAL", "CH"))){
    stop("Process distribution should be one of the following: 'NIG', 'Normal', 'GAL', 'CH'")
  }
  if(!(reffects %in% c("NIG", "Normal", "tdist"))){
    stop("Random-effects distribution should be one of the following: 'NIG', 'Normal', 'tdist'")
  }
  if(!(error %in% c("NIG", "Normal", "tdist"))){
    stop("Measurement error distribution should be one of the following: 'NIG', 'Normal', 'tdist'")
  }
  
  # correct input for timeVar
  if(use.process == TRUE & is.null(location.names) == TRUE){
    stop("'location.names' should be specified")
  }
  
  # alpha and alpha.init are in the correct interval?
  if(controls$alpha < 0 | controls$alpha > 1){
    stop("alpha should be in (0, 1]")
  }
  if(controls.init$alpha.init < 0 | controls.init$alpha.init > 1){
    stop("alpha.init should be in (0, 1]")
  }
  
  # extract id variable
  if(is.null(random)){
    use.random = FALSE
  } else {
    use.random = TRUE
  }
  
  if(use.random){
    idname <- rev(unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE)))[1]
    id <- data[, idname]  
  } else if(!is.null(group.id)){
    idname = group.id
    id <- data[, idname]  
  } else {
    idname = NULL
  }
  effects <- extract.effects(data = data, fixed = fixed,random=random, idname = idname)
  Y<-effects$Y
  B_fixed <- effects$B_fixed
  x_fixed <- effects$x_fixed
  x_random <- effects$x_random
  x_fixed_f <- effects$x_fixed_f
  to_del_x_fixed <- effects$to_del_x_fixed
  
  #if we have a bivariate model, extract effect matrices for second dimension. 
  bivariate = FALSE
  if(!is.null(fixed2)){ 
    bivariate = TRUE
    effects2 <- extract.effects(data = data, fixed = fixed2,random=random2, idname = idname)
    x_random <- cbind(x_random,effects2$x_random)
    to_del_x_fixed <- cbind(to_del_x_fixed,effects2$to_del_x_fixed)
    x_fixed <- cbind(x_fixed,effects2$x_fixed)
    x_fixed_f <- cbind(x_fixed_f,effects2$x_fixed_f)
  
    #combine fixed effect matrices and data
    for(i in 1:length(B_fixed)){
      B_fixed[[i]] <- as.matrix(bdiag(B_fixed[[i]],effects2$B_fixed[[i]]))
      Y[[i]] <- cbind(Y[[i]],effects2$Y[[i]])
      if(use.random){
        B_random[[i]] <- as.matrix(bdiag(B_random[[i]],effects2$B_random[[i]]))
      }
    }
  }
  
  # extract variables for process
  if(use.process == TRUE){
    if(is.null(idname)){
      locs <- list(as.matrix(data[, location.names]))  
    } else {
      locs <- split(data[, location.names], id, function(x) x) 
      locs <- lapply(locs, function(x) as.matrix(x))
    }
  }
  
  nsubj <- length(Y)
  
  # if pSubsampling not set
  if(is.null(controls.init$pSubsample.init)){
    if(nsubj < 100){
      controls.init$pSubsample.init = 1
    }else if(nsubj < 500){
      controls.init$pSubsample.init = 0.2
      warning("pSubsample.init not provided. Since there are >=100 and <500 replicates, p for subsampling is set to 0.2")
    }else{
      controls.init$pSubsample.init = 0.1
      warning("pSubsample.init not provided. Since there >=500 replicates, p for subsampling is set to 0.1")
    }
  }
  if(is.null(controls$pSubsample)){
    if(nsubj < 100){
      controls$pSubsample = 1
    }else if(nsubj < 500){
      controls$pSubsample = 0.2
      warning("pSubsample not provided. Since there are >=100 and <500 replicates, p for subsampling is set to 0.2")
    }else{
      controls$pSubsample = 0.1
      warning("pSubsample not provided. Since there >=500 replicates, p for subsampling is set to 0.1")
    }
  }
  
  ## Vin is needed even if init.fit is not NULL
  Vin <- lapply(Y, function(x) rep(1, length(x)))
  
  ## Obtain starting values - if init.fit is not supplied or everything is Gaussian
  
  if(is.null(init.fit) == TRUE ||
     (reffects == "Normal" & (use.process == TRUE & process[1] == "Normal") & error == "Normal") ||
     (reffects == "Normal" & use.process == FALSE & error == "Normal")){
    
    # setup the lists that are going to be passed into the
    # functions that will obtain the starting value
    if(!silent){
      cat("Setup lists\n")
    }
    
    if(bivariate){
      Be <- list()
      for(i in 1:length(Y)){
        Be[[i]] = kronecker(diag(2),matrix(rep(1,dim(Y[[i]])[1])))
      }
      measurement_list <- list(Vs = Vin,
                               noise = "nsNormal",
                               B = Be,
                               sigma = c(0.1,0.1))
    } else {
      measurement_list <- list(Vs = Vin,
                               noise = "Normal",
                               sigma = 0.1)  
    }
    
    
    mixedEffect_list <- list(B_fixed  = B_fixed,
                             noise = "Normal",
                             Sigma_epsilon = 1)
    if(use.random){
      mixedEffect_list$B_random = B_random
    }
    
    if(use.process){
      if(bivariate){
        operator_list <- create_operator_matern2Dbivariate(mesh)
        n.grid <- length(operator_list$h[[1]])/2
        
        process_list = list(noise = "Normal", 
                            mu = as.matrix(c(0,0)), 
                            nu = as.matrix(c(1,1)))
        process_list$V <- list()
        process_list$X <- list()
        process_list$Bmu <- list()
        process_list$Bnu <- list()
      } else {
        operator_list <- create_operator_matern2D(mesh)  
        process_list = list(noise = "Normal",
                            nu  = 1,
                            mu  = 0)
        process_list$V <- list()
        process_list$X <- list()
      }
      for(i in 1:length(locs))
      {
        if(length(operator_list$h)==1){
          h_in <- operator_list$h[[1]]
        }else{
          h_in <- operator_list$h[[i]]
        }
        process_list$X[[i]] <- rep(0, length(h_in))
        process_list$V[[i]] <- h_in
        if(bivariate){
          process_list$Bmu[[i]] = kronecker(diag(2),matrix(rep(1, n.grid)))
          process_list$Bnu[[i]] = kronecker(diag(2),matrix(rep(1, n.grid)))
        }
      }
    }
     
    # starting values for measurement error and mixed effects using OLS:
    if(!silent){
      cat("Calculate starting values\n")
    }
    if(bivariate){
      mixedEffect_list       <- ME.startvalues.bivariate(Y, mixedEffect_list)
      measurement_list$theta <- mixedEffect_list$theta
    } else {
      mixedEffect_list       <- ME.startvalues(Y, mixedEffect_list)
      measurement_list$sigma <- mixedEffect_list$sigma  
    }
  }
  
  # Fitting the models: the case where the model formulation consists of W(t)
  debug.output <- list()
  if(use.process){
    
    if(is.null(init.fit) == TRUE ||
       (reffects == "Normal" & process[1] == "Normal" & error == "Normal")){
      
      #starting values for process
      #operator_list$type  <- process[2]
      if(bivariate){
        operator_list <- operator.startvalues.bivariate(Y,
                                              locs,
                                              mixedEffect_list,
                                              operator_list,
                                              measurement_list)
      } else {
        operator_list <- operator.startvalues(Y,
                                              locs,
                                              mixedEffect_list,
                                              operator_list,
                                              measurement_list)  
      }
      
    }
    
    # at least one of random effects, process or measurement error non-Gaussian
    if(reffects != "Normal" || process[1] != "Normal" || error != "Normal"){
      
      if(is.null(init.fit) == TRUE){
        
        # first fit the Gaussian model to obtain initials
        if(!silent){
          cat("Estimate Gaussian")
        }
        
        if(debug){
          debug.output$gaussian.input = list(Y =Y,
                                          locs = locs,
                                          mixedEffect_list = mixedEffect_list,
                                          measurement_list = measurement_list,
                                          process_list = process_list,
                                          operator_list = operator_list)
        }
        fit <- estimateLong(Y,
                            locs,
                            mixedEffect_list,
                            measurement_list,
                            process_list,
                            operator_list,
                            nIter = controls.init$nIter.init,
                            silent = silent,
                            learning_rate = controls.init$learning.rate.init,
                            polyak_rate = controls.init$polyak.rate.init,
                            nBurnin = controls.init$nBurnin.init,
                            nSim = controls.init$nSim.init,
                            pSubsample = controls.init$pSubsample.init,
                            nPar_burnin = controls.init$nPar.burnin.init,
                            step0 = controls.init$step0.init,
                            alpha = controls.init$alpha.init,
                            nBurnin_learningrate = controls.init$nBurnin.learningrate.init,
                            nBurnin_base = controls.init$nBurnin.base.init,
                            subsample.type = controls.init$subsample.type.init,
                            pSubsample2 = controls.init$pSubsample2.init,
                            seed = controls.init$seed.init)
        
      }else{
        fit <- init.fit
      }
      
      # then fit the non-Gaussian model
      if(!silent){
        cat("Estimate non-Gaussian")
      }
      if(use.random){
        if(fit$mixedEffect_list$noise == "Normal" && reffects != "Normal"){
          if(dim(fit$mixedEffect_list$U)[1]==1){
            lambda <- function(x) -sum(dnig(c(fit$mixedEffect_list$U),0,0,exp(x[1]),exp(x[2]),log = T))
            res <- optim(c(
              0,
              log(fit$mixedEffect_list$Sigma)/2),
              lambda)
            fit$mixedEffect_list$mu <- matrix(0, dim(B_random[[1]])[2], 1)
            fit$mixedEffect_list$nu <- max(as.matrix(exp(res$par[1])),600)
            fit$mixedEffect_list$Sigma <- as.matrix(exp(2*res$par[2]))
            if(fit$mixedEffect_list$nu< 10){
              lambda <- function(x) -sum(dnig(c(fit$mixedEffect_list$U),-x[1],x[1],exp(x[2]),exp(x[3]),log = T))
              res <- optim(c(0,
                             0,
                             log(fit$mixedEffect_list$Sigma)/2),
                           lambda)
              fit$mixedEffect_list$mu <- res$par[1]
              fit$mixedEffect_list$nu <- max(as.matrix(exp(res$par[2])),600)
              fit$mixedEffect_list$Sigma <- as.matrix(exp(2*res$par[3]))
            }
          }else{
            fit$mixedEffect_list$nu <- as.matrix(3.)
            fit$mixedEffect_list$mu <- matrix(0, dim(B_random[[1]])[2], 1)
          }
        }
        fit$mixedEffect_list$noise <- reffects  
      }
      
      if(fit$processes_list$noise == "Normal" && process[1] != "Normal"){
        if(bivariate){
          fit$processes_list$mu <- as.matrix(c(0,0))
          fit$processes_list$nu <- as.matrix(c(1,1))  
        } else {
          fit$processes_list$mu <- 0
          fit$processes_list$nu <- 1      
        }
        
      }
      fit$processes_list$noise <- process[1]
      
      if(length(process) < 3){
        nu_limit = 0
      }else {
        nu_limit = as.double(process[3])
      }
      
      fit$processes_list$nu_limit <- nu_limit
      
      
      if(fit$measurementError_list$noise == "Normal" && error != "Normal"){
        fit$measurementError_list$nu  <-  3.
        fit$measurementError_list$mu  <-  0
        fit$measurementError_list$Vs  <- Vin
      }
      fit$measurementError_list$noise      <- error
      fit$measurementError_list$assymetric <- error_assymetric
      
      fit$measurementError_list$common_V <- controls$individual.sigma
      
      if(debug){
        debug.output$nongaussian.input = list(Y =Y,
                                           locs = locs,
                                           mixedEffect_list = fit$mixedEffect_list,
                                           measurement_list = fit$measurement_list,
                                           process_list = fit$process_list,
                                           operator_list = fit$operator_list)
      }
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
                          seed = controls$seed)
    } else { # random effects, process and measurement error are all Gaussian
      
      if(!silent){
        cat("Estimate Model")
      }
      
      # Obtain parameter estimates
      # Obtain parameter estimates
      if(debug){
        debug.output$gaussian.input = list(Y =Y,
                                           locs = locs,
                                           mixedEffect_list = mixedEffect_list,
                                           measurement_list = measurement_list,
                                           process_list = process_list,
                                           operator_list = operator_list)
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
                          seed = controls$seed)
    }
  } else {## Fitting the models: the case where W(t) is excluded
    
    # either random effects or measurement error is non-Gaussian
    if(reffects != "Normal" || error != "Normal"){
      
      if(is.null(init.fit) == TRUE){
        
        if(!silent){
          cat("Estimate Gaussian")
        }
        
        # first fit the Gaussian model to obtain the initials
        fit <- estimateLong(Y,
                            locs,
                            mixedEffect_list,
                            measurement_list,
                            nIter = nIter,
                            silent = silent,
                            learning_rate = controls.init$learning.rate.init,
                            polyak_rate = controls.init$polyak.rate.init,
                            nBurnin = controls.init$nBurnin.init,
                            nSim = controls.init$nSim.init,
                            pSubsample = controls.init$pSubsample.init,
                            nPar_burnin = controls.init$nPar.burnin.init,
                            step0 = controls.init$step0.init,
                            alpha = controls.init$alpha.init,
                            nBurnin_learningrate = controls.init$nBurnin.learningrate.init,
                            nBurnin_base = controls.init$nBurnin.base.init,
                            subsample.type = controls.init$subsample.type.init,
                            pSubsample2 = controls.init$pSubsample2.init,
                            seed = controls.init$seed.init)
        
      }else{
        
        fit <- init.fit
        
      }
      
      # then fit the non-Gaussian model
      if(!silent){
        cat("Estimate non-Gaussian")
      }
      
      fit$mixedEffect_list$noise <- reffects
      fit$mixedEffect_list$nu    <- as.matrix(10)
      fit$mixedEffect_list$mu    <- matrix(0, dim(B_random[[1]])[2], 1)
      
      fit$measurementError_list$noise    <- error
      fit$measurementError_list$nu       <- 3.
      fit$measurementError_list$assymetric <- error_assymetric
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
                          seed = controls$seed)
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
                          seed = controls$seed)
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
  if(use.random){
    fixed_est1 <- as.numeric(fit$mixedEffect_list$beta_fixed)
    fixed_est2 <- as.numeric(fit$mixedEffect_list$beta_random)
    names(fixed_est1) <- colnames(x_fixed)
    names(fixed_est2) <- to_del_x_fixed
    fixed_est <- rep(NA, ncol(x_fixed_f))
    names(fixed_est) <- colnames(x_fixed_f)
    
    index_fixed  <- which(names(fixed_est) %in% names(fixed_est1))
    index_random <- which(names(fixed_est) %in% names(fixed_est2))
    
    fixed_est[index_fixed]  <- fixed_est1
    fixed_est[index_random] <- fixed_est2
    
    fixed_est1_vec <- fit$mixedEffect_list$betaf_vec
    fixed_est2_vec <- fit$mixedEffect_list$betar_vec
    colnames(fixed_est1_vec) <- colnames(x_fixed)
    colnames(fixed_est2_vec) <- to_del_x_fixed
    fixed_est_vec <- matrix(NA, ncol = ncol(x_fixed_f), nrow = nIter)
    colnames(fixed_est_vec) <- colnames(x_fixed_f)
    fixed_est_vec[, index_fixed]  <- fixed_est1_vec
    fixed_est_vec[, index_random] <- fixed_est2_vec
    
    # random effects
    ranef_Sigma           <- fit$mixedEffect_list$Sigma
    colnames(ranef_Sigma) <- rownames(ranef_Sigma) <- colnames(x_random)
    ranef_Sigma_vec       <- fit$mixedEffect_list$Sigma_vec
    ranef_sigma_epsilon   <- fit$mixedEffect_list$Sigma_epsilon  
  } else {
    fixed_est <- as.numeric(fit$mixedEffect_list$beta_fixed)
    names(fixed_est) <- colnames(x_fixed)
    fixed_est_vec <- fit$mixedEffect_list$betaf_vec
    colnames(fixed_est_vec) <- colnames(x_fixed)
    index_fixed  <- which(names(fixed_est) %in% names(fixed_est))
    index_random <- NA
    ranef_Sigma       <- NA
    ranef_Sigma_vec       <- NA
    ranef_sigma_epsilon   <- NA
    
  }
  
  if(reffects %in% c("NIG","tdist")){
    ranef_mu     <- fit$mixedEffect_list$mu
    ranef_mu_vec <- fit$mixedEffect_list$mu_vec
    
    ranef_nu <- fit$mixedEffect_list$nu
    ranef_nu_vec <- fit$mixedEffect_list$nu_vec
    
    if(reffects %in% c("NIG") && rev(ranef_nu_vec)[1] == 100){
      warning("nu = 100 indicates NIG has converged to Normal")
    }
    
  }else{
    ranef_mu <- ranef_mu_vec <- ranef_nu <- ranef_nu_vec <- NA
  }
  A = NULL
  # process and operator
  if(use.process == TRUE){
    A <- fit$A
    #operator
    operator_tau     <- fit$operator_list$tau
    operator_tau_vec <- fit$operator_list$tauVec
    operator_kappa <- fit$operator_list$kappa
    operator_kappa_vec <- fit$operator_list$kappaVec
      
    #process
    if(process[1] %in% c("NIG", "GAL")){
      process_nu <- fit$processes_list$nu
      process_nu_vec <- fit$processes_list$nu_vec
      
      if(process[1] == "NIG" && rev(process_nu_vec)[1] == 100){
        warning("nu = 100 indicates NIG has converged to Normal")
      }
      
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
  if(bivariate){
    meas_error_sigma <- exp(fit$measurementError_list$theta)
    meas_error_sigma_vec <- exp(fit$measurementError_list$theta_vec)
  } else {
    meas_error_sigma <- fit$measurementError_list$sigma
    meas_error_sigma_vec <- fit$measurementError_list$sigma_vec  
  }
  
  
  if(error %in% c("NIG", "tdist")){
    meas_error_nu <- fit$measurementError_list$nu
    meas_error_nu_vec <- fit$measurementError_list$nu_vec
    
    if(error == "NIG" && rev(meas_error_nu_vec)[1] == 100){
      warning("nu = 100 indicates NIG has converged to Normal")
    }
    
  }else{
    meas_error_nu <- meas_error_nu_vec <- NA
  }
  
  fisher_est <- NA
  
  out <- list(
    use_process = use.process,
    x_fixed_f = x_fixed_f,
    x_random = x_random,
    random_distr = reffects,
    process_distr = process[1],
    operator_type = process[2],
    error_distr = error,
    Y = Y,
    locs = locs,
    A = A,
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
    fisher_est = NA,
    call = match.call(),
    index_fixed = index_fixed,
    index_random = index_random,
    debug = debug.output
  )
  
  class(out) <- "ngme.spatial"
  out
  
}
