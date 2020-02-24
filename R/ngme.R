
#' @title Parameter estimation.
#'
#' @description Estimates model parameters for longitudianl models using maximum likelihood
#'   implemented by a computationally efficient stochastic gradient algorithm. 
#'   See \code{\link{ngme.spatial}} for estimation of spatial models. 
#'
#' @param fixed A two-sided formula to specify the fixed effects design matrix.
#' @param random A one-sided formula to specify the random effects design matrix.
#' @param use.process A logical variable for inclusion of the stochastic process in
#'   the mixed model: \code{"TRUE"} indicates inclusion, \code{"FALSE"} exclusion.
#' @param reffects A character string that indicates the distribution of the
#'   random effects. Available options are:
#'   \code{"Normal"} for Normal distribution, and
#'   \code{"NIG"} for Normal-inverse Gaussian.
#' @param process A character vector with two elements to specify
#'   the process. Whilst the first element is for the covariance structure, the
#'   second element for the process distribution. Available options for the
#'   first are:
#'   \code{"fd2"} for integrated Random-Walk (integrated Brownian Motion),
#'   \code{"matern"} for Matern covarince with smoothnes 3/2,
#'   \code{"exponential"} for exponential covarince (Matern with smoothnes 1/2);
#'   for the second element are:
#'   \code{"Normal"} for Normal distribution,
#'   \code{"NIG"} for Normal-inverse Gaussian,
#'   \code{"GAL"} for generalised-asymmetric Laplace, and
#'   \code{"CH"} for Cauchy.
#'   \code{process} is ignored when \code{use.process} is set to FALSE.
#' @param error A character string to specify the distribution of the error term.
#'   Available options are:
#'   \code{"Normal"} for Normal distribution,
#'   \code{"NIG"} for Normal-inverse Gaussian,
#'   \code{"tdist"} for t.
#' @param error_assymetric if true the non-Gaussian error is assymetric
#' @param data A data-frame from which the response and covariates to be extracted.
#' @param timevar A character string that indicates the column name of the time variable
#'   in \code{data}.
#' @param silent A logical value for printing the details of the iterations;
#'   \code{"TRUE"} indicates do not print, \code{"FALSE"} print.
#' @param nIter A numeric value for the number of iteration that will be
#'   used by the stochastic gradient.
#' @param mesh A list of control variables for creating mesh.
#'  \itemize{
#'    \item \code{max.dist}
#'    \item \code{cutoff} A numeric value. All time points that are separated less
#'    than \code{cutoff} will be merged to one node.
#'    \item \code{commond.grid} A logical value for creating same grids for different
#'    subjects. \code{"TRUE"} indicates common grid, \code{"FALSE"} uncommon grid.
#'    \item \code{extend} A numeric value or two element numeric vector for
#'    extending the meshes beyond the measurement locations with the specified
#'    percentage(s). If \code{extend}
#'    is specified by a single value, it indicates extending the mesh towards
#'    left and right by the same amount. If specified by a two element vector,
#'    whilst the first element is for extending towards the left, the second is
#'    for extending towards the right. If you want to extend by a fixed amount instead of with a percentage of the 
#'    interval length, use a negative value of \code{extend}. So for example, \code{extend = c(-1,0.1)} means that 
#'    the interval is extend by 1 to the left and by 10 percent to the right. 
#'    \item \code{n.cores} A numeric value for the number of cores to be used to create the
#'    mesh.
#'  }
#' @param controls A list of control variables for parameter estimation.
#'  \itemize{
#'     \item \code{learning.rate} A numeric value for the parameter of stochastic gradient.
#'     \item \code{polyak.rate} A numeric value for moving average of parameters;
#'       -1: inactive, 0: pure mean.
#'     \item \code{nBurnin} A numeric value for the number of steps before starting
#'       gradient estimation.
#'     \item \code{nSim} A numeric value for the number of samples of the Gibbs sampler
#'       to estimate the gradient.
#'     \item \code{pSubsample} A numeric value for the portion of data to be used in each
#'       gradient iteration. \code{pSubsample = 1} indicates use of all subjects' data.
#'     \item \code{nPar.burnin} A numeric value; "M-step" updates will be used until this
#'       iteration.
#'     \item \code{nIter.fisher} A numeric value for the number of iterations to be used to
#'       obtain the Fisher-Information matrix.
#'     \item \code{nSim.fisher} A numeric value for the number of samples of the Gibbs sampler
#'       to obtain the Fisher-Information matrix.
#'     \item \code{step0} A numeric value for step-size of the optimizer;
#'       where step-size is defined as step0 / i^alpha, i being the iteration,
#'       \code{alpha} is another tuning parameter specified next.
#'     \item \code{alpha} A numeric value for stepsize of the optimizer; step0 / i^alpha.
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
#'     \item \code{standardize.mixedEffects} A logical variable for standardising the covariates;
#'       \code{"FALSE"} indicates no standardisation, \code{"TRUE"} standardisation.
#'     \item \code{estimate.fisher} A logical variable or numeric scalar
#'      for whether Fisher-Information matrix
#'       to be obtained; \code{"FALSE"} indicates do not obtain, \code{1} expected Fisher
#'       matrix, \code{2} for observed.
#'     \item \code{individual.sigma} A logical variable for specifying patient-specific mixture
#'       random variable for the error-term; \code{"FALSE"} indicates do not obtain,
#'       \code{"TRUE"} obtain.
#'  }
#' @param controls.init A list of control variables to be used to fit the normal model
#'      to get the initial values for fitting a model with at least one of random effects,
#'      process and error being non-Gaussian.
#'  \itemize{
#'  \item \code{learning.rate.init} See \code{learning.rate} in \code{controls}.
#'  \item \code{polyak.rate.init} See \code{polyak.rate} in \code{controls}.
#'  \item \code{nBurnin.init} See \code{nBurnin} in \code{controls}.
#'  \item \code{nSim.init} See \code{nSim} in \code{controls}.
#'  \item \code{nIter.init} See \code{nIter} in \code{controls}.
#'  \item \code{pSubsample.init} See \code{pSubsample} in \code{controls}.
#'  \item \code{nPar.burnin.init} See \code{nPar.burnin} in \code{controls}.
#'  \item \code{step0.init} See \code{step0} in \code{controls}.
#'  \item \code{alpha.init} See \code{alpha} in \code{controls}.
#'  \item \code{nBurnin.learningrate.init} See \code{nBurnin.learningrate} in \code{controls}.
#'  \item \code{nBurnin.base.init} See \code{nBurnin.base} in \code{controls}.
#'  \item \code{subsample.type.init} See \code{subsample.type} in \code{controls}.
#'  \item \code{pSubsample2.init} See \code{pSubsample2} in \code{controls}.
#'  \item \code{standardize.mixedEffects.init} See \code{standardize.mixedEffects} in \code{controls}.
#'  \item \code{individual.sigma.init = FALSE} See \code{individual.sigma} in \code{controls}.
#'  }
#' @param init.fit A fitted \code{ngme} object with normal distribution for random effects,
#'  process and error.
#' @details This function is a user-friendly wrapper that calls the \code{estimateLong} function.
#'     Generic functions \code{summary}, \code{print} and \code{plot} are available for the
#'     output returned by the function \code{ngme}. For Matern covariance function,
#'     currently the shape parameter is set to 0.5 which corresponds to exponential correlation
#'     function.
#' @return A list of outputs.
#' @seealso \code{\link{ngme.spatial}} 
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'
#'   # transform pwl to decrese it correlation with bage
#'   # then center all the covariates for better convergence
#'   srft_data$pwl2 <- with(srft_data, pwl - bage/1.5)
#'   srft_data[, c("sex_cent", "bage_cent", "fu_cent", "pwl2_cent")] <-
#'     scale(srft_data[, c("sex", "bage", "fu", "pwl2")], scale = FALSE)
#'
#'   # fit the model with normal assumption for random effects, process and error
#'   # covariance function is integrated random walk
#'   set.seed(123)
#'   fit_normal_normal_normal_fd2 <- ngme(fixed = log(egfr) ~ sex_cent + bage_cent + fu_cent + pwl2_cent,
#'                                        random = ~ 1|id,
#'                                        data = srft_data,
#'                                        reffects = "Normal",
#'                                        process = c("Normal", "fd2"),
#'                                        error = "Normal",
#'                                        timeVar = "fu",
#'                                        nIter = 20000,
#'                                        use.process = TRUE,
#'                                        silent = FALSE,
#'                                        mesh = list(cutoff = 1/365,
#'                                                    max.dist = 1/12,
#'                                                    extend = 0.01
#'                                                    ),
#'                                        controls = list(pSubsample = 0.025,
#'                                                        step0 = 1,
#'                                                        estimate.fisher = FALSE,
#'                                                        subsample.type = 1,
#'                                                        polyak.rate = 0.01,
#'                                                        alpha = 0.01
#'                                                        )
#'                                        )
#' # fit the model with NIG assumption for all the random components
#' set.seed(123)
#' fit_nig_nig_nig_fd2 <- update(fit_normal_normal_normal_fd2,
#'                               reffects = "NIG",
#'                               process = c("NIG", "fd2"),
#'                               error = "NIG",
#'                               init.fit = fit_normal_normal_normal_fd2
#'                               )
#' }

ngme <- function(fixed,
                 random = NULL,
                 use.process = FALSE,
                 reffects = "Normal",
                 process = c("Normal", "fd2"),
                 error = "Normal",
                 error_assymetric = FALSE,
                 data,
                 timeVar = NULL,
                 silent = TRUE,
                 nIter = 1000,
                 mesh = list(max.dist = NULL,
                             cutoff = NULL,#1e-10,
                             common.grid = FALSE,
                             extend = NULL,
                             n.cores = 1),
                 controls = list(learning.rate = 0.2,
                                 polyak.rate = -1,
                                 nBurnin = 100,
                                  nSim = 2,
                                 pSubsample = NULL,
                                 nPar.burnin = 0,
                                 nIter.fisher = 1000,
                                 nSim.fisher = 1000,
                                 step0 = 1,
                                 alpha = 0.6,
                                 nBurnin.learningrate = NULL,
                                 nBurnin.base = 0,
                                 subsample.type = 4,
                                 pSubsample2 = 0.3,
                                 standardize.mixedEffects = FALSE,
                                 estimate.fisher = FALSE,
                                 individual.sigma = FALSE,
                                 iter.start = 0),
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
                                      standardize.mixedEffects.init = FALSE,
                                      individual.sigma.init = FALSE),
                 init.fit = NULL
                 )
{
  
  # being sure that controls includes everything
  if(length(controls) < 18){
    controls.full <- list(learning.rate = 0.2,
                          polyak.rate = -1,
                          nBurnin = 100,
                          nSim = 2,
                          pSubsample = NULL,
                          nPar.burnin = 0,
                          nIter.fisher = 1000,
                          nSim.fisher = 1000,
                          step0 = 1,
                          alpha = 0.6,
                          nBurnin.learningrate = NULL,
                          nBurnin.base = 0,
                          subsample.type = 4,
                          pSubsample2 = 0.3,
                          standardize.mixedEffects = FALSE,
                          estimate.fisher = FALSE,
                          individual.sigma = FALSE,
                          iter.start = 0
                          )
    for(i in 1:length(controls.full)){
      if(!(names(controls.full)[i] %in% names(controls))){
        controls[names(controls.full)[i]] <- controls.full[i]
      }
    }

  }
  # generate a seed
  gen_seed <- ceiling(10^8 * runif(1))
  controls$seed <-  gen_seed

  # check for controls.init
  if(is.null(controls.init) || length(controls.init) < 16){
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
                                 standardize.mixedEffects.init = FALSE,
                                 individual.sigma.init = FALSE)
      for(i in 1:length(controls.init.full)){
        if(!(names(controls.init.full)[i] %in% names(controls.init))){
          controls.init[names(controls.init.full)[i]] <- controls.init.full[i]
        }
      }
      controls.init$seed.init <- gen_seed
    }

  

  ## check mesh
  if(length(mesh) < 5){
    mesh.full = list(max.dist = NULL,
                     cutoff = NULL,
                     common.grid = FALSE,
                     extend = NULL,
                     n.cores = 1
                     )
    for(i in 1:length(mesh.full)){
      if(!(names(mesh.full)[i] %in% names(mesh))){
        mesh[names(mesh.full)[i]] <- mesh.full[i]
      }
    }
  }

  # return an error if max.dist and cutoff not provided
  if(use.process == TRUE){

    if((reffects == "Normal" & process[1] == "Normal" & error == "Normal") ||
       (is.null(init.fit) == TRUE & (reffects != "Normal" & process[1] != "Normal" & error != "Normal"))){

      if(is.null(mesh$max.dist) == TRUE & is.null(mesh$cutoff) == TRUE){
        stop("Provide 'max.dist' and 'cutoff' for creating mesh")
      }

    }

  }else{
    locs = NULL
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
  if(use.process == TRUE & is.null(timeVar) == TRUE){
    stop("'timeVar' should be specified, since the model consists of process")
  }

  # alpha and alpha.init are in the correct interval?
  if(controls$alpha < 0 | controls$alpha > 1){
    stop("alpha should be in (0, 1]")
  }
  if(controls.init$alpha.init < 0 | controls.init$alpha.init > 1){
    stop("alpha.init should be in (0, 1]")
  }
  
  # extract id variable
  idname <- rev(unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE)))[1]
  id <- data[, idname]

  effects <- extract.effects(data = data, 
                              fixed = fixed, 
                              random = random,
                              idname = idname)
  Y = effects$Y
  B_random  = effects$B_random
  B_fixed   = effects$B_fixed
  x_fixed   = effects$x_fixed
  x_random  = effects$x_random
  x_fixed_f = effects$x_fixed_f
  to_del_x_fixed = effects$to_del_x_fixed
    
  # extract variables for process
  if(use.process == TRUE){
    locs <- tapply(as.matrix(data[, timeVar]), id, function(x) x)  
  }

  nsubj <- length(Y)
  
  # if pSubsampling not set
  if(is.null(controls.init$pSubsample.init)){
    if(nsubj < 100){
      controls.init$pSubsample.init = 1
      warning("pSubsample.init not provided. Since there are <100 subjects, p for subsampling is set to 1")
    }else if(nsubj < 500){
      controls.init$pSubsample.init = 0.2
      warning("pSubsample.init not provided. Since there are >=100 and <500 subjects, p for subsampling is set to 0.2")
    }else{
      controls.init$pSubsample.init = 0.1
      warning("pSubsample.init not provided. Since there >=500 subjects, p for subsampling is set to 0.1")
    }
  }
  if(is.null(controls$pSubsample)){
    if(nsubj < 100){
      controls$pSubsample = 1
      warning("pSubsample not provided. Since there are <100 subjects, p for subsampling is set to 1")
    }else if(nsubj < 500){
      controls$pSubsample = 0.2
      warning("pSubsample not provided. Since there are >=100 and <500 subjects, p for subsampling is set to 0.2")
    }else{
      controls$pSubsample = 0.1
      warning("pSubsample not provided. Since there >=500 subjects, p for subsampling is set to 0.1")
    }
  }
  
  
  ###
  # TODO add something if init.fit is non Gaussian
  ###
  ## Vin is needed even if init.fit is not NULL
  Vin <- lapply(Y, function(x) rep(1, length(x)))

  ## a warning message
  if(use.process == TRUE & is.null(init.fit) == FALSE){
    warning("'cutoff', 'max.dist' and 'extend' for 'mesh' were inherited from the 'init.fit'")
  }

  ## Obtain starting values - if init.fit is not supplied or everything is Gaussian

  if(is.null(init.fit) == TRUE ||
     (reffects == "Normal" & (use.process == TRUE & process[1] == "Normal") & error == "Normal") ||
     (reffects == "Normal" & use.process == FALSE & error == "Normal")){

    # setup the lists that are going to be passed into the
    # functions that will obtain the starting value
    if(!silent){
      cat("Setup lists\n")
    }
    if(!is.null(init.fit)){
      measurement_list <- init.fit$measurementError_list
      mixedEffect_list <- init.fit$mixedEffect_list
      
      if(use.process){
        operator_list <- init.fit$operator_list
        process_list <- init.fit$processes_list
      }
    } else {
      measurement_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1)
      mixedEffect_list <- list(B_random = B_random,
                               B_fixed  = B_fixed,
                               noise = "Normal",
                               Sigma_epsilon = 1)
      
      if(use.process){
        operator_list <- create_operator(locs,
                                         name = process[2],
                                         common.grid = mesh$common.grid,
                                         extend  = mesh$extend,
                                         max.dist = mesh$max.dist,
                                         cutoff = mesh$cutoff,
                                         n.cores = mesh$n.cores)
        
        process_list = list(noise = "Normal", nu  = 1, mu  = 0)
        process_list$V <- list()
        process_list$X <- list()
        
        for(i in 1:length(locs))
        {
          process_list$X[[i]] <- rep(0, length(operator_list$h[[i]]))
          process_list$V[[i]] <- operator_list$h[[i]]
        }  
      }
      # starting values for measurement error and mixed effects using OLS:
      if(!silent){
        cat("Calculate starting values\n")
      }
      
      mixedEffect_list       <- ME.startvalues(Y, mixedEffect_list)
      measurement_list$sigma <- mixedEffect_list$sigma
    }
    
  }

  # Fitting the models: the case where the model formulation consists of W(t)
  if(use.process){

    if(is.null(init.fit) == TRUE){
      operator_list <- operator.startvalues(Y,
                                            locs,
                                            mixedEffect_list,
                                            operator_list,
                                            measurement_list)
    }

    # at least one of random effects, process or measurement error non-Gaussian
    if(reffects != "Normal" || process[1] != "Normal" || error != "Normal"){

      if(is.null(init.fit) == TRUE){

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
                            seed = controls.init$seed.init,
                            standardize.mixedEffects = controls.init$standardize.mixedEffects.init,
                            estimate_fisher = FALSE)

      }else{
        fit <- init.fit
      }

      # then fit the non-Gaussian model
      if(!silent){
        cat("Estimate non-Gaussian")
      }

      if(fit$mixedEffect_list$noise == "Normal" && reffects != "Normal"){
        if(dim(fit$mixedEffect_list$U)[1]==1){
          lambda <- function(x) -sum(dnig(c(fit$mixedEffect_list$U),0,0,exp(x[1]),exp(x[2]),log = T))
          res <- optim(c(
                  0,
                  log(fit$mixedEffect_list$Sigma)/2),
                lambda)
          fit$mixedEffect_list$mu <- matrix(0, dim(B_random[[1]])[2], 1)
          fit$mixedEffect_list$nu <- max(as.matrix(exp(res$par[1])),10)
          fit$mixedEffect_list$Sigma <- as.matrix(exp(2*res$par[2]))
          if(fit$mixedEffect_list$nu< 10){
            lambda <- function(x) -sum(dnig(c(fit$mixedEffect_list$U),-x[1],x[1],exp(x[2]),exp(x[3]),log = T))
            res <- optim(c(0,
              0,
              log(fit$mixedEffect_list$Sigma)/2),
              lambda)
            fit$mixedEffect_list$mu <- res$par[1]
            fit$mixedEffect_list$nu <- max(as.matrix(exp(res$par[2])),10)
            fit$mixedEffect_list$Sigma <- as.matrix(exp(2*res$par[3]))
          }
        }else{
          fit$mixedEffect_list$nu <- as.matrix(3.)
          fit$mixedEffect_list$mu <- matrix(0, dim(B_random[[1]])[2], 1)
        }
      }
      fit$mixedEffect_list$noise <- reffects
      
      if(fit$processes_list$noise == "Normal" && process[1] != "Normal"){
        if(0){
          cat("\n", "Estimate starting values for process parameters.\n")
          like.proc <- function(x,y,h)
          {
            nu = x[1]
            mu = x[2]
            tau = x[3]
            if(nu<0){
              return(Inf)
            } else {
              return(-sum(dnig(y, -mu*h/tau, mu/tau, nu*h, 1/tau,log=TRUE)))  
            }
          }
          
          noise <- unlist(fit$processes_list$W)
          h <- unlist(fit$operator_list$h)
          pars <- optim(c(1,0,1),like.proc,y=noise,h=h)
          cat("initial process parameters: nu = ", pars$par[1], ", mu = ", pars$par[2], "tau_fix = ", pars$par[3],"\n")
          if(pars$par[1]<10 && abs(pars$par[2]) < 10){
            fit$processes_list$nu <- pars$par[1]  
            fit$processes_list$mu <- pars$par[2]
            fit$operator_list$tau <- fit$operator_list$tau*pars$par[3]
          } else {
            fit$processes_list$mu <- 0
            fit$processes_list$nu <- 10  
          }  
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
      }
      fit$measurementError_list$noise      <- error
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
                          iter_start = controls$iter.start,
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
      if(controls$estimate.fisher > 0){
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
                              estimate_fisher = controls$estimate.fisher
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
                          iter_start = controls$iter.start,
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
      if(controls$estimate.fisher > 0){
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
                              estimate_fisher = controls$estimate.fisher
                              )
        fit$FisherMatrix <- fit.f$FisherMatrix
      }
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
                            seed = controls.init$seed.init,
                            standardize.mixedEffects = controls.init$standardize.mixedEffects.init,
                            estimate_fisher = FALSE)

      }else{

        fit <- init.fit

      }

      # then fit the non-Gaussian model
      if(!silent){
        cat("Estimate non-Gaussian")
      }
      if(fit$mixedEffect_list$noise != reffects){
        fit$mixedEffect_list$noise <- reffects
        fit$mixedEffect_list$nu    <- as.matrix(10)
        fit$mixedEffect_list$mu    <- matrix(0, dim(B_random[[1]])[2], 1)
        
        fit$measurementError_list$noise    <- error
        fit$measurementError_list$nu       <- 3.
        fit$measurementError_list$assymetric <- error_assymetric
        fit$measurementError_list$common_V <- controls$individual.sigma
        fit$measurementError_list$Vs       <- Vin  
      }
      

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
                          iter_start = controls$iter.start,
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

      if(controls$estimate.fisher > 0){
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
                              estimate_fisher = controls$estimate.fisher
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
                          iter_start = controls$iter.start,
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
      if(controls$estimate.fisher > 0){
        fit.f <- estimateLong(Y,
                              locs,
                              fit$mixedEffect_list,
                              fit$measurementError_list,
                              learning_rate = controls$learning.rate,
                              nIter = controls$nIter.fisher,
                              silent = silent,
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
                              estimate_fisher = controls$estimate.fisher
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
  fixed_est1    <- as.numeric(fit$mixedEffect_list$beta_fixed)
  fixed_est2    <- as.numeric(fit$mixedEffect_list$beta_random)
  index_fixed   = 1:length(fixed_est1)
  index_random  = length(fixed_est1) + (1:length(fixed_est2))
  names(fixed_est1) <- colnames(x_fixed)
  names(fixed_est2) <- colnames(x_random)
  fixed_est <- c(fixed_est1,fixed_est2)

  fixed_est1_vec <- fit$mixedEffect_list$betaf_vec
  fixed_est2_vec <- fit$mixedEffect_list$betar_vec
  colnames(fixed_est1_vec) <- colnames(x_fixed)
  colnames(fixed_est2_vec) <- colnames(x_random)
  fixed_est_vec <- cbind(fixed_est1_vec,fixed_est2_vec)

  # random effects
  ranef_Sigma           <- fit$mixedEffect_list$Sigma
  colnames(ranef_Sigma) <- rownames(ranef_Sigma) <- colnames(x_random)
  ranef_Sigma_vec       <- fit$mixedEffect_list$Sigma_vec
  ranef_sigma_epsilon   <- fit$mixedEffect_list$Sigma_epsilon

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

    if(process[2] %in% c("fd2")){
      operator_tau_vec <- operator_tau_vec[-1] #exclude the first element
    }

    #operator - matern
    if(process[2] %in% c("matern", "exponential","matern.asym")){
      operator_kappa <- fit$operator_list$kappa
      operator_kappa_vec <- fit$operator_list$kappaVec

    }else{
      operator_kappa <- operator_kappa_vec <- NA
    }

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
  meas_error_sigma <- fit$measurementError_list$sigma
  meas_error_sigma_vec <- fit$measurementError_list$sigma_vec

  if(error %in% c("NIG", "tdist")){
    meas_error_nu <- fit$measurementError_list$nu
    meas_error_nu_vec <- fit$measurementError_list$nu_vec

    if(error == "NIG" && rev(meas_error_nu_vec)[1] == 100){
      warning("nu = 100 indicates NIG has converged to Normal")
    }

  }else{
    meas_error_nu <- meas_error_nu_vec <- NA
  }

  # checking Fisher Matrix is estimated?
  if(controls$estimate.fisher > 0){

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
    fisher_est = fisher_est,
    call = match.call(),
    index_fixed = index_fixed,
    index_random = index_random
  )
  
  class(out) <- "ngme"
  out

}
