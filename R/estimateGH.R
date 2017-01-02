estimate.wrapper <- function(Y,
                             locs,
                             B_random,
                             B_fixed,
                             use.process = TRUE,
                             operator.type = "fd2",
                             n.process = NULL,
                             measurement.distribution = "Normal",
                             random.effect.distribution= "Normal",
                             process.distribution= "Normal",
                             individual.sigma = FALSE,
                             silent = FALSE,
                             estimation.options = NULL,
                             estimate_fisher = FALSE,
                             ...)
{

  estimation.controls = list(learning.rate = 0,
                             polyak_rate = 0.1,
                             nBurnin = 100,
                             nSim = 2,
                             nIter.gauss = 1000,
                             nIter = 10000,
                             pSubsample = 0.1,
                             nPar_burnin = 0)
  if(!missing(estimation.options) && !is.null(estimation.options)){
    for(i in 1:length(estimation.options)){
      estimation.controls[names(estimation.options)[i]] = estimation.options[i]
    }
}
  if(!silent)
    cat("Setup lists\n")
  Vin <- list()
  for(i in 1:n.pers)
  {
    Vin[[i]] <- rep(1, length(Y[[i]]))
  }
  measurement_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1)
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            noise = "Normal",
                            Sigma_epsilon=1)
  if(use.process){
    if(is.null(n.process)){
      n.process = max(round(mean(unlist(lapply(locs,length)))),1)
    }
    operator_list <- create_operator(locs, n.process, name = operator.type)

    process_list = list(noise = "Normal",
                        nu  = 1,
                        mu  = 0)
    process_list$V <- list()
    process_list$X <- list()
    for(i in 1:length(locs))
    {
      process_list$X[[i]] <- rep(0,length(operator_list$h[[1]]))
      process_list$V[[i]] <- operator_list$h[[1]]
    }
  }

  #starting values for measurement error and mixed effects using OLS:
  if(!silent)
    cat("Calculate starting values\n")
  mixedEffect_list <- ME.startvalues(Y,mixedEffect_list)
  measurement_list$sigma = mixedEffect_list$sigma

  if(use.process){
    #starting values for process:
    operator_list <- operator.startvalues(Y,locs,mixedEffect_list,operator_list,measurement_list)
    operator_list$type  <- operator.type

    if(random.effect.distribution != "Normal" || process.distribution != "Normal" || measurement.distribution != "Normal"){
      #estimate Gaussian process model
      if(!silent)
        cat("Estimate Gaussian")

        res <- estimateLong(Y, locs,
                            mixedEffect_list,
                            measurement_list,
                            process_list,
                            operator_list,
                            learning_rate = estimation.controls$learning.rate,
                            nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                            polyak_rate = 0,
                            nSim = 2,
                            nBurnin = estimation.controls$nBurnin,
                            nIter = estimation.controls$nIter,
                            nPar_burnin = estimation.controls$nPar_burnin,
                            pSubsample = estimation.controls$pSubsample,
                            silent = silent,
                            estimate_fisher = FALSE,
                            ...)

      if(!silent)
        cat("Estimate non-Gaussian")

      res$mixedEffect_list$noise = random.effect.distribution
      res$mixedEffect_list$nu = as.matrix(10)
      res$mixedEffect_list$mu = matrix(0,dim(B_random[[1]])[2],1)

      res$processes_list$noise = process.distribution
      res$processes_list$mu = 0
      res$processes_list$nu = 10

      res$measurementError_list$noise = measurement.distribution
      res$measurementError_list$nu = 10
      res$measurementError_list$common_V = individual.sigma
      res$measurementError_list$Vs = Vin

      res <- estimateLong(Y, locs,
                          res$mixedEffect_list,
                          res$measurementError_list,
                          res$processes_list,
                          res$operator_list,
                          learning_rate = estimation.controls$learning.rate,
                          nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                          polyak_rate = 0,
                          nBurnin = estimation.controls$nBurnin,
                          nIter = estimation.controls$nIter,
                          nPar_burnin = estimation.controls$nPar_burnin,
                          pSubsample = estimation.controls$pSubsample,
                          silent = silent,
                          estimate_fisher = estimate_fisher,
                          ...)

    } else {
      if(!silent)
        cat("Estimate Model")

      res <- estimateLong(Y, locs,
                          mixedEffect_list,
                          measurement_list,
                          process_list,
                          operator_list,
                          learning_rate = estimation.controls$learning.rate,
                          nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                          polyak_rate = 0,
                          nSim = 2,
                          nBurnin = estimation.controls$nBurnin,
                          nIter = estimation.controls$nIter,
                          nPar_burnin = estimation.controls$nPar_burnin,
                          pSubsample = estimation.controls$pSubsample,
                          silent = silent,
                          estimate_fisher = estimate_fisher,
                          ...)

      }
    } else {

      if(random.effect.distribution != "Normal" || measurement.distribution != "Normal"){
        if(!silent)
          cat("Estimate Gaussian")
        res <- estimateLong(Y, locs,
                            mixedEffect_list,
                            measurement_list,
                            learning_rate = estimation.controls$learning.rate,
                            nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                            polyak_rate = 0,
                            nSim = 2,
                            nBurnin = estimation.controls$nBurnin,
                            nIter = estimation.controls$nIter.gauss,
                            nPar_burnin = estimation.controls$nPar_burnin,
                            pSubsample = estimation.controls$pSubsample,
                            silent = silent,
                            estimate_fisher = FALSE,
                            ...)
          if(!silent)
          cat("Estimate non-Gaussian")

        res$mixedEffect_list$noise = random.effect.distribution
        res$mixedEffect_list$nu = as.matrix(10)
        res$mixedEffect_list$mu = matrix(0,dim(B_random[[1]])[2],1)

        res$measurementError_list$noise = measurement.distribution
        res$measurementError_list$nu = 10
        res$measurementError_list$common_V = individual.sigma
        res$measurementError_list$Vs = Vin

        res <- estimateLong(Y, locs,
                            res$mixedEffect_list,
                            res$measurementError_list,
                            res$processes_list,
                            res$operator_list,
                            learning_rate = estimation.controls$learning.rate,
                            nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                            polyak_rate = 0,
                            nBurnin = estimation.controls$nBurnin,
                            nIter = estimation.controls$nIter,
                            nPar_burnin = estimation.controls$nPar_burnin,
                            pSubsample = estimation.controls$pSubsample,
                            silent = silent,
                            estimate_fisher = estimate_fisher,
                            ...)

      } else {
        res <- estimateLong(Y, locs,
                            mixedEffect_list,
                            measurement_list,
                            learning_rate = estimation.controls$learning.rate,
                            nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                            polyak_rate = 0,
                            nSim = 2,
                            nBurnin = estimation.controls$nBurnin,
                            nIter = estimation.controls$nIter,
                            nPar_burnin = estimation.controls$nPar_burnin,
                            pSubsample = estimation.controls$pSubsample,
                            silent = silent,
                            estimate_fisher = estimate_fisher,
                            ...)
      }
  }

  return(res)
}

#' @param   Y           - list with the observations
#' @param   locs        - list with position of the observations (Y)
#' @param mixedEffect_list -
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
#'
#' @param processes_list  - for the stochastic noise driving the
#' @param noise           - either Normal, NIG or GAL (change name to type rather then noise)
#' @param nu              - shape parameter for NIG or GAL
#' @param mu              - assymetric parameter for NIG or GAL
#'
#' @param learning_rate   - parameter for sthocastic gradient
#' @param nBurnin_learningrate - don't start learning before
#' @param nPar_burnin - use "M-step" updates until this iteration.
#' @param polyak_rate     - taking moving average of parameters (-1 means inactive, 0 mean pure mean)
#' @param step0           - stepsize for optimizer is step0 / i^alpha
#' @param alpha           - stepsize for optimizer is step0 / i^alpha
#' @param pSubsample      - precentage of data used in each gradient subsampling
#' @param subsample.type  - Type of subsampling: 0 - uniform without replacement
#'                                               1 - sample size weighted
#'                                               3 - weighted sampling by gradient size
#' @param pSubsample2     - precentage of data used in each gradient subsampling weighted by gradient
#' @param nIter           - number of iteration of the stochastic gradient
#' @param nSim            - number of samples of the gibbs sampler to estimate the gradient
#' @parma silent          - print iteration info
#' @param seed            - (unsinged int) seed for debuging
estimateLong <- function(Y,
                         locs,
                         mixedEffect_list,
                         measurment_list,
                         processes_list,
                         operator_list,
                         step0 = 0.3,
                         alpha = 0.3,
                         learning_rate = 0,
                         nBurnin_learningrate = NULL,
                         nBurnin_base = 0,
                         pSubsample = 1.,
                         polyak_rate = -1.,
                         subsample.type = 1,
                         nPar_burnin = 0,
                         pSubsample2 = 0.3,
                         nIter = 10,     # iterations to run the stochastic gradient
                         nSim  = 1,
                         nBurnin = 10,   # steps before starting gradient estimation
                         silent  = FALSE, # print iteration info
                         seed    = NULL,
                         standardize.mixedEffects = FALSE,
                         estimate_fisher  = FALSE
                         )
{
  obs_list <- list()
  use.process = TRUE
  if(missing(processes_list) || is.null(processes_list)){
    use.process = FALSE
  }
  if(use.process){
    common.grid = FALSE
    if(length(operator_list$loc)==1){
      common.grid = TRUE
    }
  }
  for(i in 1:length(locs)){
    if(length(Y[[i]]) != length(locs[[i]])){
      stop("Length of Y and locs differ.")
    }
    obs_list[[i]] <- list(Y=Y[[i]],
                          locs = locs[[i]])
    if(use.process){
      if(common.grid){
        obs_list[[i]]$A = spde.A(locs[[i]],operator_list$loc[[1]],
                                 right.boundary = operator_list$right.boundary,
                                 left.boundary = operator_list$left.boundary)
      } else {
        obs_list[[i]]$A = spde.A(locs[[i]],operator_list$loc[[i]],
                                 right.boundary = operator_list$right.boundary,
                                 left.boundary = operator_list$left.boundary)
      }
    }
  }

  if(standardize.mixedEffects){
    Bf.list <- standardize.covariates(mixedEffect_list$B_fixed)
    if(!is.null(Bf.list)){
      mixedEffect_list$B_fixed <- Bf.list$B
      mixedEffect_list$beta_fixed <- scale.beta(mixedEffect_list$beta_fixed,Bf.list)
    }
    Br.list <- standardize.covariates(mixedEffect_list$B_random)
    if(!is.null(Br.list)){
      mixedEffect_list$B_random <- Br.list$B
      mixedEffect_list$beta_random <- scale.beta(mixedEffect_list$beta_random,Br.list)
      mixedEffect_list$Sigma <- scale.sigma(mixedEffect_list$Sigma,Br.list)
    }
  }


  input <- list( obs_list         = obs_list,
                 measurementError_list  = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 pSubsample       = pSubsample,
                 pSubsample2      = pSubsample2,
                 subsample_type   = subsample.type,
                 nIter            = nIter,     # iterations to run the stochastic gradient
                 nSim             = nSim,
                 nBurnin          = nBurnin,   # steps before starting gradient estimation
                 silent           = silent, # print iteration info)
                 step0            = step0,
                 nBurnin_base     = nBurnin_base,
                 nPar_burnin      = nPar_burnin,
                 alpha            = alpha,
                 learning_rate    = learning_rate,
                 polyak_rate      = polyak_rate,
                 estimate_fisher  = estimate_fisher)
  if(use.process){
    input$processes_list   = processes_list
    input$operator_list    = operator_list
    input$common.grid      = common.grid
  }

  if(is.null(nBurnin_learningrate) == FALSE)
    input$nBurnin_learningrate =  nBurnin_learningrate
  if(is.null(seed) == FALSE)
    input <- setseed_ME(input, seed)

  output <- estimateLong_cpp(input)

  if(standardize.mixedEffects){
    if(!is.null(Bf.list)){
      output$mixedEffect_list$beta_fixed <- scale.beta(output$mixedEffect_list$beta_fixed,Bf.list,inv=TRUE)
      output$mixedEffect_list$betaf_vec <- scale.beta(output$mixedEffect_list$betaf_vec,Bf.list,inv=TRUE)
    }
    if(!is.null(Br.list)){
      output$mixedEffect_list$beta_random <- scale.beta(output$mixedEffect_list$beta_random,Br.list,inv=TRUE)
      output$mixedEffect_list$betar_vec <- scale.beta(output$mixedEffect_list$betar_vec,Br.list,inv=TRUE)
      mixedEffect_list$Sigma <- scale.sigma(mixedEffect_list$Sigma,Br.list,inv=TRUE)
    }
  }

  if(use.process){
    output$operator_list$left.boundary <- operator_list$left.boundary
    output$operator_list$right.boundary <- operator_list$right.boundary
    output$operator_list$type <- operator_list$type
    if(operator_list$type == "Matern"){
      output$operator_list$G <- operator_list$G
      output$operator_list$C <- operator_list$C
      output$operator_list$loc = operator_list$loc
      output$operator_list$h = operator_list$h
    } else if(operator_list$type == "fd2"){
      output$operator_list$Q <- operator_list$Q
    }
  }

  return(output)
}

#'
#' estimating mixed effect model
#'
#'
#' @param   Y           - list with the observations
#' @param   locs        - list with position of the observations (Y)
#' @param mixedEffect_list -
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
#'
#' @param learning_rate   - parameter for sthocastic gradient
#' @param polyak_rate     - taking moving average of parameters (-1 means inactive, 0 mean pure mean)
#' @param step0           - stepsize for optimizer is step0 / i^alpha
#' @param alpha           - stepsize for optimizer is step0 / i^alpha
#' @param pSubsample      - precentage of data used in each gradient subsampling
#' @param subsample.type  - Type of subsampling: 0 - uniform without replacement
#'                                               1 - sample size weighted
#' @param nIter           - number of iteration of the stochastic gradient
#' @param nSim            - number of samples of the gibbs sampler to estimate the gradient
#' @parma silent          - print iteration info
#' @param seed            - (unsinged int) seed for debuging
estimateME <- function(Y,
                         mixedEffect_list,
                         measurment_list,
                         step0 = 0.3,
                         alpha = 0.3,
                         learning_rate = 0,
                         nBurnin_learningrate = NULL,
                         pSubsample = 1.,
                         pSubsample2 = 0.3,
                         polyak_rate = -1.,
                         subsample.type = 1,
                         nPar_burnin = 0,
                         nBurnin_base = 0,
                         nIter = 10,     # iterations to run the stochastic gradient
                         nSim  = 1,
                         nBurnin = 10,   # steps before starting gradient estimation
                         silent  = FALSE, # print iteration info
                         seed = NULL,
                         estimate_fisher =FALSE
)
{
  obs_list <- list()
  for(i in 1:length(Y)){
      obs_list[[i]] <- list(Y = Y[[i]])
  }



  input <- list( obs_list         = obs_list,
                 measurementError_list  = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 pSubsample       = pSubsample,
                 pSubsample2      = pSubsample2,
                 subsample_type   = subsample.type,
                 nIter            = nIter,     # iterations to run the stochastic gradient
                 nSim             = nSim,
                 nBurnin          = nBurnin,   # steps before starting gradient estimation
                 silent           = silent, # print iteration info)
                 step0            = step0,
                 alpha            = alpha,
                 nPar_burnin      = nPar_burnin,
                 nBurnin_base = nBurnin_base,
                 learning_rate    = learning_rate,
                 polyak_rate      = polyak_rate,
                 estimate_fisher  = estimate_fisher
  )

  if(is.null(nBurnin_learningrate) == FALSE)
    input$nBurnin_learningrate =  nBurnin_learningrate

  if(is.null(seed) == FALSE)
    input <- setseed_ME(input, seed)

  output <- estimateLong_cpp(input)

  return(output)
}

setseed_ME <- function(input, seed)
{
  seed.old <- sample.int(10^6, 1)
  set.seed(seed)
  input$seed                  <- sample.int(10^6, 1)
  input$mixedEffect_list$seed <- sample.int(10^6, 1)
  set.seed(seed.old)
  return(input)
}

# Treat random effects as fixed effects to obtain start values using OLS
ME.startvalues <- function(Y,mixedEffect_list)
{
  n = length(mixedEffect_list$B_fixed)
  nc.f= dim(mixedEffect_list$B_fixed[[1]])[2]
  nc.r = 0
  if(!is.null(mixedEffect_list$B_random)){
    nc.r= dim(mixedEffect_list$B_random[[1]])[2]
  }
  nc = nc.f + nc.r
  BB = matrix(0,nc,nc)
  BY = matrix(0,nc,1)

  if(nc.r  > 0){
    beta_r = matrix(0,n,nc.r)
    I = diag(1,nc)

    for(i in 1:n){

      Bi = cBind(mixedEffect_list$B_fixed[[i]],mixedEffect_list$B_random[[i]])
      BB = BB + t(Bi)%*%Bi
      BY = BY + t(Bi)%*%Y[[i]]
    }
    beta = solve(BB,BY)
    mixedEffect_list$beta_fixed = beta[1:nc.f]
    mixedEffect_list$beta_random = beta[(nc.f+1):nc]
    res = NULL
    Sigma = matrix(0,nc.r,nc.r)
    br = matrix(0,n,nc.r)
    for(i in 1:n){
      BB = solve(t(mixedEffect_list$B_random[[i]])%*%mixedEffect_list$B_random[[i]])
      br[i,] = BB%*%t(mixedEffect_list$B_random[[i]])%*%(Y[[i]] - mixedEffect_list$B_fixed[[i]]%*%mixedEffect_list$beta_fixed)
      res <- c(res,Y[[i]] - mixedEffect_list$B_fixed[[i]]%*%mixedEffect_list$beta_fixed - mixedEffect_list$B_random[[i]]%*%br[i,])
    }
    m = br - colMeans(br)
    mixedEffect_list$Sigma = t(m)%*%m/n
    mixedEffect_list$sigma = sqrt(var(res))
  } else {
    for(i in 1:n){
      BB = BB + t(mixedEffect_list$B_fixed[[i]])%*%mixedEffect_list$B_fixed[[i]]
      BY = BY + t(mixedEffect_list$B_fixed[[i]])%*%Y[[i]]
    }
    beta_fixed = solve(BB,BY)
    for(i in 1:n){
      res <- c(res,Y[[i]] - mixedEffect_list$B_fixed[[i]]%*%beta_fixed)
    }
    mixedEffect_list$sigma = sqrt(var(res))
  }
  return(mixedEffect_list)
}

operator.startvalues <- function(Y,locs,mixedEffect_list,operator_list,measurement_list)
{
  if(operator_list$type == "fd2"){
    operator_list$tau = 1/measurement_list$sigma
  } else if(operator_list$type == "matern"){
    operator_list$tau = 1/measurement_list$sigma
    m = min(unlist(lapply(lapply(locs,range),min)))
    M = max(unlist(lapply(lapply(locs,range),max)))
    range = 0.5*(M-m)
    operator_list$kappa = sqrt(8)/range
  }
  return(operator_list)
}
