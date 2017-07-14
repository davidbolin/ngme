#'
#' @title Wrapper for parameter estimation. 
#' 
#' @description A wrapper function for parameter estimation. 
#' 
#' @param Y A numeric list that contains outcome values. 
#' @param locs A numeric list that contains the timings at which the outcomes 
#'    are collected. 
#' @param B_random A numeric list of random effects covariate matrices.
#' @param B_fixed A numeric list of fixed effects covariate matrices.
#' @param operator.type A character string for specifying the operator type. 
#'   Available options are: 
#'   \code{"fd2"} for integrated Random-Walk (integrated Brownian Motion), 
#'   and \code{"matern"} for Matern family. 
#' @param n.process A numerical value for the number of basis functions to
#'     approximate the stochastic process.
#' @param measurement.distribution A character string to specify the distribution 
#'   of the error term. Available options are: \code{"Normal"} for Normal, \code{"NIG"} 
#'   for Normal-inverse Gaussian, \code{"tdist"} for t-distribution.
#' @param random.effect.distribution A character string that indicates the distribution 
#'   of the random effects. Available options are:  \code{"Normal"} for Normal,
#'   and \code{"NIG"} for Normal-inverve Gaussian distributions.
#' @param process.distribution A character vector that indicates the distribution of 
#'   the process. Available options are: 
#'   \code{"Normal"} for Normal, \code{"NIG"} for Normal-inverse Gaussian,
#'   \code{"GAL"} for generalised-asymmetric Laplace, and \code{"CH"} for Cauchy
#'   distributions.
#' @param individual.sigma A logical variable for specifying patient-specific mixture
#'     random variable for the error-term; \code{"FALSE"} indicates do not obtain,
#'     \code{"TRUE"} obtain.
#' @param silent A logical value for printing the details of the iterations;
#'      \code{"TRUE"} indicates do not print, \code{"FALSE"} indicates print.
#' @param estimation.options A list of control inputs. See \code{"estimation.controls"} 
#'     for the \code{"nglda_est"} function. 
#' @param estimate_fisher A logical variable for whether Fisher-Information matrix
#'     to be obtained; \code{"FALSE"} indicates do not obtain, \code{"TRUE"} obtain.
#' @inheritParams  
#' @param ... Additional arguments.     
#' 
#' @return A list of fitted results.
#'
#' @details This function is a wrapper function (wraps \code{"estimateLong"}) 
#'    for parameter estimation. It internally selects the initial values to 
#'    start the stochastic gradient algorithn. The function is not advised to 
#'    be used. It is indeed called within wrapped by \code{"nglda_est"} that is a 
#'    more user-friendly function for parameter estimation. 
#'    
#' @seealso \code{\link{nglda_est}}, \code{\link{estimateLong}}    
#'
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'   estimate.wrapper(...)
#'   }

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
                             subsample.type = 4,
                             nPar_burnin = 0,
                             nIter.fisher = 1000,
                             nSim.fisher = 1000)
  if(!missing(estimation.options) && !is.null(estimation.options)){
    for(i in 1:length(estimation.options)){
      estimation.controls[names(estimation.options)[i]] = estimation.options[i]
    }
}
  if(!silent)
    cat("Setup lists\n")
  Vin <- list()
  n.pers = length(Y)
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
    operator_list$type  <- operator.type
    operator_list <- operator.startvalues(Y,locs,mixedEffect_list,operator_list,measurement_list)


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
                            polyak_rate = estimation.controls$polyak_rate,
                            nSim = estimation.controls$nSim,
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
                          polyak_rate = estimation.controls$polyak_rate,
                          nBurnin = estimation.controls$nBurnin,
                          nIter = estimation.controls$nIter,
                          nPar_burnin = estimation.controls$nPar_burnin,
                          pSubsample = estimation.controls$pSubsample,
                          silent = silent,
                          estimate_fisher = FALSE,
                          ...)
      if(estimate_fisher){
        res.f <- estimateLong(Y, locs,
                          res$mixedEffect_list,
                          res$measurementError_list,
                          res$processes_list,
                          res$operator_list,
                          learning_rate = estimation.controls$learning.rate,
                          nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                          polyak_rate = -1,
                          nBurnin = estimation.controls$nBurnin,
                          nIter = estimation.controls$nIter.fisher,
                          nSim = estimation.controls$nSim.fisher,
                          nPar_burnin = estimation.controls$nPar_burnin,
                          pSubsample = estimation.controls$pSubsample,
                          silent = silent,
                          estimate_fisher = estimate_fisher,
                          ...)
       res$FisherMatrix <- res.f$FisherMatrix
      }
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
                          polyak_rate = estimation.controls$polyak_rate,
                          nSim = estimation.controls$nSim,
                          nBurnin = estimation.controls$nBurnin,
                          nIter = estimation.controls$nIter,
                          nPar_burnin = estimation.controls$nPar_burnin,
                          pSubsample = estimation.controls$pSubsample,
                          subsample.type = estimation.controls$subsample.type,
                          silent = silent,
                          estimate_fisher = FALSE,
                          ...)
      if(estimate_fisher){
        res.f <- estimateLong(Y, locs,
                              res$mixedEffect_list,
                              res$measurementError_list,
                              res$processes_list,
                              res$operator_list,
                              learning_rate = estimation.controls$learning.rate,
                              nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                              polyak_rate = -1,
                              nIter = estimation.controls$nIter.fisher,
                              nSim = estimation.controls$nSim.fisher,
                              nBurnin = estimation.controls$nBurnin,
                              nPar_burnin = estimation.controls$nPar_burnin,
                              pSubsample = estimation.controls$pSubsample,
                              silent = silent,
                              estimate_fisher = estimate_fisher,
                              ...)
        res$FisherMatrix <- res.f$FisherMatrix
      }
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
                            polyak_rate = estimation.controls$polyak_rate,
                            nSim = estimation.controls$nSim,
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
                            polyak_rate = estimation.controls$polyak_rate,
                            nBurnin = estimation.controls$nBurnin,
                            nIter = estimation.controls$nIter,
                            nPar_burnin = estimation.controls$nPar_burnin,
                            pSubsample = estimation.controls$pSubsample,
                            silent = silent,
                            estimate_fisher = FALSE,
                            ...)
        if(estimate_fisher){
          res.f <- estimateLong(Y, locs,
                                res$mixedEffect_list,
                                res$measurementError_list,
                                res$processes_list,
                                res$operator_list,
                                learning_rate = estimation.controls$learning.rate,
                                nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                                polyak_rate = -1,
                                nSim = estimation.controls$nSim.fisher,
                                nBurnin = estimation.controls$nBurnin,
                                nIter = estimation.controls$nIter.fisher,
                                nPar_burnin = estimation.controls$nPar_burnin,
                                pSubsample = estimation.controls$pSubsample,
                                silent = silent,
                                estimate_fisher = estimate_fisher,
                                ...)
          res$FisherMatrix <- res.f$FisherMatrix
        }
      } else {
        res <- estimateLong(Y, locs,
                            mixedEffect_list,
                            measurement_list,
                            learning_rate = estimation.controls$learning.rate,
                            nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                            polyak_rate = estimation.controls$polyak_rate,
                            nSim = estimation.controls$nSim,
                            nBurnin = estimation.controls$nBurnin,
                            nIter = estimation.controls$nIter,
                            nPar_burnin = estimation.controls$nPar_burnin,
                            pSubsample = estimation.controls$pSubsample,
                            silent = silent,
                            estimate_fisher = FALSE,
                            ...)
        if(estimate_fisher){
          res.f <- estimateLong(Y, locs,
                                res$mixedEffect_list,
                                res$measurementError_list,
                                learning_rate = estimation.controls$learning.rate,
                                nBurnin_learningrate = estimation.controls$nBurnin_learningrate,
                                polyak_rate = -1,
                                nSim = estimation.controls$nSim.fisher,
                                nBurnin = estimation.controls$nBurnin,
                                nIter = estimation.controls$nIter.fisher,
                                nPar_burnin = estimation.controls$nPar_burnin,
                                pSubsample = estimation.controls$pSubsample,
                                silent = silent,
                                estimate_fisher = estimate_fisher,
                                ...)
          res$FisherMatrix <- res.f$FisherMatrix
        }
      }
  }

  return(res)
}

#'
#' @title Estimate parameters.
#' 
#' @description A function that estimates parameters by 
#'    calling the \code{"estimateLong_cpp()"} function.
#'
#' @inheritParams estimate.wrapper
#' @param mixedEffect_list A list of inputs for random effects. 
#'   \itemize{
#'   \item \code{noise} The distribution of the mixed effects.
#'   \item \code{B_random} A list that contains the random effect 
#'      covariates (needs to be matrix, can be NULL).
#'   \item \code{B_fixed} A list that contains the fixed effect 
#'      covariates (needs to be matrix, can be NULL).
#'   \item \code{beta_random} Initial values for the parameters of the 
#'      random effects (mean parameter) (if not specified set to zero).
#'   \item \code{beta_fixed} Initial  values for the parameters of the 
#'      fixed effects (if not specified set to zero).
#'   \item \code{Sigma} Initial values for the parameters of the 
#'      variance-covariance matrix of the random effects 
#'      (if not specified set to I ).
#'   \item \code{nu} Shape parameter for noise (NIG only)
#'   \item \code{mu} Shift parameter for noise (NIG only)
#'   \item \code{U} A list of inital values of the random effects.
#'   \item \code{V} A list of inital values of the variance effects.
#'   }
#' @param measurment_list A list of inputs for measurement error.
#'   \itemize{
#'   \item \code{sigma} Measurement noise variance parameter.
#'   \item \code{nu}    Shape parameter for noise (NIG only).
#'   \item \code{Vs}    A list of inital values for the noise of the measurement.
#'   }
#' @param processes_list A list of inputs for the process.
#'   \itemize{
#'   \item \code{noise} Distribution of the process.
#'   \item \code{nu}    Shape parameter (for NIG or GAL).
#'   \item \code{mu}    Asymmetry parameter (for NIG or GAL).
#'   }
#' @param learning_rate A numeric value for the parameter of stochastic gradient.
#' @param nBurnin_learningrate A numeric value until which the learning will
#'     not be started.
#' @param nPar_burnin A numeric value; "M-step" updates will be used until this
#'     iteration.
#' @param polyak_rate A numeric value for moving average of parameters;
#'     -1: inactive, 0: pure mean.
#' @param step0 A numeric value for stepsize for the optimizer; step0 / i^alpha.
#' @param alpha A numeric value for stepsize for the optimizer; step0 / i^alpha.
#' @param pSubsample A numeric value for the portion of data to be used in each
#'     gradient iteration.
#' @param subsample.type A numeric value for the type of subsampling;
#'     0: uniform without sampling, 1: sample size weighted,
#'     3: weighted sampling by gradient size.
#' @param pSubsample2 A numeric value for the portion of the data
#'     to be used in each gradient subsampling weighted by gradient.
#' @param nIter A numeric value for the number of iteration that will be
#'     used by the stochastic gradient.
#' @param nSim A numeric value for the number of samples of the Gibbs sampler
#'     to estimate the gradient.
#' @param silent A logical value for printing the details of the iterations;
#'      \code{"TRUE"} indicates do not print, \code{"FALSE"} indicates print.
#' @param seed A numerical value for starting the Gibbs samplers from fixed seed.
#' 
#' @return A list of fitted results.
#'
#' @details This function calls \code{"estimateLong_cpp()"} internally. 
#'    It is wrapped by \code{"estimate.wrapper"}), and is not advised to 
#'    be used. 
#'    
#' @seealso \code{\link{nglda_est}}, \code{\link{estimate.wrapper}}    
#'
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'   estimateLong(...)
#'   }

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
    obs_list[[i]] <- list(Y=Y[[i]], locs = locs[[i]])
    if(use.process){
      obs_list[[i]]$A = build.A.matrix(operator_list,locs,i)
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
  free.samples = 0
  groups <- list()
  if(subsample.type == 4){
    group <- group.fixed(mixedEffect_list$B_fixed)
    groups <- group$groups
    free.samples = group$free
  }

  input <- list( obs_list         = obs_list,
                 measurementError_list  = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 pSubsample       = pSubsample,
                 pSubsample2      = pSubsample2,
                 free_samples     = free.samples,
                 group_list       = groups,
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


  return(output)
}

#'
#' @title estimating mixed effect model
#'
#' @description
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
#' @param silent          - print iteration info
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
    res.list = list()
    for(i in 1:n){
      BB = ginv(t(mixedEffect_list$B_random[[i]])%*%mixedEffect_list$B_random[[i]])
      br[i,] = BB%*%t(mixedEffect_list$B_random[[i]])%*%(Y[[i]] - mixedEffect_list$B_fixed[[i]]%*%mixedEffect_list$beta_fixed)
      res.list[[i]] =  Y[[i]] - mixedEffect_list$B_fixed[[i]]%*%mixedEffect_list$beta_fixed - mixedEffect_list$B_random[[i]]%*%br[i,]
      res <- c(res,res.list[[i]])
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
    res.list = list()
    for(i in 1:n){
      res.list[[i]] = Y[[i]] - mixedEffect_list$B_fixed[[i]]%*%beta_fixed
      res <- c(res,res.list[[i]])
    }
    mixedEffect_list$sigma = sqrt(var(res))
  }
  #mixedEffect_list$res = res.list
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
    range = min(4*operator_list$h[[1]][1],0.1*(M-m))
    operator_list$kappa = sqrt(8)/range
  }
  return(operator_list)
}
