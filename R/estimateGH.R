library(GHmixedeffect)

#' @param   Y           - list with the observations
#' @param   locs        - list with position of the observations (Y)
#' @param mixedEffect_list -
#' @param   meas_noise  - the aviable noise classes: Normal or NIG
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
#' @param processes_list  - for the stochastic noise driving the
#' @param noise           - either Normal, NIG or GAL (change name to type rather then noise)
#' @param nu              - shape parameter for NIG or GAL
#' @param mu              - assymetric parameter for NIG or GAL
#'
#' @param step0           - stepsize for optimizer is step0 / i^alpha
#' @param alpha           - stepsize for optimizer is step0 / i^alpha
#' @param pSubsample      - precentage of data used in each gradient subsampling
#' @param nIter           - number of iteration of the stochastic gradient
#' @param nSim            - number of samples of the gibbs sampler to estimate the gradient
#' @parma silent          - print iteration info
estimateLong <- function(Y,
                         locs,
                         mixedEffect_list,
                         measurment_list,
                         processes_list,
                         operator_list,
                         step0 = 0.5,
                         alpha = 0.1,
                         pSubsample = 1.,
                         nIter = 10,     # iterations to run the stochastic gradient
                         nSim  = 1,
                         nBurnin = 10,   # steps before starting gradient estimation
                         silent  = FALSE # print iteration info
                         )
{
  obs_list <- list()
  for(i in 1:length(locs))
    obs_list[[i]] <- list(A = spde.A(locs[[i]],
                                     operator_list$loc,
                                     right.boundary = operator_list$right.boundary,
                                     left.boundary = operator_list$left.boundary),
                          Y=Y[[i]],
                          locs = locs[[i]])

  input <- list( obs_list         = obs_list,
                 operator_list    = operator_list,
                 measurementError_list  = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 processes_list   = processes_list,
                 pSubsample       = pSubsample,
                 nIter            = nIter,     # iterations to run the stochastic gradient
                 nSim             = nSim,
                 nBurnin          = nBurnin,   # steps before starting gradient estimation
                 silent           = silent, # print iteration info)
                 step0            = step0,
                 alpha            = alpha
              )

  output <- estimateLong_cpp(input)
  return(output)
}


#' Estimating longitudal model
#'
#' @param Y list of observations
#' @param locs list of location of observations
#' @param theta list with  kappa, tau
#' @param stepsize length to step in gradient dir
#' @param n discretization size
#' @param Niter  - max number of iterations
#' @param burnin - sampling prior to stating calculating the gradient
#' @param long.percentage percentage of longitudinal samples to use in each iteration.
#' @param nsim - how many samples to calculate the gradient
#' @param B list of matrix with covariates
#' @param noise - either Normal, Generalized assymetric Laplace, or Normal inverse Gaussian
#' @param silent - print redsults
#' @param V - list of inital guess of Variance components
#' @param U - initaul guess of random intercept
#' @param mixedEffectList - list contaning $B (mixed effects)
#' @param Bmixed - the covariates for mixed effects (if mixeEffectList is not null, otherwise uses only mixed effects)
#' @param operatorType - type of operator: 1) Matern 2) Finite Difference
#'
#'
estimateLongGH <- function(Y,
                           loc,
                           Bfixed       = NULL,
                           Brandom      = NULL,
                           theta = NULL,
                           stepsize = 1e-2,
                           n = 100,
                           Niter = 100,
                           burnin = 10,
                           long.percentage = 100,
                           commonsigma = TRUE, # only implimented right know
                           nsim = 1,
                           noise  = c("Normal","NIG","GAL"),
                           mNoise = c("Normal", "NIG"),
                           mNoiseList = NULL,
                           silent=FALSE,
                           V = NULL,
                           mixedEffect = FALSE,
                           mixedType   = c("Normal", "NIG"),
                           mixedEffectList = NULL,
                           operatorType = "Matern"
                           )
{
  noise <- match.arg(noise)
  mixedType <- match.arg(mixedType)
  mixedType <- match.arg(mixedType)
  mNoise    <- match.arg(mNoise)
  if(missing(Y))
    stop("Must supply data")

  if(missing(loc))
    stop("Must supply list with locations")

  if(is.null(mNoiseList))
    mNoiseList     <- MeasurementErrorInit(Y, mNoise)

  if(is.null(mixedEffectList))
    mixedEffectList <- MixedInit(mixedType, B_fixed = Bfixed, B_random = Brandom)

  operator_List <- create_operator(locs, n, name = operatorType)
  mesh1d <- operator_List$mesh1d

  obs_ <- list()
  for(i in 1:length(locs))
  {
    obs_[[i]] <- list(A = inla.mesh.1d.A(mesh1d, locs[[i]]), Y=Y[[i]], locs = locs[[i]])
  }
  Nlong = length(locs)
  Nlong = max(min(round(long.percentage*Nlong/100),Nlong),1)

  if(is.null(V))
  {
    if(noise == "Normal"){
      for(i in 1:length(locs))
      {
        V[[i]] <- operator_List$h
      }
    }else{
      for(i in 1:length(locs))
      {
        V[[i]] <- c(rep(1, length(operator_List$h)))
      }
    }

  }

  result$mixedEffect_list <- mixedEffectList
  result$MeasureError_list      <- mNoiseList

  result$obs_list            <- obs_
  result$operator_list   <- operator_list

  result$stepsize <- stepsize
  result$Niter    <- Niter
  result$nsim     <- nsim
  result$Nlong    <- Nlong
  result$nsim     <- nsim
  result$burnin   <- burnin
  result$silent   <- silent
  result$processes_list <- list()
  result$processes_list$noise <- noise
  result$processes_list$tau   <- 1.
  result$processes_list$kappa <- 1.
  if(noise == 'NIG')
  {
    result$processes_list$tau <- 10.
    result$processes_list$mu  <- 0.
  }
  result$processes_list$V_list     <- V

  output <- estimateLong_cpp(result)

  return(output)
}

#' posterior samples generates mean and variance of X, and V
#'
#' @param output estimateLongGH
#'
sample_posteriror<- function(output, sim = 100)
{
  noise.i = 0
  if(output$noise == "Normal"){
    noise.i = -1
  }else if(output$noise == "GAL"){
    noise.i = 1
  }
  sim_res <- samplePosteriorGH(output$obs_, output$operator_List, output$theta,output$mixedEffectList, output$theta$V, sim, noise.i,output$commonsigma)
  return(sim_res)
}


#' plot parameter tracjetories
#'
#' @param output estimateLongGH
#'
plotLongGh <-function(output, beta = c(1))
{
  theta <- output$theta
  par(mfrow=c(3,2))
  plot(theta$beta_vec[,beta], main=expression(beta))
  if(output$operator_List$type == 'matern')
    plot(theta$kappa_vec, main = expression(kappa))
  if(output$commonsigma == 1){
    plot(theta$sigma_vec, main = expression(sigma))
  }else{
    plot(theta$asigma_vec, main = expression(asigma))
    plot(theta$bsigma_vec, main = expression(bsigma))
  }

  plot(theta$tau_vec,   main = expression(tau))

  if(output$noise != "Normal")
    plot(theta$lambda_vec,   main = expression(lambda))

}

plottraj <- function(nr, output, samples)
{
  X_mean <- samples$X_mean[[nr]]
  plot(output$operator_List$mesh1d$loc, X_mean, xlab = 't', ylab = 'X(t)', type = 'l')

}
#' simualte from the prior model using output from parameter estimate
#' @param output estimateLongGH
sample_prior_from_result_LongGH <- function(output)
{
  theta_in <- output$theta
  theta_in$GAL <- result$GAL
  return(simulateLongGH_cpp(output$obs_, output$operator_List, theta_in))
}