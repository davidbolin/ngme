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
#' @param processes_list  - for the stochastic noise driving the
#' @param noise           - either Normal, NIG or GAL (change name to type rather then noise)
#' @param nu              - shape parameter for NIG or GAL
#' @param mu              - assymetric parameter for NIG or GAL
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
estimateLong <- function(Y,
                         locs,
                         mixedEffect_list,
                         measurment_list,
                         processes_list,
                         operator_list,
                         step0 = 0.5,
                         alpha = 0.1,
                         learning_rate = 0,
                         pSubsample = 1.,
                         polyak_rate = -1.,
                         subsample.type = 1,
                         nIter = 10,     # iterations to run the stochastic gradient
                         nSim  = 1,
                         nBurnin = 10,   # steps before starting gradient estimation
                         silent  = FALSE # print iteration info
                         )
{
  obs_list <- list()
  common.grid = FALSE
  if(length(operator_list$loc)==1){
    common.grid = TRUE
  }
  for(i in 1:length(locs)){
    if(length(Y[[i]]) != length(locs[[i]])){
      stop("Length of Y and locs differ.")
    }
    if(common.grid){
      obs_list[[i]] <- list(A = spde.A(locs[[i]],
                                       operator_list$loc[[1]],
                                       right.boundary = operator_list$right.boundary,
                                       left.boundary = operator_list$left.boundary),
                            Y=Y[[i]],
                            locs = locs[[i]])
    } else {
      obs_list[[i]] <- list(A = spde.A(locs[[i]],
                                       operator_list$loc[[i]],
                                       right.boundary = operator_list$right.boundary,
                                       left.boundary = operator_list$left.boundary),
                            Y=Y[[i]],
                            locs = locs[[i]])
    }
  }


  input <- list( obs_list         = obs_list,
                 operator_list    = operator_list,
                 measurementError_list  = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 processes_list   = processes_list,
                 pSubsample       = pSubsample,
                 subsample_type   = subsample.type,
                 nIter            = nIter,     # iterations to run the stochastic gradient
                 nSim             = nSim,
                 nBurnin          = nBurnin,   # steps before starting gradient estimation
                 silent           = silent, # print iteration info)
                 step0            = step0,
                 alpha            = alpha,
                 common.grid      = common.grid,
                 learning_rate    = learning_rate,
                 polyak_rate      = polyak_rate
              )

  output <- estimateLong_cpp(input)

  output$operator_list$left.boundary <- operator_list$left.boundary
  output$operator_list$right.boundary <- operator_list$right.boundary
  output$operator_list$type <- operator_list$type
  output$operator_list$Q <- operator_list$Q

  output$measurementError_list$Vs <- measurment_list$Vs
  output$measurementError_list$common_V <- measurment_list$common_V


    return(output)
}


