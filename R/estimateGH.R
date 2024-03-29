#' @title Estimate parameters.
#'
#' @description A function that estimates parameters by
#'    calling the \code{"estimateLong_cpp()"} function.
#'
#' @param Y A numeric list that contains outcome values.
#' @param locs A numeric list that contains the timings at which the outcomes
#'    are collected.
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
#' @param step0 A numeric value for stepsize for the optimizer; step0 / i^alpha.
#' @param alpha A numeric value for stepsize for the optimizer; step0 / i^alpha.
#' @param learning_rate A numeric value for the parameter of stochastic gradient.
#' @param nBurnin_learningrate A numeric value until which the learning will
#'     not be started.
#' @param nBurnin_base A numerical value for burn-in simulations that are performed
#'       for a subject that is sampled for the first time in the estimation method.
#' @param pSubsample A numeric value for the portion of data to be used in each
#'     gradient iteration.
#' @param polyak_rate A numeric value for moving average of parameters;
#'     -1: inactive, 0: pure mean.
#' @param subsample.type A numeric value for the type of subsampling;
#'       1: uniform sampling,
#'       2: sample size weighted,
#'       3: weighted sampling by gradient size,
#'       4: grouped sub-sampler.
#' @param nPar_burnin A numeric value; "M-step" updates will be used until this
#'     iteration.
#' @param pSubsample2 A numeric value for the portion of the data
#'     to be used in each gradient subsampling weighted by gradient.
#' @param nIter A numeric value for the number of iteration that will be
#'     used by the stochastic gradient.
#' @param nSim A numeric value for the number of samples of the Gibbs sampler
#'     to estimate the gradient.
#' @param silent A logical value for printing the details of the iterations;
#'      \code{"TRUE"} indicates do not print, \code{"FALSE"} indicates print.
#' @param seed A numerical value for starting the Gibbs samplers from fixed seed.
#' @param standardize.mixedEffects A logical variable for standardising the covariates;
#'       \code{"FALSE"} indicates no standardisation, \code{"TRUE"} standardisation.
#' @param estimate_fisher A logical variable for whether Fisher-Information matrix
#'     to be obtained; \code{"FALSE"} indicates do not obtain, \code{"TRUE"} obtain.
#' @return A list of output.
#'
#' @details This function calls \code{"estimateLong_cpp()"} internally.
#'    It is wrapped by \code{"ngme"}, and is not advised to
#'    be used alone.
#'
#' @seealso \code{\link{ngme}}
#'
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'   fit <- estimateLong(...)
#'   }

estimateLong <- function(Y,
                         locs = NULL,
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
                         iter_start = 0, #starting iteration to use when continuing estimation
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
  bivariate <- FALSE
  if(use.process){
    if(operator_list$type == "matern bivariate")
      bivariate <- TRUE
      if(processes_list$noise == "Normal"){
        operator_list$estimate_theta = 0 #theta not identifiable if the model is Gaussian
      } else {
        operator_list$estimate_theta = 1
      }
  }
  for(i in 1:length(Y)){
    if(use.process){
      A <- build.A.matrix(operator_list,locs,i)
      Yi <- Y[[i]]
      locsi = locs[[i]]
      if(bivariate){ #for bivariate fields, stack observations 
        A1 <- A[!is.na(Yi[,1]),]
        A2 <- A[!is.na(Yi[,2]),]
        A <- bdiag(A1, A2)
        Yi <- c(Y[[i]])
        Yi <- Yi[!is.na(Yi)]
        locs1 = locs[[i]][!is.na(Y[[i]][,1]),]
        locs2 = locs[[i]][!is.na(Y[[i]][,2]),]
        locsi <- rbind(locs1,locs2)
      }
      obs_list[[i]] <- list(A = A, Y=Yi, locs = locsi)
    } else {
      obs_list[[i]] <- list(Y=Y[[i]])  
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
                 iter_start       = iter_start,
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
  }

  if(is.null(nBurnin_learningrate) == FALSE)
    input$nBurnin_learningrate =  nBurnin_learningrate
  if(is.null(seed) == FALSE)
    input <- setseed_ME(input, seed)

  #input <- asS4(input)
  output <- estimateLong_cpp(input)

  if(use.process){
    output$A <- lapply(1:length(obs_list), function(i) obs_list[[i]]$A)
  }
    
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


#' @title STUFF.
#'
#' @description STUFF
#' @param input STUFF
#' @param seed STUFF
#'
#' @return STUFF
#'
#' @details STUFF
#'
#' @seealso \code{\link{estimateLong}}
#'
#' @examples
#'   \dontrun{
#'   setseed_ME(...)
#'   }
#'

setseed_ME <- function(input, seed)
{
  seed.old <- sample.int(10^6, 1)
  set.seed(seed)
  input$seed                  <- sample.int(10^6, 1)
  input$mixedEffect_list$seed <- sample.int(10^6, 1)
  set.seed(seed.old)
  return(input)
}


#'
#' @title Obtain initials for random effects.
#'
#' @description A function to obtain initial values for the random effects.
#'
#' @inheritParams estimateLong
#' @param mixedEffect_list A list for random effects.
#'
#' @return Returns a list for the mixed effects.
#'     See e.g. \code{mixedEffect_list} in the \code{"EstimateLong"} function.
#'
#' @details The function treats random effects as fixed effects to
#'    obtain start values using ordinary least squares (OLS).
#'
#' @seealso \code{\link{ngme}}
#'
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'   ME.startvalues(...)
#'   }
#'
ME.startvalues <- function(Y, mixedEffect_list)
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
  res = NULL
  if(nc.r  > 0){
    beta_r = matrix(0,n,nc.r)
    I = diag(1,nc)

    for(i in 1:n){

      Bi = cbind(mixedEffect_list$B_fixed[[i]],mixedEffect_list$B_random[[i]])
      BB = BB + t(Bi)%*%Bi
      BY = BY + t(Bi)%*%Y[[i]]
    }
    beta = solve(BB,BY)
    mixedEffect_list$beta_fixed = beta[1:nc.f]
    mixedEffect_list$beta_random = beta[(nc.f+1):nc]
    
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
    mixedEffect_list$beta_fixed = solve(BB,BY)
    res.list = list()
    for(i in 1:n){
      res.list[[i]] = Y[[i]] - mixedEffect_list$B_fixed[[i]]%*%mixedEffect_list$beta_fixed
      res <- c(res,res.list[[i]])
    }
    mixedEffect_list$sigma = sqrt(var(res))
  }
  #mixedEffect_list$res = res.list
  return(mixedEffect_list)
}

#'
#' @title Obtain initials for random effects.
#'
#' @description A function to obtain initial values for random effects for bivariate models.
#'
#' @inheritParams estimateLong
#' @param mixedEffect_list A list for random effects.
#'
#' @return Returns a list for the mixed effects.
#'     See e.g. \code{mixedEffect_list} in the \code{"EstimateLong"} function.
#'
#' @details The function treats random effects as fixed effects to
#'    obtain start values using ordinary least squares (OLS).
#'
#' @seealso \code{\link{ngme}}
#'
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'   ME.startvalues(...)
#'   }
#'

ME.startvalues.bivariate <- function(Y, mixedEffect_list)
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
  res = NULL
  if(nc.r  > 0){
    beta_r = matrix(0,n,nc.r)
    I = diag(1,nc)
    
    for(i in 1:n){
      cY = cbind(Y[[i]][,1],Y[[i]][,2])
      ind <- !is.na(cY)
      Bi = cbind(mixedEffect_list$B_fixed[[i]],mixedEffect_list$B_random[[i]])
      Bi = Bi[ind,]
      BB = BB + t(Bi)%*%Bi
      BY = BY + t(Bi)%*%cY[ind]
    }
    beta = solve(BB,BY)
    mixedEffect_list$beta_fixed = beta[1:nc.f]
    mixedEffect_list$beta_random = beta[(nc.f+1):nc]
    
    Sigma = matrix(0,nc.r,nc.r)
    br = matrix(0,n,nc.r)
    res.list = list()
    for(i in 1:n){
      cY = cbind(Y[[i]][,1],Y[[i]][,2])
      ind <- !is.na(cY)
      BB = ginv(t(mixedEffect_list$B_random[[i]])%*%mixedEffect_list$B_random[[i]])
      br[i,] = BB%*%t(mixedEffect_list$B_random[[i]][ind,])%*%(cY[ind] - mixedEffect_list$B_fixed[[i]]%*%mixedEffect_list$beta_fixed)
      res.list[[i]] =  cY[ind] - mixedEffect_list$B_fixed[[i]][ind,]%*%mixedEffect_list$beta_fixed - mixedEffect_list$B_random[[i]][ind,]%*%br[i,]
      res <- c(res,res.list[[i]])
    }
    m = br - colMeans(br)
    mixedEffect_list$Sigma = t(m)%*%m/n
    mixedEffect_list$sigma = sqrt(var(res))
  } else {
    for(i in 1:n){
      cY = c(Y[[i]][,1],Y[[i]][,2])
      ind <- !is.na(cY)
      BB = BB + t(mixedEffect_list$B_fixed[[i]])%*%mixedEffect_list$B_fixed[[i]]
      BY = BY + t(mixedEffect_list$B_fixed[[i]])%*%cY[ind]
    }
    mixedEffect_list$beta_fixed = solve(BB,BY)
    res.1 = c()
    res.2 <- c()
    for(i in 1:n){
      cY = c(Y[[i]][,1],Y[[i]][,2])
      Y_ <- cbind(Y[[i]][,1],Y[[i]][,2])
      ind.2 <- !is.na(Y_)
      ind <- !is.na(cY)
      res_ <- 0*Y_
      res_[ind.2]   <- cY[ind] - mixedEffect_list$B_fixed[[i]]%*%mixedEffect_list$beta_fixed
      res.1  = c(res.1, res_[!is.na(Y_[,1]), 1])
      res.2  = c(res.2, res_[!is.na(Y_[,2]), 2])
    }
    mixedEffect_list$theta = log(c(sd(res.1), sd(res.2)))
  }
  #mixedEffect_list$res = res.list
  return(mixedEffect_list)
}



#' @title Obtain initial values for the operator.
#'
#' @description A function to obtain initials for the operator.
#'
#' @inheritParams estimateLong
#' @param mixedEffect_list A list for random effects.
#' @param operator_list A list for operator.
#' @param measurement_list A list for measurement error.
#'
#' @return A list for operator.
#'    See e.g. \code{"operator_list"} in \code{"estimateLong"}.
#'
#' @examples
#'   \dontrun{
#'   operator.startvalues(...)
#'   }
#'
operator.startvalues <- function(Y, locs, mixedEffect_list, operator_list, measurement_list)
{
  if(operator_list$type == "fd2"){
    operator_list$tau = 1/measurement_list$sigma
  } else if(operator_list$type == "matern" || operator_list$type == "exponential"){
    operator_list$tau = 1/measurement_list$sigma
    m = min(unlist(lapply(lapply(locs,range),min)))
    M = max(unlist(lapply(lapply(locs,range),max)))
    h.min = min(unlist(lapply(operator_list$h,min)))
    range = max(4*h.min,0.1*(M-m))
    if(operator_list$type == "matern"){
      operator_list$kappa = sqrt(8*3/2)/(0.2*range  )
    } else {
      operator_list$kappa = sqrt(8*1/2)/(0.2 *range)
    }
    
  }
  return(operator_list)
}

#' @title Obtain initial values for the operator.
#'
#' @description A function to obtain initials for the bivariate operator.
#'
#' @inheritParams estimateLong
#' @param mixedEffect_list A list for random effects.
#' @param operator_list A list for operator.
#' @param measurement_list A list for measurement error.
#'
#' @return A list for operator.
#'    See e.g. \code{"operator_list"} in \code{"estimateLong"}.
#'
#' @examples
#'   \dontrun{
#'   operator.startvalues(...)
#'   }
#'
operator.startvalues.bivariate <- function(Y, locs, mixedEffect_list, operator_list, measurement_list)
{
  operator_list$tau1 = 1/exp(measurement_list$theta[1])
  operator_list$tau2 = 1/exp(measurement_list$theta[2])
  
  m = min(unlist(lapply(lapply(locs,range),min)))
  M = max(unlist(lapply(lapply(locs,range),max)))
  h.min = min(unlist(lapply(operator_list$h,min)))
  range = max(4*h.min,0.1*(M-m))
  
  operator_list$kappa1 = sqrt(8*3/2)/(0.3*range)  
  operator_list$kappa2 = sqrt(8*3/2)/(0.3*range) 
  operator_list$rho = 0
  operator_list$theta = 0
    
  return(operator_list)
}
