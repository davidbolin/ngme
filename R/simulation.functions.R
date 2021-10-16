simulate.ngme <- function(object,id=NULL)
{
  
  if(!is.null(id)){
    id_list <- as.numeric(names(object$Y))
    ind <- which(id_list %in% id)
    object <- crop.lists(object,ind)
    
  }
  
  sim_object <- simulateLongPrior(Y                 = object$Y,
                                  locs              = object$locs,
                                  mixedEffect_list  = object$mixedEffect_list,
                                  measurment_list   = object$measurementError_list,
                                  processes_list    = object$processes_list,
                                  operator_list     = object$operator_list)
  
  a <- attributes(object$Y)
  for(i in 1:length(sim_object)){
    if(names(sim_object)[i]!="U"){
      attributes(sim_object[[names(sim_object)[i]]]) <- a
    }
  }
  
  return(sim_object)
}
#' @title Simulates prior model using processes and operator
#' 
#' @description 
#' 
#' @param n number of simulations
#' @param process_list procsess list
#' @param operator_list operator list
#' @return A list of simulation
simulate.process <- function(n, operator_list, process_list){
  
  process_list$V <- list()
  process_list$X <- list()
  for(i in 1:length(operator_list$h)){
    process_list$X[[i]] <- rep(0, length(operator_list$h[[i]]))
    process_list$V[[i]] <- operator_list$h[[i]]
  }
  if(grepl('bivariate',operator_list$type)){
    process_list$Bmu <- list()
    process_list$Bnu <- list()
    n.grid <- length(operator_list$h[[1]])/2
    for(i in 1:length(operator_list$h)){
      process_list$Bmu[[i]] = kronecker(diag(2),matrix(rep(1, n.grid)))
      process_list$Bnu[[i]] = kronecker(diag(2),matrix(rep(1, n.grid)))
    }
  }
  simulateLongProcesses_cpp(list(nsim = n,
                                 operator_list=operator_list,
                                 processes_list=process_list))
}


#' @title Simulates nig distribution
#'
#' @description internal function for simulating nig iid random variables using
#'              V ~ IG(nu,h^2 nu)
#'              Z ~ N(0, 1)
#'              X ~ delta + (V - h) * mu + sqrt(V) * sigma  * Z
#'               
#' 
#' @param n number of simulations (if h is not NULL then n is ignored)
#' @param delta (real) location parameter
#' @param mu    (real) symmetric parameter
#' @param sigma (real) scale parameter
#' @param nu    (real, >0) shape parameter
#' @param h     (m x 1) (h = 1) discritization vector
#' @export
rNIG <- function(n = 1, delta, mu, sigma, nu, h = NULL){
  if(is.null(h))
    h = rep(1,n)
  n = length(h)
  if(length(nu) != n)
    nu = rep(nu, n)
  if(length(delta)!=n)
    delta = rep(delta, n)
  if(length(mu)!=n)
    mu = rep(mu, n)
  if(length(sigma)!=n)
    sigma = rep(sigma, n)
  
  V =  ngme::rGIG(p = rep(-0.5,n) , nu, h^2 * nu, sample.int(10^6, 1))
  X = h*delta + (V - h) * mu + sqrt(V) * sigma * rnorm(n)
  return(X)
}


#' @title Simulates data from the prior model.
#' 
#' @description 
#' 
#' @param Y only used to get size of objects
#' @param locs measurement locations
#' @param mixedEffect_list mixed effects list
#' @param measurement_list measurement error list
#' @param process_list procsess list
#' @param operator_list operator list
#' @details STUFF 
#' @return A list of outputs.
#' @examples
#'   \dontrun{
#'   simulateLongPrior(...)
#'   }
simulateLongPrior <- function( Y,
                               locs,
                               mixedEffect_list,
                               measurment_list,
                               processes_list,
                               operator_list)
{
  bivariate <- FALSE
  if(!missing(processes_list) && !is.null(processes_list)){
    if(operator_list$type == "matern bivariate")
      bivariate <- TRUE
  }
  
  if(missing(Y)){
    Y <- list()
    for(i in 1:length(locs)){
      if(is.matrix(locs[[i]])){
        Y[[i]] <- rep(1,dim(locs[[i]])[1])
      } else {
        Y[[i]] <- rep(1,length(locs[[i]])[1])
      }
      if(bivariate){
        Y[[i]] <- cbind(Y[[i]],Y[[i]])
      }
    }
  }

  if(!missing(processes_list) && !is.null(processes_list)){
    common.grid = FALSE
    obs_list <- list()
    for(i in 1:length(locs)){
      A <- build.A.matrix(operator_list,locs,i)
      Yi <- Y[[i]]
      locsi = locs[[i]]
      if(bivariate){ #for bivariate fields, stack observations 
        #do not and remove NaN, simulate all locations
        #A1 <- A[!is.nan(Y[[i]][,1]),]
        #A2 <- A[!is.nan(Y[[i]][,2]),]
        A <- bdiag(A,A)
        Yi <- c(Y[[i]])
        #Yi <- Yi[!is.nan(Yi)]
        #locs1 = locs[[i]][!is.nan(Y[[i]][,1]),]
        #locs2 = locs[[i]][!is.nan(Y[[i]][,2]),]
        locsi <- rbind(locsi,locsi)
        
      }
      obs_list[[i]] <- list(A = A, Y=Yi, locs = locsi)
    }

    input <- list( obs_list = obs_list,
                   operator_list = operator_list,
                   measurment_list = measurment_list,
                   mixedEffect_list = mixedEffect_list,
                   processes_list = processes_list)
    output <- simulateLongGH_cpp(input)
  } else {
    obs_list <- list()
    for(i in 1:length(locs)){
      if(length(Y[[i]]) != length(locs[[i]])){
        stop("Length of Y and locs differ.")
      }
      obs_list[[i]] <- list(Y=Y[[i]], locs = locs[[i]])
    }

    input <- list( obs_list = obs_list,
                   measurment_list = measurment_list,
                   mixedEffect_list = mixedEffect_list)

    output <- simulateLongME_cpp(input)
    }
    output$Ystar <- lapply(1:length(output$Y),function(i) output$Y[[i]]-output$E[[i]])
    if(bivariate){
      for(i in 1:length(locs)){
        #reshape Y and remove NaN
        Yi <- output$Y[[i]] 
        Yi.star <- output$Ystar[[i]] 
        n.i <- length(Yi)
        Yi <- matrix(Yi,n.i/2,2)
        Yi.star <- matrix(Yi.star,n.i/2,2)
        Yi[is.nan(Y[[i]])] <- NaN
        Yi.star[is.nan(Y[[i]])] <- NaN
        output$Y[[i]] <- Yi
        output$Y.star[[i]] <- Yi.star
      }
    }

  return(output)
}

#' @title Simulating longitudal model using only R functions
#' 
#' @description STUFF
#'
#' @param locs list of observation locations
#' @param theta list with parameters covariates mu, kappa, sigma_eps, sigma
#' @param B list of matrix with covariates
#' @return A list of outputs.
#' @examples
#'   \dontrun{
#'   simulateLong.R(...)
#'   }
simulateLong.R <- function(loc,
                           theta,
                           n = 100,
                           B = NULL,
                           noise = c("Normal","NIG","GAL"),
                           mixedEffect = FALSE,
                           Bmixed      = NULL,
                           mixedEffectList = NULL,
                           operatorType = "Matern",
                           boundary = NULL)
{
  if(is.list(locs)){
    nrep = length(loc)
  } else {
    nrep = 1
    locs = loc
    loc = list()
    loc[[1]] = locs
  }

  if(is.null(B))
  {
    B=list()
    for(i in 1:nrep)
      B[[i]] <- as.matrix(rep(1, length(loc[[i]])))
    theta$beta <- c(0)
  }

  if(nrep != length(B))
  {
    cat('locs and B should be equal length lists\n')
    return(-1)
  }


  operator_List <- create_operator(loc, n, name = operatorType)


  tau = as.double(theta$tau)
  sigma = theta$sigma
  beta = as.double(theta$beta)
  if(length(sigma) == 1){
    sigma <- rep(sigma,nrep)
  }

  if(operatorType=="matern"){
    kappa = theta$kappa
    K = tau * (operator_List$G[[1]] + kappa*operator_List$C[[1]])
    Q = (K%*%operator_List$Ci[[1]]%*%K)
    R = chol(Q)
  }else{
    K = tau*operator_List$Q[[1]]
    Ci = as(sparseMatrix(i=1:n,j=1:n,x=1/operator_List$h[[1]],dims=c(n, n)), "CsparseMatrix")
    print(dim(Ci))
    Q = Matrix::t(K) %*% Ci %*% K
    R = chol(Q)
  }

  A <- list()
  for(i in 1:length(loc))
  {
    A[[i]] <-  spde.A(x = operator_List$loc[[1]], loc = loc[[i]])
  }

  Y = list()
  X = list()
  Zs = list()
  V = list()
  h <- operator_List$h[[1]]
  for(i in 1:nrep){
    if(noise == "Normal"){
      X[[i]] = solve(R,rnorm(dim(R)[1]))
    }else if (noise == "NIG"){
      V[[i]] =  LDMod::rGIG(rep(-0.5, n),
                     rep( theta$nu, n),
                     h^2 * theta$nu, sample.int(10^6, 1))
      Z <- (- h  + V[[i]]) * theta$mu + sqrt(V[[i]]) * rnorm(n)
      X[[i]] <- solve(K, Z)
      Zs[[i]] <- Z
    }else if( noise == "GAL"){
      V[[i]] =  rgamma(n, h * theta$nu, rep(theta$nu, n)) + 10e-14
      Z <- (- h  + V[[i]]) * theta$mu + sqrt(V[[i]]) * rnorm(n)
      Zs[[i]] <- Z
      X[[i]] <- solve(K, Z)
    }else if( noise == "CH"){
      V[[i]] =  1/rgamma(n, 0.5, 0.25 * h^2)
      Z <-  sqrt(V[[i]]) * rnorm(n)
      Zs[[i]] <- Z
      X[[i]] <- solve(K, Z)
    }

    Y[[i]] = (B[[i]]%*%beta + A[[i]]%*%X[[i]] + sigma[i]*rnorm(dim(A[[i]])[1]))@x
  }
  return(list(Y=Y, X=X, xloc = operator_List$mesh1d$loc, A=A, V= V, Z = Zs))
}
