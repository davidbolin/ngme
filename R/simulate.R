# simulate data from the prior model
# @param Y only used to get size of objects
simulateLongPrior <- function( Y,
                               locs,
                               mixedEffect_list,
                               measurment_list,
                               processes_list,
                               operator_list)
{
  common.grid = FALSE
  if(length(operator_list$loc)==1){
    common.grid = TRUE
  }
  obs_list <- list()
  for(i in 1:length(locs)){
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



  input <- list( obs_list = obs_list,
                 operator_list = operator_list,
                 measurment_list = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 processes_list = processes_list)

  output <- simulateLongGH_cpp(input)
  return(output)
}


#' Simulating longitudal model
#'
#' @param locs list of location of observations
#' @param theta lost with covariates mu, kappa, sigma_eps, sigma
#' @param B list of matrix with covariates
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

  if(operatorType=="Matern"){
    kappa = theta$kappa
    K = sqrt(tau) * (operator_List$G + kappa*operator_List$C)
    Q = (K%*%operator_List$Ci%*%K)
    R = chol(Q)
  }else{
    K = tau*operator_List$Q[[1]]

    Ci = as(sparseMatrix(i=1:n,j=1:n,x=1/operator_List$h[[1]],dims=c(n, n)), "CsparseMatrix")
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
      V[[i]] =  rGIG(rep(-0.5, n),
                     rep( theta$nu, n),
                     h^2 * theta$nu)
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
