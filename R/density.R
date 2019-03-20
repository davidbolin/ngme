#' @title Density function of Normal inverse Gaussian distribution.
#'
#' @description A function to calculate the value of the probability density 
#'    function of Normal inverse Gaussian distribution.
#' @param x     A numeric vector for quantiles.
#' @param delta A numeric value for the location parameter.
#' @param mu    A numeric value for the shift parameter.
#' @param nu    A numeric value for the shape parameter.
#' @param sigma A numeric value for the scaling parameter.
#' @details \eqn{f(x|\delta, \mu, \nu, \sigma) =} STUFF 
#' @return A list of outputs.
#' @examples
#'   \dontrun{
#'   dnig(...)
#'   }

dnig <- function(x, delta, mu, nu, sigma,log=F)
{
  c0 <- sqrt(nu) * sqrt( mu^2/sigma^2 + nu) / pi
  f <- nu + mu * (x - delta) /sigma^2
  coeff <- sqrt( nu * sigma^2 + (x - delta)^2)
  f <- f + log(c0) -  log(coeff)
  f <- f + log(besselK(coeff  * sqrt(mu^2/sigma^4 + nu/sigma^2), -1,TRUE)) -coeff  * sqrt(mu^2/sigma^4 + nu/sigma^2)
  if(log==F)
    return(exp(f))
  
  return(f)
}

#dgig <-function(x, p, a, b, log=T){
#  f <- 0.5 * p * log(a/b) +(p-1)*log(x) - a*x/2 - b*x/2
#  f <- f - log(2 * besselI
#}

#' @title  Density for Generalized hyperbolic distribution
#' 
#' @description A function for evaluating the GH density.
#' 
#' @param x vector of points to be evaluated
#' @param delta location
#' @param mu    assymetric
#' @param sigma sscale
#' @param a     distribution parameter
#' @param b     distribution parameter
#' @param p     distribution parameter
#'
#'
#'
#'
dGH <- function(x, delta, mu, sigma, a, b, p, logd = TRUE){
  
  x_d <- x - delta
  x_sd <- x_d/sigma^2
  c1 <- a + x_sd*x_d
  c1 <- c1 * (b + mu^2/sigma^2)
  c1 <- sqrt(c1)
  
  loglik = mu*x_sd
  loglik <- loglik + (p - 0.5) * log(c1)
  loglik <- loglik + log(besselK(c1, p - 0.5, TRUE)) - c1
  
  c0 = p* log(b) + (0.5 -p) * log(b + mu^2/sigma^2)
  c0 = c0 - (0.5* p ) * log(a*b) - log(sigma) * log(besselK(sqrt(a*b), p, FALSE))
  loglik <- loglik + c0
  if(logd)
    return(loglik)
  
  return(exp(loglik))
}
#' @title  Density for t-distribution
#' @description density of t distribution
dtv2 <- function(x, delta, mu, sigma, nu, logd=TRUE){
  
  p = -nu
  a = 2 * (nu + 1)
  b = 0
  x_d <- x - delta
  x_sd <- x_d/sigma^2
  c1 <- a + x_sd*x_d
  c1 <- c1 * (b + mu^2/sigma^2)
  c1 <- sqrt(c1)
  c2 <- (a + x_sd*x_d)/(b + mu^2/sigma^2)
  
  loglik = mu*x_sd
  loglik <- loglik + 0.5*(p - 0.5) * log(c2)
  loglik <- loglik + log(besselK(c1, p - 0.5, TRUE)) - c1 #last term is for expontially scaled bessel
  c0 =  - lgamma(-p) - 0.5*log(pi) - (-p-0.5) * log(2) - p*log(a)
  c0 <- c0  - log(sigma) 
  loglik <- loglik + c0
  if(logd)
    return(loglik)
  
  return(exp(loglik))
}

dvgamma <- function(x, delta, mu, sigma, nu, logd=TRUE){
  
  p = nu
  a = 0
  b = 2*nu
  x_d <- x - delta
  x_sd <- x_d/sigma^2
  c1 <- a + x_sd*x_d
  c1 <- c1 * (b + mu^2/sigma^2)
  c1 <- sqrt(c1)
  c2 <- (a + x_sd*x_d)/(b + mu^2/sigma^2)
  
  loglik = mu*x_sd
  loglik <- loglik + 0.5*(p - 0.5) * log(c2)
  loglik <- loglik + log(besselK(c1, p - 0.5, TRUE)) - c1 #last term is for expontially scaled bessel
  c0 =  - lgamma(p) - 0.5*log(pi) + (1-p) * log(2) - 0.5*log(a/2)
  c0 <- c0  - log(sigma) 
  loglik <- loglik + c0
  if(logd)
    return(loglik)
  
  return(exp(loglik))
}

#' @title Plot for mixed effects model fit.
#'
#' @description A function to plot the results of mixed effects model fit.
#' @param res   A fitted object from mixed model fit.
#' @param dRE   A numeric vector for dimension to plot.
#' @details STUFF.
#' @return A plot.
#' @examples
#'   \dontrun{
#'   plot_GH_noise(...)
#'   }

plot_GH_noise  <- function(res, dRE = c())
{
  resid <- c()
  for(i in 1:length(res$Y))
  {
    Y <- res$Y[[i]] - res$mixedEffect_list$B_fixed[[i]]%*%as.vector(res$mixedEffect_list$beta_fixed)
    beta <- res$mixedEffect_list$beta_random + res$mixedEffect_list$U[,i]
    resid <- rbind(resid, Y-res$mixedEffect_list$B_random[[i]]%*%beta)
  }
  
  range_y <- c(min(resid), max(resid))
  x_ <- seq(range_y[1], range_y[2],length=100)
  if(res$measurementError_list$noise == 'NIG')
    f <- dnig(x_, 0, 0, res$measurementError_list$nu, res$measurementError_list$sigma)
  else
    f <- dnorm(x_,sd = res$measurementError_list$sigma)
  par(mfrow= c( ceiling((length(dRE) + 1)/2), 2) )
  hist(resid,50,prob=T, main='residual', xlab='x')
  lines(x_,
        f,
        col='red')
  if( length(dRE) > 0){
      for(i in 1:length(dRE))
      {
        j = dRE[i]
        U = res$mixedEffect_list$U[j,]
        x_ <- seq(min(U), max(U),length=100)
        hist(U,30,prob=T, main=paste('sample of centered RE ',j,sep=""), xlab='x')
        if(res$mixedEffect_list$noise == 'NIG'){
          f <- dnig(x_,
                   -res$mixedEffect_list$mu[j],
                   res$mixedEffect_list$mu[j],
                   res$mixedEffect_list$nu,
                   sqrt(res$mixedEffect_list$Sigma[j, j]))
        }else{
          f <- dnorm(x_,sd = sqrt(res$mixedEffect_list$Sigma[j, j]))
        }
        lines(x_,
              f,
              col='red')
      }
    }
}
