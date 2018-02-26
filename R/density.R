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

dnig <- function(x, delta, mu, nu, sigma)
{
  c0 <- sqrt(nu) * sqrt( mu^2/sigma^2 + nu) / pi
  f <- exp(nu + mu * (x - delta) /sigma^2)
  coeff <- sqrt( nu * sigma^2 + (x - delta)^2)
  f <- f * c0/ coeff
  f <- f * besselK(coeff  * sqrt(mu^2/sigma^4 + nu/sigma^2), 1)
  return(f)
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
  for(i in 1:length(res$obs_list))
  {
    Y <- res$obs_list[[i]]$Y - res$mixedEffect_list$B_fixed[[i]]%*%res$mixedEffect_list$beta_fixed
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
