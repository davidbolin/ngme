EiV_NIG <- function(U, Sigma, mu, nu)
{
  p <- -0.5*(1+length(U))
  b <- t(U + mu)%*%solve(Sigma, U + mu) + nu
  a <- mu%*%solve(Sigma, mu) + nu
  sqrt_ab = sqrt(a * b)
  K1 <- besselK(sqrt_ab, p, expon.scaled=T)
  K0 <- besselK(sqrt_ab, p+1, expon.scaled=T)
  
  sqrt_a_div_b <- sqrt(a/b)
  EiV = K0 / K1
  EiV = EiV * sqrt_a_div_b - (2 * p) * 1/b
  return(EiV)
}
dEiV_NIG <- function(U, Sigma, mu, nu)
{
  p <- -0.5*(1+length(U))
  b <- t(U + mu)%*%solve(Sigma, U + mu) + nu
  a <- mu%*%solve(Sigma, mu) + nu
  sqrt_ab = sqrt(a * b)
  K1 <- besselK(sqrt_ab, p, expon.scaled=F)
  K0 <- besselK(sqrt_ab, p+1, expon.scaled=F)
  
  sqrt_a_div_b <- sqrt(a/b)
  EiV = K0# / K1
  #EiV = EiV * sqrt_a_div_b - (2 * p) * 1/b
  
  K0dK1 = K0# / K1
  dEiV = 0
  dEiV = 1 - (p+1) * K0dK1 / sqrt_ab
  #dEiV = dEiV - (K0^2 + (p/sqrt_ab) *K1 * K0)/K1^2
  dEiV = dEiV * 0.5 * sqrt_a_div_b
  #dEiV = dEiV  * sqrt_a_div_b
  #dEiV = dEiV - 0.5 * K0dK1 * sqrt_a_div_b / b
  #dEiV = dEiV +  (2 * p) / b^2
  #dEiV = dEiV * 2 * solve(Sigma, U + mu)
 
  
  eps <- 10
  b_ <- t(U + mu)%*%solve(Sigma, U + mu) + nu + eps
  #sqrt_ab = sqrt(a * b_)
  #K1 <- besselK(sqrt(a * b_), p, expon.scaled=F)
  K0 <- besselK(sqrt(a * b_), p+1, expon.scaled=F)
  
  #sqrt_a_div_b <- sqrt(a/b_)
  EiV2 = K0# / K1
  EiV2 = EiV2 #* sqrt_a_div_b - (2 * p) * 1/b_
  # - 0.5*K0/K1 * sqrt{a/b} * 1/b  * 2*solve(Sigma, U + mu)
  # K'0 = (K1 - p K0 /sqrt{a*b}) *  2*solve(Sigma, U + mu)
  # 1/K1 = - (K0 + p K1 /sqrt{a*b})/K1^2
  print("***")
  print(dEiV)
  print((EiV2 - EiV )/eps)
  return(dEiV)
}

##
# computes dU and E[ddU] and E[ddU |U, Y] 
#
#
##
dU_NIG <- function(U, res, B, sigma, Sigma, mu, nu)
{
  # - (res-BU)^T(res-BU )/(2*\sigma^2) - \frac{1}{2V}(U - \mu V + \mu)^T \Sigma^{-1}(U - \mu V + \mu)
  #  dU = - (res - BU)^T B/sigma^2 -  \frac{1}{V}(U - \mu V + \mu)^T \Sigma^{-1}
  #    = (res - BU)^T B/sigma^2 -  (U V^{-1} - \mu  + \mu V^{-1})^T \Sigma^{-1}
  # E[V^{-1} | U ] = 
  # ddU = - B^TB /sigma^2 - \frac{\Sigma^{-1}}{V}
  # 
  dU = - sigma^(-2)*t(res - B%*%U)%*%B
  EiV = EiV_NIG(U, Sigma, mu, nu)
  dU = dU + solve(Sigma,EiV * U - mu + EiV * mu)
  ddU = - sigma^(-2)*t(B)%*%B - c(EiV) * solve(Sigma)
  EddU = - sigma^(-2)*t(B)%*%B - (1+2/nu) * solve(Sigma)
  return(list(dU = dU, ddU = -ddU, EddU = -EddU))
}

MALA <- function(U, res, B, sigma, Sigma, mu, nu)
{
  dU_list <- dU_NIG(U, res, B, sigma, Sigma, mu, nu)
  dU_list$ddU <- dU_list$EddU*4
  R <- chol(dU_list$ddU)
  U_mean <- U + solve(dU_list$ddU, t(dU_list$dU))/2
  Ustar = U_mean + solve(R, rnorm(length(U)))
  Ustar = c(Ustar)
  dU_star_list <- dU_NIG(Ustar, res, B, sigma, Sigma, mu, nu)
  dU_star_list$ddU <- dU_star_list$EddU*4
  U_star_mean <- Ustar + solve(dU_star_list$ddU, t(dU_star_list$dU))/2
  Rstar <- chol(dU_star_list$ddU)
  qstar <- sum( log( diag( R))) - 0.5* t(Ustar - U_mean)%*%dU_list$ddU%*%(Ustar - U_mean)
  q     <-sum( log( diag( Rstar)))  - 0.5* t(U - U_star_mean)%*%dU_star_list$ddU%*%(U - U_star_mean)
  loglik <- lNIG(U, res, B, sigma, solve(Sigma), mu, nu)
  loglikstar <- lNIG(Ustar, res, B, sigma, solve(Sigma), mu, nu)
  alpha <- loglikstar - qstar + q - loglik
  if(log(runif(1)) < alpha)
    U <- Ustar
  
  return(U)
}

lNIG <- function(U, res, B, sigma, iSigma, mu, nu)
{
  U_ <- U + mu
  b <- t(U_)%*%iSigma%*%U_ + nu
  a <- mu%*%iSigma%*%mu + nu
  
  logf = t(U_)%*%iSigma%*%mu
  logf = logf - 0.75 * log(b)
  # -1.5 * iSigma%*%U_/b # db /dU = 2* iSigma%*%U_
  sqrt_ab = sqrt(a * b)
  K1 <- besselK(sqrt_ab, -1.5, expon.scaled=T)
  # K_{-0.5} - 1.5 * K_{-1.5}/ sqrt{a*b}
  # / K_{-1.5} 
  # a * iSigma*(U_+mu) /\sqrt{a*b} =   iSigma*(U_+mu) * \sqrt{a}/\sqrt{b}
  
  # =iSigma*(U_+mu) * \sqrt{a}/\sqrt{b} * K_{-0.5}/K_{-1.5}
  # = -  1.5 * iSigma*(U_+mu) * /b #second bessel term
  # = -  1.5 * iSigma%*%U_/b
  # = iSigma*(U_+mu) *( \sqrt{a}/\sqrt{b} * 0.5 *K_{-0.5}/K_{-1.5}) - 3 /b)
  logf = logf + log(K1) - sqrt_ab
  logf = logf - sum((res - B%*%U)^2)/(2*sigma^2)
  return(logf)
}