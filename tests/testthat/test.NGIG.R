library(testthat)
library(LDMod)
context("N_GIG")
lNIG <- function(U, res, B, sigma, iSigma, mu, nu)
{
  p <- -0.5*(1+length(U))
  U_ <- U + mu
  b <- t(U_)%*%iSigma%*%U_ + nu
  a <- mu%*%iSigma%*%mu + nu
  logf = t(U_)%*%iSigma%*%mu
  logf = logf - 0.75 * log(b)
  sqrt_ab = sqrt(a * b)
  K1 <- besselK(sqrt_ab, p, expon.scaled=T)

  logf = logf + log(K1) - sqrt_ab
  logf = logf - sum((res - B%*%U)^2)/(2*sigma^2)
  return(logf)
}


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
dU_NIG <- function(U, res, B, sigma, Sigma, mu, nu)
{
  # 
  dU = - sigma^(-2)*t(res - B%*%U)%*%B
  EiV = EiV_NIG(U, Sigma, mu, nu)
  dU = dU + solve(Sigma,EiV * U - mu + EiV * mu)
  ddU = + sigma^(-2)*t(B)%*%B + c(EiV) * solve(Sigma)
  EddU = - sigma^(-2)*t(B)%*%B - (1+2/nu) * solve(Sigma)
  dEiV <- dEiV_NIG(U, Sigma, mu, nu)
  d_ <- dEiV%*%t(solve(Sigma,(U + mu)))
  d_ <- (d_ + t(d_))/2
  EddU = EddU - d_
  ddU  = ddU + d_
  return(list(dU = dU, ddU = ddU, EddU = -EddU))
}
dEiV_NIG <- function(U, Sigma, mu, nu)
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
  
  K0dK1 = K0 / K1
  dEiV = 0
  dEiV = -1 - (p+1) * K0dK1 / sqrt_ab
  dEiV = dEiV - (-K0^2 + (p/sqrt_ab) *K1 * K0)/K1^2
  dEiV = dEiV * 0.5 * sqrt_a_div_b
  dEiV = dEiV  * sqrt_a_div_b
  dEiV = dEiV - 0.5 * K0dK1 * sqrt_a_div_b / b
  dEiV = dEiV +  (2 * p) / b^2
  dEiV = c(dEiV * 2) * solve(Sigma, U + mu)
  return(dEiV)
}

test_that("EiV_NIG_post", {
  nu <- runif(1)+0.1
  p <- -0.5
  a <- nu
  b <- nu
  Sigma <- matrix(c(2,1,1,2), nrow=2)*(1 + runif(1))
  U     <- c(1.,2.)*runif(1)
  mu    <- c(2,5)*runif(1)
  delta <- -mu
  expect_equal(test_EiV_NGIG(U, Sigma, delta, mu, p, a, b),
               EiV_NIG(U, Sigma, mu, nu)[1],
               tolerance  = 10^-6)
})
test_that("dU_NIG_post", {
  res <- c(1, 2, 3.)
  sigma <- .2
  nu <- runif(1)+0.1
  p <- -0.5
  a <- nu
  b <- nu
  Sigma <- matrix(c(2,1,1,2), nrow=2)*(1 + runif(1))
  U     <- c(1.,2.)*runif(1)
  mu    <- c(2,5)*runif(1)
  B <- matrix(runif(3*2), ncol = 2)
  delta <- -mu
  dU_NIG <- dU_NIG(U, res, B, sigma, Sigma, mu, nu)
  

  res <- test_dU_EiV( U,
                      Sigma,
                      delta,
                      mu,
                      p,
                      a,
                      b,
                      res,
                      diag(length(res))/sigma^2,
                      B)
  expect_equal(res$dU, c(dU_NIG$dU), tolerance = 10^-6)
  expect_equal(res$ddU, dU_NIG$ddU, tolerance = 10^-6)
})


test_that("log(f_NIG)", {
  res <- c(1, 2, 3.)
  sigma <- .2
  B <- matrix(runif(3*2), ncol = 2)
  nu <- runif(1)+0.1
  p <- -0.5
  a <- nu
  b <- nu
  Sigma <- matrix(c(2,1,1,2), nrow=2)*(1 + runif(1))
  U     <- c(1.,2.)*runif(1)
  mu    <- c(2,5)*runif(1)
  expect_equal(c(lNIG(U, res, B, sigma, solve(Sigma), mu, nu)),
               test_logf_NIG(U, mu, -mu, solve(Sigma), nu) - sum((res - B%*%U)^2)/(2*sigma^2),
               tolerance = 10^-6)
})