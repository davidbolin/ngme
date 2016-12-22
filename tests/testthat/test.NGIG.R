library(testthat)
library(LDMod)
context("N_GIG")
test_that("EiV_NIG_post", {
  
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