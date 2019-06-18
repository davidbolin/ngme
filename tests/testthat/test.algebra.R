##
# testing third and fourth moment
#
##
rm(list=ls())
library(MASS)
library(testthat)
library(ngme)
library(mvtnorm)
library(numDeriv)
library(mvtnorm)
graphics.off()
set.seed(2)
sim = 5000
d = 2
d2 = 4
d_s <- 2
mu = rnorm(n=d)
mu2 = rnorm(n=d2)
A  = matrix(rnorm(n=d*d), ncol = d)
A2  = matrix(rnorm(n=d2*d2), ncol = d2)


Sigma = A%*%t(A)
Sigma <- diag(sqrt(1/diag(Sigma))) %*% Sigma %*% diag(sqrt(1/diag(Sigma)))
xKx <- matrix(0, nrow=sim, ncol= d*d)
for(i in 1:sim){
  x <- mvrnorm(1, mu = mu, Sigma = Sigma)
  xKx[i, ] <- kronecker(x,x)
}
ExKx <- colMeans(xKx)
VxKx <- cov(xKx)
E2xKx <- VxKx + ExKx%*%t(ExKx)

E2xKx_true <- getNormalouterKron(Sigma,mu)
test_that("Forth moment", {
  expect_equal(c(E2xKx_true),
               c(E2xKx),
               tolerance = 0.1)
})

mu = mu2
A  = A2
Sigma = A%*%t(A)
Sigma <- diag(sqrt(1/diag(Sigma))) %*% Sigma %*% diag(sqrt(1/diag(Sigma)))
xKx <- matrix(0, nrow=sim, ncol= d_s*d_s*d2)
for(i in 1:sim){
  x <- mvrnorm(1, mu = mu, Sigma = Sigma)
  xKx[i, ] <- c(kronecker(kronecker(x[1:d_s],x[1:d_s]), t(x)))
}
ExKx <- colMeans(xKx)

ExKx_true <- c(getXXty(Sigma[1:d_s,1:d_s],Sigma[1:d_s,],mu[1:d_s],mu))

test_that("third moment", {
  expect_equal(c(E2xKx_true),
               c(E2xKx),
               tolerance = 0.1)
})
n <- 3
Sigma <- Sigma[1:4,1:4]
Z <- rmvnorm(n, sigma = Sigma)
ZZt <- matrix(0, nrow=ncol(Sigma), ncol=ncol(Sigma))
for(i in 1:n){
  ZZt <- ZZt  + Z[i,]%*%t(Z[i,])
}
loglik <- function(sigma, Z, ZZt){
  d <- ncol(Z)
  n <- nrow(Z)
  Sigma <- matrix(0, ncol = d, nrow = d)
  Sigma[] <- sigma
  iSigma <- solve(Sigma)
  res <- -n*d*0.5*log(2*pi) -n*0.5 * determinant(Sigma)$modulus -0.5* sum(diag(iSigma%*%ZZt))
  #res <- res - 0.5 * sum(apply(Z, 1, function(x){t(x)%*%iSigma%*%x}))
  return(c(res))
}
Hnum <- matrix(hessian(loglik, c(Sigma), Z=Z, ZZt = ZZt,
                       method.args=list(eps=1e-4,d=0.00001)) , nrow=nrow(Sigma)^2)
iSigma <- solve(Sigma)
H <-HeissanSigmaMVN(ZZt, iSigma, Sigma,n)

test_that("MVN Hessian", {
  expect_equal(c(H),
               c(Hnum),
               tolerance = 0.001)
})