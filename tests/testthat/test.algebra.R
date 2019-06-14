##
# testing third and fourth moment
#
##
rm(list=ls())
library(MASS)
library(testthat)
library(ngme)
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