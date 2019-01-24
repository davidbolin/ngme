##
#  test mixed effect p measurement error
#  D:2019-01-15
##
rm(list=ls())
graphics.off()
library(testthat)
library(ngme)
set.seed(4)

n_iter <- 1000

nindv <- 10
n     <- 110
nu      <- 0.6
mu      <- 1.
sigma        <- 0.5


Y        <- list()
V        <- list()
Vs       <- list()
for(indv in 1:nindv){
  V[[indv]]  <- rGIG(rep(-0.5,n),
                     rep( nu, n),
                     rep( nu, n),
                     as.integer(1000 * runif(1) ))
  E      = sigma * sqrt(V[[indv]]) * rnorm(n) + (V[[indv]]-1)*mu
  Y[[indv]]  <-   E
  Vs[[indv]] <- rep(1,n)
  
}
context("error")
test_that("NIG", {
  error_list <- list(name = 'NIG',
                     mu   = 0.,
                     nu   = 1.,
                     assymetric = 1,
                     Vs   = Vs
                    )
  res <- test_error(n_iter, 
                    Y, 
                    error_list)
  expect_equal(c(res$measurementError_list$nu,
                 res$measurementError_list$mu),
               c(nu,mu),
               tolerance  = 10^-1)
})