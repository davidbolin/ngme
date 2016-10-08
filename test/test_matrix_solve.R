###
# simple function to ensure that prediagonalised solver works correct!
# 
###


rm(list=ls())
library(testthat)
library(LDMod)
#library(rGIG)
library(methods)
nobs <- 10
locs <- list()
locs[[1]]   <- seq(0, 1, length = nobs)
n <- 12
operator_list <- create_operator(locs, n, name = "fd2")
z <- rnorm(n = n)
b <- rnorm(n = n) + 1
input_list <- list(Q=  operator_list$Q[[1]]%*% t(operator_list$Q[[1]]), z = z, b = b)
fail <- test_PreDiagsolver(input_list)
test_that("testing diagonalization stabilization for sampling!",{
  expect_equal(0, fail, tolerance  = 1e-8)
})
