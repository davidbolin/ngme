library(testthat)
library(LDMod)
context("GIG")
test_that("EiV_GIG", {
  p <- runif(1) + 1
  a <- runif(1) + 1
  b <- runif(1) + 1
  EiV <- test_EiV_GIG(p, a, b)
  dGIG <- function(x, p, a ,b)
  {
    c <- (a/b)^(p/2) / (2 * besselK(sqrt(a*b),p) )
    f <- c*x^(p-1) * exp(-(a*x + b/x)/2)
    return(f)
  }
  
  
  EiVR <- integrate(function(x){ x^-1 * dGIG(x, p, a, b)}, 0, Inf)
  expect_equal(EiV, EiVR$value, tolerance  = 10^-6)
})
test_that("db_EiV_GIG", {
  p <- runif(1) + 1
  a <- runif(1) + 1
  b <- runif(1) + 1
  EiV   <- test_EiV_GIG(p, a, b)
  EiV_e <- test_EiV_GIG(p, a, b + 10^-6)
  db_EiV_num <- (EiV_e - EiV)/10^-6 
  db_EiV <- test_db_EiV_GIG(p, a, b)
  expect_equal(db_EiV, db_EiV_num, tolerance  = 10^-6)
})