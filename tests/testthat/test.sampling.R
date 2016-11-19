#seem to work

rm(list = ls())
set.seed(1)
library(testthat)
library(LDMod)
n <- 1000
sim <- 1000

m <- 100
w <- rgamma(n, 10, 1)
w <- w/sum(w)
Y <- c()
for(i in 1:sim)
  Y <- c(Y, sampleR(m, w))
res<-hist(Y,breaks = 0:999, plot = F)
Fn <- ecdf(Y)
#nice brownian bridge
#plot(Fn(0:999)-cumsum(w))
EY = length(Y)*w
chi2 <- sum((EY-length(Y)*(Fn(0:(n-1))-c(0,Fn(0:(n-2)))))^2/EY)
qchisq(0.95,n)  < chi2

test_that("sample chisquare", {
  expect_lt(chi2, qchisq(0.95,n))
})