##
#  testing convergence of nonGaussian mixed effect
#  comparing to the case when the variance components are known
#  D:2019-03-16
##
library(testthat)
library(ngme)
rm(list=ls())
graphics.off()
set.seed(2)

n_iter <- 2000

nindv <- 100
n     <- 110

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
mu_mixed     <- -1.5
nu_mixed     <- 0.75#0.5


B_fixed  <- list()
B_random <- list()
B_sigma  <- list()
Y        <- list()
V_mixed  <- rep(nindv)


for(indv in 1:nindv){
  B_fixed[[indv]]  <- cbind(rep(1, n), rnorm(n))
  B_random[[indv]] <- as.matrix(1:n)
  V_mixed[indv] <-rGIG(rep(-0.5,1),
                 rep( nu_mixed, 1),
                 rep( nu_mixed, 1),
                 as.integer(1000 * runif(1) ))
  V_   <- V_mixed[indv]
  Y[[indv]]        <- B_fixed[[indv]]%*%beta_fixed +
                   B_random[[indv]]%*%(beta_random + mu_mixed*(-1+V_) + sqrt(V_)*sigma_random*rnorm(1)) +
                   sigma * rnorm(n)
  
  B_sigma[[indv]]  <- as.matrix(rep(1, n))
  
}

context("mixed with fixed V")
mixed_list <- list(B_fixed    = B_fixed, 
                   B_random   = B_random,
                   V          = V_mixed,
                   fixedV     = 1,
                   name       = 'NIG')
error_list <- list(name       = 'Normal',
                   B          = B_sigma)
res_Vknown <- test_mixed(n_iter, 
                   Y, 
                   mixed_list,
                   error_list)
mixed_list <- list(B_fixed    = B_fixed, 
                   B_random   = B_random,
                   V          = V_mixed,
                   fixedV     = 0,
                   name       = 'NIG')
res_Vunkown <- test_mixed(n_iter, 
                   Y, 
                   mixed_list,
                   error_list)



context("mixed")
test_that(" V known NIG,", {
  expect_equal(as.vector(res_Vknown$mixedEffect_list$beta_fixed),
               beta_fixed, 
               tolerance  = 10^-1)
  
  expect_equal(as.vector(res_Vknown$mixedEffect_list$beta_random),
               beta_random, 
               tolerance  = 10^-1)
  
})

test_that(" V unkknown, NIG", {
  expect_equal(as.vector(res_Vunkown$mixedEffect_list$beta_fixed),
               beta_fixed, 
               tolerance  = 10^-1)
  
  expect_equal(as.vector(res_Vunkown$mixedEffect_list$beta_random),
               beta_random, 
               tolerance  = 5*10^-1)
  
})
###
#plots for debugging
###
if(0){
x11()
par(mfrow=c(2,1))
plot(sqrt(res_Vknown$mixedEffect_list$Sigma_vec))
plot(sqrt(res_Vunkown$mixedEffect_list$Sigma_vec))
cat('sigma (V known)= ',sqrt(res_Vknown$mixedEffect_list$Sigma),'\n')
cat('sigma = ',sqrt(res_Vunkown$mixedEffect_list$Sigma),'\n')
x11()
x <-seq(-10,5,length.out = 2000)
dens <-dnig(x, as.vector(res_Vknown$mixedEffect_list$beta_random)-as.vector(res_Vknown$mixedEffect_list$mu), as.vector(res_Vknown$mixedEffect_list$mu), res_Vknown$mixedEffect_list$nu,  sqrt(res_Vknown$mixedEffect_list$Sigma[1,1]))
dens2 <-dnig(x, beta_random -mu_mixed,mu_mixed, nu_mixed, sigma_random)
dens3 <-dnig(x, as.vector(res_Vunkown$mixedEffect_list$beta_random)-as.vector(res_Vunkown$mixedEffect_list$mu), as.vector(res_Vunkown$mixedEffect_list$mu), res_Vunkown$mixedEffect_list$nu,  sqrt(res_Vunkown$mixedEffect_list$Sigma[1,1]))

plot(x,dens2,type='l')
lines(x, dens, col='red')
lines(x, dens3, col='blue')
x11()
par(mfrow=c(3,1))
hist(res_Vknown$mixedEffect_list$U,20)
hist(res_Vunkown$mixedEffect_list$U,20)
hist(res_Vunkown$mixedEffect_list$U-res_Vknown$mixedEffect_list$U)
}


