##
#  checking estimation again
#  D:2019-03-16
##
library(testthat)
library(ngme)
rm(list=ls())
graphics.off()
set.seed(3)

n_iter <- 2*1500

nindv <- 1200
n     <- 30

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
mu_mixed     <- 2
nu_mixed     <- 1#0.5


B_fixed  <- list()
B_random <- list()
B_sigma  <- list()
Y        <- list()
V_mixed  <- rep(nindv)

data <- c()
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
  id <- rep(indv, n)
  data <- rbind(data, cbind(B_fixed[[indv]],
                            B_random[[indv]],
                            Y[[indv]],
                            id))
}
dimnames(data)[[2]] <- c('B1','B2','B3','Y','id')

NIGMVD_ass <- ngme( fixed       = Y ~ B1 + B2,
                    random      = ~ -1+B3|id,
                    data        = as.data.frame(data),
                    reffects    = 'NIG',
                    use.process = F,
                    silent      = T,
                    controls.init = list(nIter.init=100),
                    nIter  = n_iter,
                    controls    = list(estimate.fisher = FALSE,
                                       pSubsample = 0.5,
                                       subsample.type  = 0,
                                       nSim  =2,
                                       nBurnin = 2,
                                       alpha = 0.01,
                                       step0 = 0.99))


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
x11()
par(mfrow=c(3,2))
plot(res_Vknown$mixedEffect_list$betar_vec,type='l')
plot(res_Vknown$mixedEffect_list$mu_vec,type='l')
plot(res_Vunkown$mixedEffect_list$betar_vec,type='l',col='red')
plot(res_Vunkown$mixedEffect_list$mu_vec,type='l',col='red')
plot(NIGMVD_ass$mixedEffect_list$betar_vec,type='l',col='blue')
plot(NIGMVD_ass$mixedEffect_list$mu_vec,type='l',col='blue')
