##
# simple test that verifies that the model can correctly idenitfy the parameters, 
# Normal noise,
# NIG   mixed effect.
# from simulated data
#
##
rm(list=ls())
library(testthat)
graphics.off()
library(LDMod)
library(MASS)
seed     <- 5
silent   <- 0
plotflag <- 1

nIter <- 200
pSubsample <- 0.1
learning_rate <- 0.
n.pers <- 5000 #number of patients
n.obs  <- 50 #number of obs per patient

COV_beta <- matrix(c(0.2,0.1,0.1,0.2), ncol = 2, nrow = 2)
sd_Y    <- 0.1 # error of the noise

Br_list <- list()
betar <- c(0.9,0.4)
betaf <- c(1.)
mu   <- c(0.2, -0.2)
nu <- 1
betar_list <- list()
Bf_list    <- list()
V_list     <- list()
Y_list     <- list()
set.seed(seed)
for(i in 1:n.pers)
{
  Bf_list[[i]]    <- as.matrix(runif(n = n.obs))
  Br_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  V <- LDMod::rGIG(-0.5, nu, nu, sample.int(10^6,1))
  V_list[[i]] <- V
  betar_list[[i]] <- betar - mu * 1  + V * mu +
                     sqrt(V) * mvrnorm(n = 1, mu  =c(0,0), Sigma = COV_beta)
  Y_list[[i]]        <- rnorm(n = n.obs,
                              Br_list[[i]]%*%betar_list[[i]]
                            + Bf_list[[i]]%*%betaf, sd = sd_Y)

}


meas_list <- list(sigma_eps = sd_Y, noise = "Normal")
mixedEffect_list <- list(B_random = Br_list, 
                         B_fixed  = Bf_list,
                         Sigma = COV_beta, 
                         beta_random = c(0.,0.),
                         beta_fixed  = c(0.),  
                         mu          = 0*as.matrix(mu),
                         nu          = as.matrix(1.) + nu,
                         noise = "NIG")


beta_mat <- t(matrix(unlist(betar_list), nrow= 2, ncol = n.pers))
if(0){
x11()
par(mfrow=c(2,1))
hist(beta_mat[,1],100)
hist(beta_mat[,2],100)
}
res <- estimateME(Y = Y_list, 
                  mixedEffect_list = mixedEffect_list,
                  measurment_list = meas_list,
                  nSim = 2,
                  alpha = 0.3,
                  pSubsample = pSubsample,
                  pSubsample2 = pSubsample,
                  subsample.type = 3,
                  step0 = 0.3,
                  nIter = nIter,
                  silent = silent,
                  learning_rate = learning_rate,
                  polyak_rate = -1,
                  seed = seed)
if(plotflag){
  x11()
  par(mfrow=c(3,1))
  mu_vec = res$mixedEffect_list$mu_vec
  plot(mu_vec[,2], type='l', col='red', ylim=c(min(mu_vec), max(mu_vec)))
  lines(mu_vec[,1], col='red')
  n_ <- length(mu_vec[,1])
  lines(c(1, n_), c(mu[1],mu[1]))
  lines(c(1, n_), c(mu[2],mu[2]))
  betar_vec = res$mixedEffect_list$betar_vec
  plot(betar_vec[,1], type='l', col='red', ylim=c(min(betar_vec), max(betar_vec)))
  lines(betar_vec[,2], col='red')
  lines(c(1, n_), c(betar[1], betar[1]))
  lines(c(1, n_), c(betar[2], betar[2]))
  
  plot(res$mixedEffect_list$nu_vec, type='l', col='red' )
  lines(c(1, n_), c(nu, nu))
}
res <- estimateME(Y = Y_list, 
                  mixedEffect_list = mixedEffect_list,
                  measurment_list = meas_list,
                  nSim = 2,
                  alpha = 0.3,
                  pSubsample = 2*pSubsample,
                  pSubsample2 = pSubsample,
                  subsample.type = 1,
                  step0 = 0.3,
                  nIter = nIter,
                  silent = silent,
                  learning_rate = learning_rate,
                  polyak_rate = -1,
                  seed = seed)
if(plotflag){
  x11()
  par(mfrow=c(3,1))
  mu_vec = res$mixedEffect_list$mu_vec
  plot(mu_vec[,2], type='l', col='red', ylim=c(min(mu_vec), max(mu_vec)))
  lines(mu_vec[,1], col='red')
  n_ <- length(mu_vec[,1])
  lines(c(1, n_), c(mu[1],mu[1]))
  lines(c(1, n_), c(mu[2],mu[2]))
  betar_vec = res$mixedEffect_list$betar_vec
  plot(betar_vec[,1], type='l', col='red', ylim=c(min(betar_vec), max(betar_vec)))
  lines(betar_vec[,2], col='red')
  lines(c(1, n_), c(betar[1], betar[1]))
  lines(c(1, n_), c(betar[2], betar[2]))
  
  plot(res$mixedEffect_list$nu_vec, type='l', col='red' )
  lines(c(1, n_), c(nu, nu))
}

test_that("NIG-Gaussian random",
{
  expect_equal(c(res$mixedEffect_list$beta_random), betar, tolerance  = 0.1)
})
test_that("NIG-Gaussian fixed",
{
  expect_equal(c(res$mixedEffect_list$beta_fixed), betaf, tolerance  = 0.1)
})
test_that("NIG-Gaussian mu",
{
  expect_equal(c(res$mixedEffect_list$mu), mu, tolerance  = 0.2)
})
test_that("NIG-Gaussian sigma",
{
  expect_equal(c(res$measurementError_list$sigma), sd_Y, tolerance  = 0.01)
})
test_that("NIG-Gaussian nu",
{
            expect_equal(c(res$mixedEffect_list$nu), nu, tolerance  = 0.2)
})