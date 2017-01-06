require(testthat)
library(LDMod)
library(MASS)
graphics.off()
context("Normal Mixed effect")
seed     <- 3
silent   <- 0
plotflag <- 1

nIter <- 150
pSubsample <- 0.5
learning_rate <- 0.
n.pers <- 500 #number of patients
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
test_that("Normal measurement error", {
  for(i in 1:n.pers)
  {
    Bf_list[[i]]    <- as.matrix(runif(n = n.obs))
    Br_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
    V <- 1
    V_list[[i]] <- 1
    betar_list[[i]] <- betar + sqrt(V) * mvrnorm(n = 1, mu  =c(0,0), Sigma = COV_beta)
    Y_list[[i]]        <- rnorm(n = n.obs,
                                Br_list[[i]]%*%betar_list[[i]]
                                + Bf_list[[i]]%*%betaf, sd = sd_Y)
    
  }
  
  
  meas_list <- list(sigma_eps = 1, noise = "Normal")
  mixedEffect_list <- list(B_random = Br_list,
                           B_fixed  = Bf_list,
                           Sigma = COV_beta,
                           beta_random = c(0.,0.),
                           beta_fixed  = c(0.),
                           noise = "Normal")
  
  
  beta_mat <- t(matrix(unlist(betar_list), nrow= 2, ncol = n.pers))

  res <- estimateME(Y = Y_list,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list = meas_list,
                    nSim = 2,
                    alpha = 0.3,
                    pSubsample = pSubsample,
                    step0 = 0.3,
                    nIter = nIter,
                    silent = silent,
                    learning_rate = learning_rate,
                    polyak_rate = -1,
                    seed = seed)
  if(plotflag){
    x11()
    par(mfrow=c(2,2))
    betar_vec = res$mixedEffect_list$betar_vec
    plot(betar_vec[,1], type='l', col='red', ylim=c(min(betar_vec), max(betar_vec))
         ,main="beta random", xlab='i', ylab='beta')
    lines(betar_vec[,2], col='red')
    n_ <- length(betar_vec[,1])
    lines(c(1, n_), c(betar[1], betar[1]))
    lines(c(1, n_), c(betar[2], betar[2]))
    betaf_vec = res$mixedEffect_list$betaf_vec
    plot(betaf_vec[,1], type='l', col='red', ylim=c(min(betaf_vec), max(betaf_vec))
         ,main="beta fixed", xlab='i', ylab='beta')
    lines(c(1, n_), c(betaf[1], betaf[1]))
    
    sigma_vec = res$measurementError_list$sigma_vec
    plot(sigma_vec, type='l', col='red', ylim=c(min(sigma_vec), max(sigma_vec))
         ,main="sigma error", xlab='i', ylab='sigma')
    lines(c(1, n_), c(sd_Y[1], sd_Y[1]))
  }
  expect_equal(c(res$mixedEffect_list$beta_random), betar, tolerance  = 0.05)
  expect_equal(c(res$mixedEffect_list$beta_fixed), betaf, tolerance  = 0.05)
  expect_equal(c(res$measurementError_list$sigma), sd_Y, tolerance  = 0.01)
  
  
})



nIter   <- 700
#test_that("NIG measurement error", {
  for(i in 1:n.pers)
  {
    Bf_list[[i]]    <- as.matrix(runif(n = n.obs))
    Br_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
    V <- LDMod::rGIG(rep(-0.5,n.obs), rep(nu,n.obs), rep(nu,n.obs), sample.int(10^6,1))
    V_list[[i]] <- V
    betar_list[[i]] <- betar + mvrnorm(n = 1, mu  =c(0,0), Sigma = COV_beta)
    Y_list[[i]]        <- Br_list[[i]]%*%betar_list[[i]] 
    Y_list[[i]]        <-  Y_list[[i]]+ Bf_list[[i]]%*%betaf 
    Y_list[[i]]        <-  Y_list[[i]]   + sqrt(V) * rnorm(n = n.obs, 0, sd = sd_Y)
    
  }
  
  
  meas_list <- list(sigma_eps = 2*0.1,Vs = V_list, nu = 2*nu, noise = "NIG")
  mixedEffect_list <- list(B_random = Br_list,
                           B_fixed  = Bf_list,
                           Sigma = COV_beta,
                           beta_random = c(0.,0.),
                           beta_fixed  = c(0.),
                           noise = "Normal")
  
  
  beta_mat <- t(matrix(unlist(betar_list), nrow= 2, ncol = n.pers))
  
  res <- estimateME(Y = Y_list,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list = meas_list,
                    nSim = 2,
                    alpha = 0.3,
                    pSubsample = pSubsample,
                    step0 = 0.3,
                    nIter = nIter,
                    silent = silent,
                    learning_rate = learning_rate,
                    polyak_rate = -1,
                    seed = seed,
                    nBurnin = 10)
  if(plotflag){
    x11()
    par(mfrow=c(2,2))
    betar_vec = res$mixedEffect_list$betar_vec
    plot(betar_vec[,1], type='l', col='red', ylim=c(min(betar_vec), max(betar_vec))
         ,main="beta random", xlab='i', ylab='beta')
    lines(betar_vec[,2], col='red')
    n_ <- length(betar_vec[,1])
    lines(c(1, n_), c(betar[1], betar[1]))
    lines(c(1, n_), c(betar[2], betar[2]))
    betaf_vec = res$mixedEffect_list$betaf_vec
    plot(betaf_vec[,1], type='l', col='red', ylim=c(min(betaf_vec), max(betaf_vec))
         ,main="beta fixed", xlab='i', ylab='beta')
    lines(c(1, n_), c(betaf[1], betaf[1]))
    
    sigma_vec = res$measurementError_list$sigma_vec
    plot(sigma_vec, type='l', col='red', ylim=c(min(sigma_vec), max(sigma_vec))
         ,main="sigma error", xlab='i', ylab='sigma')
    lines(c(1, n_), c(sd_Y[1], sd_Y[1]))
    nu_vec = res$measurementError_list$nu_vec
    plot(nu_vec, type='l', col='red', ylim=c(min(nu_vec), max(nu_vec))
         ,main="nu error", xlab='i', ylab='nu')
    lines(c(1, n_), c(nu[1], nu[1]))
  }
  expect_equal(c(res$mixedEffect_list$beta_random), betar, tolerance  = 0.05)
  expect_equal(c(res$mixedEffect_list$beta_fixed), betaf, tolerance  = 0.05)
  expect_equal(c(res$measurementError_list$sigma), sd_Y, tolerance  = 0.01)
  expect_equal(c(res$measurementError_list$nu), nu, tolerance  = 0.3)
  
  
#})

