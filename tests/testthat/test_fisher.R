
library(testthat)

  seed     <- 5
  silent   <- FALSE

  library(LDMod)

  nIter <- 5000
  n.pers <- 5
  nSim  <- 2
  n.obs  <- 10 + 0*(1:n.pers)
  n <- 100
  error.dist = "Normal"
  mixed.dist = "Normal"
  nBurnin = 50
  pSubsample = 1
  nPar_burnin = 100
  Y <- list()
  locs <- list()
  B_random <- list()
  B_fixed  <- list()
  Vin <- list()


  for(i in 1:n.pers)
  {
    B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))

    Y[[i]] <- rep(1,n.obs[i])
    locs[[i]] <- 1:n.obs[i]
    Vin[[i]] <- rep(1, n.obs[i])

    B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  }
  mError_list <- list(Vs = Vin, noise = error.dist, sigma = 0.1, nu = 1)
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = as.matrix(c(0.1,0.2)),
                            beta_fixed  = as.matrix(c(0.1)),
                            Sigma = diag(c(0.1, 0.2)),
                            noise = mixed.dist,
                            Sigma_epsilon=1,
                            nu = 1,
                            mu = matrix(c(2,2),2,1))


  res <- estimateME(Y = Y,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list = mError_list,
                    nSim = nSim,
                    alpha = 0.3,
                    pSubsample = pSubsample,
                    step0 = 0.3,
                    nIter = nIter,
                    silent = silent,
                    polyak_rate = -1,
                    seed = seed,
                    estimate_fisher = TRUE)

  print(max(abs(B_sum/res$FisherMatrix[1:2,1:2]-1)))




context("Fisher")

test_that("Fisher, Gaussian fixed", {

graphics.off()
library(LDMod)
library(MASS)
seed     <- 5
silent   <- 1
plotflag <- 1

nIter <- 5000
pSubsample <- 1
nSim <- 2
n.pers <- 5 #number of patients
n.obs  <- 10 #number of obs per patient

sd_Y    <- 1 # error of the noise

betaf <- c(1.,-1)
nu <- 1
betar_list <- list()
Bf_list    <- list()
V_list     <- list()
Y_list     <- list()
set.seed(seed)
B_sum <- 0
for(i in 1:n.pers)
{
  Bf_list[[i]]    <- cbind(1:n.obs,rep(1, n.obs))
  Y_list[[i]]        <- rnorm(n = n.obs, Bf_list[[i]]%*%betaf, sd = sd_Y)
  B_sum <- B_sum + t(Bf_list[[i]])%*%Bf_list[[i]]/sd_Y^2
}


meas_list <- list(sigma_eps = sd_Y, noise = "Normal")
mixedEffect_list <- list(B_fixed  = Bf_list,
                         beta_fixed  = betaf,
                         noise = "normal")



res <- estimateME(Y = Y_list,
                  mixedEffect_list = mixedEffect_list,
                  measurment_list = meas_list,
                  nSim = nSim,
                  alpha = 0.3,
                  pSubsample = pSubsample,
                  step0 = 0.3,
                  nIter = nIter,
                  silent = silent,
                  polyak_rate = -1,
                  seed = seed,
                  estimate_fisher = TRUE)

expect_equal(max(abs(B_sum/res$FisherMatrix[1:2,1:2]-1)),0,tolerance=0.01)
})


test_that("Fisher, Gaussian mixed", {
  seed     <- 5
  silent   <- FALSE

  library(LDMod)

  nIter <- 5000
  n.pers <- 5
  nSim  <- 2
  n.obs  <- 10 + 0*(1:n.pers)
  n <- 100
  error.dist = "Normal"
  mixed.dist = "Normal"
  nBurnin = 50
  pSubsample = 1
  nPar_burnin = 100
  Y <- list()
  locs <- list()
  B_random <- list()
  B_fixed  <- list()
  Vin <- list()


  for(i in 1:n.pers)
  {
    B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))

    Y[[i]] <- rep(1,n.obs[i])
    locs[[i]] <- 1:n.obs[i]
    Vin[[i]] <- rep(1, n.obs[i])

    B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  }
  mError_list <- list(Vs = Vin, noise = error.dist, sigma = 0.1, nu = 1)
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = as.matrix(c(0.1,0.2)),
                            beta_fixed  = as.matrix(c(0.1)),
                            Sigma = diag(c(0.1, 0.2)),
                            noise = mixed.dist,
                            Sigma_epsilon=1,
                            nu = 1,
                            mu = matrix(c(2,2),2,1))


  res <- estimateME(Y = Y,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list = mError_list,
                    nSim = nSim,
                    alpha = 0.3,
                    pSubsample = pSubsample,
                    step0 = 0.3,
                    nIter = nIter,
                    silent = silent,
                    polyak_rate = -1,
                    seed = seed,
                    estimate_fisher = TRUE)

  print(max(abs(B_sum/res$FisherMatrix[1:2,1:2]-1)))

})
