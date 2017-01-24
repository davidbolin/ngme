
rm(list=ls())
library(LDMod)
library(MASS)
seed     <- 3
silent   <- 1
plotflag <- 1

nIter <- 1
n.pers <- 100
nSim  <- 2
n.obs  <- 10 + 0*(1:n.pers)

nBurnin = 40
pSubsample = 1
nPar_burnin = 100
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()
Vin <- list()

betaf <- as.matrix(c(2.1))
betar = as.matrix(c(0.8,0.5))

Sigma <- matrix(c(0.02,0,0,0.01),2,2)
sd_Y = 0.1
set.seed(seed)

Ve <- list()
for(i in 1:n.pers)
{
  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],rnorm(n.obs[i]))

  Y[[i]] <- rnorm(n = n.obs[i], B_fixed[[i]]%*%betaf, sd = sd_Y) + B_random[[i]]%*%mvrnorm(n = 1,mu = betar, Sigma = Sigma)
  Q =  solve(B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + sd_Y^2*diag(n.obs[i]))
  v = Y[[i]] - B_fixed[[i]]%*%betaf - B_random[[i]]%*%betar
  Ve[[i]] <- runif(n.obs[i]) + 1
  }



mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = betar,
                          beta_fixed  = betaf,
                          Sigma = Sigma,
                          noise = "NIG")


res <- estimateME(Y = Y,
                  mixedEffect_list = mixedEffect_list,
                  measurment_list = list(sigma = sd_Y, noise = "Normal"),
                  nSim = nSim,
                  alpha = 0.3,
                  pSubsample = pSubsample,
                  step0 = 0.3,
                  nIter = nIter,
                  silent = silent,
                  polyak_rate = -1,
                  seed = seed,
                  nBurnin = nBurnin,
                  estimate_fisher = 2)

print(t(res$FisherMatrix) - res$FisherMatrix)
if(1){
res <- estimateME(Y = Y,
                  mixedEffect_list = mixedEffect_list,
                  measurment_list = list(sigma = sd_Y, noise = "NIG", nu = 1.1, Vs = Ve),
                  nSim = nSim,
                  alpha = 0.3,
                  pSubsample = pSubsample,
                  step0 = 0.3,
                  nIter = nIter,
                  silent = silent,
                  polyak_rate = -1,
                  seed = seed,
                  nBurnin = nBurnin,
                  estimate_fisher = 2)
print(t(res$FisherMatrix) - res$FisherMatrix)
}
