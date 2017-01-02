library(LDMod)
library(MASS)
rm(list=ls())
seed     <- 5
silent   <- 1
plotflag <- 1

nIter <- 4000
n.pers <- 5
nSim  <- 2
n.obs  <- 5 + 0*(1:n.pers)

nBurnin = 50
pSubsample = 1
nPar_burnin = 100
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()
Vin <- list()

betaf <- as.matrix(c(2.1))
betar = as.matrix(c(0.4,0.2))

Sigma <- matrix(c(2,0,0,1),2,2)
sd_Y = 1
F_fixed <- 0
F_random <- 0
set.seed(seed)
for(i in 1:n.pers)
{
  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))

  Y[[i]] <- rnorm(n = n.obs[i], B_fixed[[i]]%*%betaf, sd = sd_Y) + B_random[[i]]%*%mvrnorm(n = 1,mu = betar, Sigma = Sigma)
  Vin[[i]] <- rep(1, n.obs[i])

  Q =  B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + sd_Y^2*diag(n.obs[i])
  F_fixed <- F_fixed + t(B_fixed[[i]])%*%solve(Q,B_fixed[[i]])
  F_random <- F_random + t(B_random[[i]])%*%solve(Q,B_random[[i]])
}


mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = betar,
                          beta_fixed  = betaf,
                          Sigma = Sigma,
                          noise = "Normal")


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
                  estimate_fisher = TRUE)

cat("Fixed: True = ",F_fixed, ", Estimated  =", res$FisherMatrix[1,1],"\n")
cat("Random:")
print(F_random)
print(res$FisherMatrix[2:3,2:3])
