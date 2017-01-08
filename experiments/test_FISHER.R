graphics.off()
if(0){
library(LDMod)
library(MASS)
seed     <- 5
silent   <- 1
plotflag <- 1

nIter <- 2
pSubsample <- 1
nSim <- 3
n.pers <- 10 #number of patients
n.obs  <- 10 #number of obs per patient

sd_Y    <- 0.5 # error of the noise

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


meas_list <- list(sigma = sd_Y, noise = "Normal")
mixedEffect_list <- list(B_fixed  = Bf_list,
                         beta_fixed  = betaf,
                         noise = "Normal")



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
print(res$FisherMatrix)
print(B_sum)
}
library(LDMod)
library(MASS)
seed     <- 5
silent   <- 1
plotflag <- 1

nIter <- 100
n.pers <- 2
nSim  <- 10
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
sd_Y = 2
F_fixed <- 0
F_random <- 0
set.seed(seed)
for(i in 1:n.pers)
{
  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))
  
  Y[[i]] <- rnorm(n = n.obs[i], B_fixed[[i]]%*%betaf, sd = sd_Y) + B_random[[i]]%*%mvrnorm(n = 1,mu = betar, Sigma = Sigma)
  Q =  B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + sd_Y^2*diag(n.obs[i])
  F_fixed <-  F_fixed  + t(B_fixed[[i]] )%*%solve(Q,B_fixed[[i]])
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
                  nBurnin = nBurnin,
                  estimate_fisher = TRUE)
print(res$FisherMatrix[1:3,1:3])
print(F_random)
print(F_fixed)