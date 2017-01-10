library(LDMod)
library(MASS)
rm(list=ls())
seed     <- 5
silent   <- 1
plotflag <- 1

nIter <- 40
n.pers <- 2
nSim  <- 10
n.obs  <- 10 + 0*(1:n.pers)

nBurnin = 40
pSubsample = 1
nPar_burnin = 100
n = 10
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()

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
  locs[[i]] <- 1:n.obs[i]
  Y[[i]] <- rep(1,n.obs[i])
}

operator_list <- create_operator(locs, n, name = "fd2")
operator_list$tau   <- 5
K = operator_list$tau*operator_list$Q[[1]]
Sigma.Z = solve(t(K)%*%diag(1/operator_list$h[[1]])%*%K)

processes_list = list(noise = "Normal")
processes_list$V <- list()
for(i in 1:length(locs))
{
  processes_list$V[[i]] <- operator_list$h
  A = spde.A(locs[[i]],operator_list$loc[[1]])
  Q =  B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + A%*%Sigma.Z%*%t(A) + sd_Y^2*diag(n.obs[i])
  F_fixed <-  F_fixed  + t(B_fixed[[i]] )%*%solve(Q,B_fixed[[i]])
  F_random <- F_random + t(B_random[[i]])%*%solve(Q,B_random[[i]])
  
}


mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = betar,
                          beta_fixed  = betaf,
                          Sigma = Sigma,
                          noise = "Normal")


sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = list(sigma = sd_Y, noise = "Normal"),
                              processes_list    = processes_list,
                              operator_list     = operator_list)

processes_list$X <- sim_res$X

res <- estimateLong(Y                = sim_res$Y,
                    locs             = locs,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list  = list(sigma = sd_Y, noise = "Normal"),
                    processes_list   = processes_list,
                    operator_list    = operator_list,
                    nSim             = nSim,
                    nBurnin_learningrate = 0,
                    alpha = 0.3,
                    pSubsample = pSubsample,
                    step0 = 0.3,
                    nIter = nIter,
                    silent = silent,
                    polyak_rate = -1,
                    seed = seed,
                    nBurnin = nBurnin,
                    estimate_fisher = TRUE)
