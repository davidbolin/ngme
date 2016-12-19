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
silent   <- 1
plotflag <- 1

nIter <- 10
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

print(B_sum)
print(res$FisherMatrix[1:2,1:2])
print(B_sum/res$FisherMatrix[1:2,1:2])