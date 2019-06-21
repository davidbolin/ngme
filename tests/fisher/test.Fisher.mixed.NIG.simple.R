###
# For debugging extended fisher infromation
#
# D: 2019-03-13
###
graphics.off()
library(ngme)
library(testthat)
library(doParallel)
library(Matrix)

#data options
n.pers <- 30
n.obs  <- rep(10,n.pers)#10 + (1:n.pers)


#Fisher options
nIter.fisher = 1
nSim.fisher = 1000
nBurnin = 10

#simulate data
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()
locs.pred <- list()
B_random.pred <- list()
B_fixed.pred <- list()
Vin <- list()
for(i in 1:n.pers)
{
  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i] #sort(1 + 9*runif(n.obs[i]))
  #random effects, 1 and t
  B_random[[i]] <- cbind(rep(1, n.obs[i]), locs[[i]]/n.obs[i])
  #fixed effects, sqrt(t) and 1/t
  B_fixed[[i]]  <- cbind(sqrt(locs[[i]]),1/locs[[i]])
}
mError_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(0.3,0.25)),
                          beta_fixed  = as.matrix(c(0.1,0.2)),
                          Sigma = diag(c(0.1, 0.2)),
                          noise = "NIG",
                          Sigma_epsilon=0,
                          nu = 1,
                          mu = as.matrix(c(1,2)))

mixedEffect_list_in = mixedEffect_list

sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list_in,
                              measurment_list   = mError_list)

res.fisher <- estimateLong(Y                = sim_res$Y,
                           locs             = locs,
                           mixedEffect_list = mixedEffect_list,
                           measurment_list  = mError_list,
                           nIter = nIter.fisher,
                           nSim             = nSim.fisher,
                           silent = T,
                           nBurnin_base = nBurnin,
                           estimate_fisher = 2)  

RES <- solve(res.fisher$FisherMatrix[c(1:6),c(1:6)])
RES2 <- solve(res.fisher$FisherMatrix[c(1:6,10),c(1:6,10)])
RES3 <- solve(res.fisher$FisherMatrix)[1:6,1:6]