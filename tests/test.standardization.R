graphics.off()
library(LDMod)

n.threads <- 1
nIter <- 100
n.pers <- 100
nSim  <- 2
n.obs  <- 10 + 0*(1:n.pers)
n <- 100
n.pred <- 100
nBurnin = 10
pSubsample = 0.5
operator.type = "matern"
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()
Vin <- list()
for(i in 1:n.pers)
{
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(1.3))

  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  Vin[[i]] <- rep(1, n.obs[i])

  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
}
mError_list <- list(Vs = Vin, noise = "NIG", sigma = 0.1, nu = 1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(0.1,0.2)),
                          beta_fixed  = as.matrix(c(0.1)),
                          Sigma = diag(c(0.1, 0.2)),
                          noise = "Normal",
                          Sigma_epsilon=1)


operator_list <- create_operator(locs, n, name = operator.type)
if(operator.type == "matern"){
  operator_list$kappa <- 2
}

operator_list$tau   <- 15



processes_list = list(noise = "Normal",
                      nu  = 0.,
                      mu  = 0.)
processes_list$V <- list()
for(i in 1:length(locs))
{
  processes_list$V[[i]] <- operator_list$h
}
###
# simulation
#
###
mixedEffect_list_in = mixedEffect_list

sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list_in,
                              measurment_list   = mError_list,
                              processes_list    = processes_list,
                              operator_list     = operator_list)



processes_list$X <- sim_res$X


res.est <- estimateLong(Y                = sim_res$Y,
                        nIter            = nIter,
                        nSim             = nSim,
                        locs             = locs,
                        nBurnin           = nBurnin,
                        mixedEffect_list = mixedEffect_list,
                        nBurnin_learningrate = 0,
                        measurment_list  = mError_list,
                        processes_list   = processes_list,
                        operator_list    = operator_list,
                        pSubsample = pSubsample,
                        silent = FALSE,
                        standardize.mixedEffects = FALSE)


res.est2 <- estimateLong(Y                = sim_res$Y,
                        nIter            = nIter,
                        nSim             = nSim,
                        locs             = locs,
                        nBurnin           = nBurnin,
                        mixedEffect_list = mixedEffect_list,
                        nBurnin_learningrate = 0,
                        measurment_list  = mError_list,
                        processes_list   = processes_list,
                        operator_list    = operator_list,
                        pSubsample = pSubsample,
                        silent = FALSE,
                        standardize.mixedEffects = TRUE)


par(mfrow = c(2,2))
matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="FE",col=1)
matplot(res.est$mixedEffect_list$betar_vec,type="l",main="RE",col=1)
matplot(res.est2$mixedEffect_list$betaf_vec,type="l",main="FE, stand",col=2)
matplot(res.est2$mixedEffect_list$betar_vec,type="l",main="RE, stand",col=1)
