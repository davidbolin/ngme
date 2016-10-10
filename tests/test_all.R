graphics.off()
library(LDMod)

nIter <- 5
n.pers <- 100
nSim  <- 3
n.obs  <- 10 + 10*(1:n.pers)
n <- 100
n.pred <- 100
nBurnin = 100
pred.type <- "Filter"
pSubsample = 0.1
subsample.type = 2
test.pred = TRUE
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
  B_random[[i]] <- cbind(rep(1, n.obs[i]), (1:n.obs[i]) / n.obs[i] )
  B_random.pred[[i]] <- cbind(rep(1, n.pred), (1:n.pred) / n.pred )

  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  locs.pred[[i]] <- seq(from = 1, to = n.obs[i], length.out = n.pred)
  Vin[[i]] <- rep(1, n.obs[i])

  B_fixed[[i]]  <- as.matrix(locs[[i]])
  B_fixed.pred[[i]]  <- as.matrix(locs.pred[[i]])
}
mError_list <- list(Vs = Vin, noise = "NIG", sigma = 0.1, nu = 1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(2,-1)),
                          beta_fixed  = as.matrix(c(.1)),
                          Sigma = diag(c(0.1, 0.2)),
                          noise = "Normal",
                          Sigma_epsilon=1)

operator_list <- create_operator(locs, n, name = "fd2")

#operator_list$type  <- "matern"
operator_list$type  <- "fd2"
#operator_list$kappa <- -2
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
                    measurment_list  = mError_list,
                    processes_list   = processes_list,
                    operator_list    = operator_list,
                    pSubsample = pSubsample,
                    subsample.type = subsample.type)


n.plots <- 3
if(n.plots <= 3) {
  par(mfrow = c(1,3))
} else {
  par(mfrow = c(3,3))
}
matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="fixed effects")
lines(rep(mixedEffect_list$beta_fixed,length(res.est$mixedEffect_list$betaf_vec)),col=nIter)
matplot(res.est$mixedEffect_list$betar_vec,type="l",main="random effects",col=1)
matplot(t(matrix(rep(mixedEffect_list$beta_random,nIter),nrow=length(mixedEffect_list$beta_random),ncol = nIter)),add=TRUE,
        col=2,type="l")

plot(res.est$operator_list$tauVec,type="l",main="process tau")
lines(rep(operator_list$tau,nIter),col=2)

if(test.pred){
  res <- predictLong( Y = sim_res$Y,
                      locs.pred = locs.pred,
                      Brandom.pred = B_random.pred,
                      Bfixed.pred = B_fixed.pred,
                      type = pred.type,
                      nSim             = nSim,
                      locs             = locs,
                      mixedEffect_list = mixedEffect_list,
                      measurment_list  = mError_list,
                      processes_list   = processes_list,
                      operator_list    = operator_list,
                      return.samples = TRUE,
                      quantiles = c(0.05,0.95),
                      max.num.threads = 1)


  #x11()
  k = 1
  plot(locs[[k]],sim_res$Y[[k]],ylim=c(min(res$X.summary[[k]]$quantiles[[1]]$field),
                                       max(res$X.summary[[k]]$quantiles[[2]]$field)))
  lines(res$locs[[k]],res$X.summary[[k]]$Mean)
  lines(res$locs[[k]],res$X.summary[[k]]$quantiles[[1]]$field)
  lines(res$locs[[k]],res$X.summary[[k]]$quantiles[[2]]$field)

}
