graphics.off()
library(LDMod)

n.threads <- 1
nIter <- 500
n.pers <- 100
nSim  <- 2
n.obs  <- 10 + 2*(1:n.pers)
n <- 100
n.pred <- 100
nBurnin = 10
pred.type <- "Filter"
pSubsample = 0.1
operator.type = "matern"
#subsample.type = 2
test.pred = FALSE
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
  B_random[[i]] <- cbind(rep(1, n.obs[i]), 1000*(1:n.obs[i]))
  B_random.pred[[i]] <- cbind(rep(1, n.pred), (1:n.pred))

  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  locs.pred[[i]] <- seq(from = 1, to = n.obs[i], length.out = n.pred)
  Vin[[i]] <- rep(1, n.obs[i])

  B_fixed[[i]]  <- cbind(sqrt(locs[[i]]),1/locs[[i]])
  B_fixed.pred[[i]]  <- cbind(sqrt(locs.pred[[i]]),1/locs.pred[[i]])
}
mError_list <- list(Vs = Vin, noise = "NIG", sigma = 0.1, nu = 1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(1,2)),
                          beta_fixed  = as.matrix(c(1,2)),
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

if(operator.type == "matern"){
  operator_list$kappa <- 1
}
#operator_list$tau   <- 15


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
                    silent = FALSE)


n.plots <- 3
if(operator.type == "matern")
  n.plots = 4

if(n.plots <= 3) {
  par(mfrow = c(1,3))
} else {
  par(mfrow = c(2,2))
}
matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1)
#matplot(t(matrix(rep(mixedEffect_list$beta_fixed,nIter),nrow=length(mixedEffect_list$beta_fixed),ncol = nIter)),add=TRUE,col=1,type="l")
matplot(res.est$mixedEffect_list$betar_vec,type="l",main="random effects",col=1)
#matplot(t(matrix(rep(mixedEffect_list$beta_random,nIter),nrow=length(mixedEffect_list$beta_random),ncol = nIter)),add=TRUE,col=1,type="l")

plot(res.est$operator_list$tauVec,type="l",main="process tau")
lines(rep(operator_list$tau,nIter),col=2)
if(operator.type == "matern"){
  plot(res.est$operator_list$kappaVec,type="l",main="process kappa")
  lines(rep(operator_list$kappa,nIter),col=2)

}
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
                      max.num.threads = n.threads)


  #x11()
  k = 1
  plot(locs[[k]],sim_res$Y[[k]],ylim=c(min(res$X.summary[[k]]$quantiles[[1]]$field),
                                       max(res$X.summary[[k]]$quantiles[[2]]$field)))
  lines(res$locs[[k]],res$X.summary[[k]]$Mean)
  lines(res$locs[[k]],res$X.summary[[k]]$quantiles[[1]]$field)
  lines(res$locs[[k]],res$X.summary[[k]]$quantiles[[2]]$field)

}
