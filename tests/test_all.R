graphics.off()
library(LDMod)

nIter <- 50
n.pers <- 1
nSim  <- 10
n.obs  <- 10
n <- 100
n.pred <- 100
pred.type <- "Filter"


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
  B_random[[i]] <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  B_random.pred[[i]] <- cbind(rep(1, n.pred), (1:n.pred) / n.pred )

  Y[[i]] <- rep(1,n.obs)
  locs[[i]] <- 1:n.obs
  locs.pred[[i]] <- seq(from = 1, to = n.obs, length.out = n.pred)
  Vin[[i]] <- rep(1, n.obs)

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

res <- estimateLong(Y                = sim_res$Y,
                    nIter            = nIter,
                    nSim             = nSim,
                    locs             = locs,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list  = mError_list,
                    processes_list   = processes_list,
                    operator_list    = operator_list)


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
                      quantiles = c(0.05,0.95))


x11()
k = 1
plot(locs[[k]],sim_res$Y[[k]],ylim=c(min(res$X.summary[[k]]$quantiles[[1]]$field),
                                     max(res$X.summary[[k]]$quantiles[[2]]$field)))
lines(res$locs[[k]],res$X.summary[[k]]$Mean)
lines(res$locs[[k]],res$X.summary[[k]]$quantiles[[1]]$field)
lines(res$locs[[k]],res$X.summary[[k]]$quantiles[[2]]$field)
