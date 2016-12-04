graphics.off()
library(LDMod)

n.threads <- 1
nIter <- 1000
n.pers <- 1000
nSim  <- 2
n.obs  <- 30 + 0*(1:n.pers)
n <- 100
nBurnin = 50
pSubsample = 0.1
nPar_burnin = 100
operator.type = "fd2"
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()
Vin <- list()

kappa_true = 2
tau_true = 50

for(i in 1:n.pers)
{
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))

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
                          noise = "NIG",
                          Sigma_epsilon=1,
                          nu = 1,
                          mu = matrix(c(1,1),2,1))


operator_list <- create_operator(locs, n, name = operator.type)
if(operator.type == "matern"){
  operator_list$kappa <- kappa_true
}

operator_list$tau   <- tau_true

processes_list = list(noise = "NIG",
                      nu  = 1,
                      mu  = 1)
processes_list$V <- list()
for(i in 1:length(locs))
{
  processes_list$V[[i]] <- operator_list$h
}

sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = mError_list,
                              processes_list    = processes_list,
                              operator_list     = operator_list)

res.est <- estimate.wrapper(Y = sim_res$Y,
                            locs = locs,
                            B_random= B_random,
                            B_fixed = B_fixed,
                            operator.type = "fd2",
                            n.process = NULL,
                            measurement.distribution = "NIG",
                            random.effect.distribution = "NIG",
                            process.distribution = "NIG",
                            estimation.options = list(nIter.gauss = 10,nIter = nIter,
                                                      pSubsample = 0.1,
                                                      nPar_burnin = nPar_burnin))




par(mfrow = c(3,3))

matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1)
matplot(res.est$mixedEffect_list$betar_vec,type="l",main="random effects",col=1)
matplot(res.est$mixedEffect_list$Sigma_vec,type="l",main="RE Sigma",col=1)
plot(res.est$operator_list$tauVec,type="l",main="process tau")
lines(rep(tau_true,nIter),col=2)
if(operator.type == "matern"){
  plot(res.est$operator_list$kappaVec,type="l",main="process kappa")
  lines(rep(kappa_true,nIter),col=2)
}
plot(res.est$measurementError_list$nu_vec,type="l",main="error nu")
plot(res.est$processes_list$nu_vec,type="l",main="process nu")
plot(res.est$mixedEffect_list$nu_vec,type="l",main="mixed nu")

matplot(res.est$mixedEffect_list$mu_vec,type="l",main="RE mu",col=1)
plot(res.est$processes_list$mu_vec,type="l",main="process mu")
truth =  c(mixedEffect_list$beta_fixed,
           mixedEffect_list$beta_random,
           mixedEffect_list$Sigma[c(1,2,4)],
           mError_list$sigma,
           operator_list$tau,
           mError_list$nu,
           mixedEffect_list$nu,
           mixedEffect_list$mu,
           processes_list$nu,
           processes_list$mu)

start.values  = c(res.est$mixedEffect_list$betaf_vec[1],
                  res.est$mixedEffect_list$betar_vec[1,],
                  res.est$mixedEffect_list$Sigma_vec[1,c(1,2,4)],
                  res.est$measurementError_list$sigma_vec[1],
                  res.est$operator_list$tauVec[1],
                  res.est$measurementError_list$nu_vec[1],
                  res.est$mixedEffect_list$nu_vec[1],
                  res.est$mixedEffect_list$mu_vec[1,],
                  res.est$processes_list$nu_vec[1],
                  res.est$processes_list$mu_vec[1])

estimates =  c(res.est$mixedEffect_list$betaf_vec[nIter],
               res.est$mixedEffect_list$betar_vec[nIter,],
               res.est$mixedEffect_list$Sigma_vec[nIter,c(1,2,4)],
               res.est$measurementError_list$sigma_vec[nIter],
               res.est$operator_list$tauVec[nIter],
               res.est$measurementError_list$nu_vec[nIter],
               res.est$mixedEffect_list$nu_vec[nIter],
               res.est$mixedEffect_list$mu_vec[nIter,],
               res.est$processes_list$nu_vec[nIter],
               res.est$processes_list$mu_vec[nIter])

result = data.frame(start = start.values,
                    estimate = estimates,
                    true = truth,
                    row.names = c("beta.fixed",
                                  "beta.random1","beta.random2",
                                  "Sigma1", "Sigma2","Sigma3",
                                  "sigma",
                                  "tau",
                                  "nu.error",
                                  "nu.mixed",
                                  "mu.mixed1",
                                  "mu.mixed2",
                                  "nu.process",
                                  "mu.process"))
print(result)
