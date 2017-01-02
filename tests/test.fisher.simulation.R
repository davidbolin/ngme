graphics.off()
library(LDMod)

n.iterations <- 100
nIter <- 1000
n.pers <- 10
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


beta_fixed = as.matrix(c(0.1))
beta_random = as.matrix(c(0.1,0.2))
Sigma = diag(c(0.1, 0.2))
sigma.e = 0.3
kappa = 2
tau = 50
nu.error = 0.4
nu.mixed = 0.6
mu.mixed = matrix(c(1,1),2,1)
nu.process = 0.7
mu.process = 3

truth =  c(beta_fixed,
           beta_random,
           Sigma[c(1,2,4)],
           sigma.e,
           tau,
           nu.error,
           nu.mixed,
           mu.mixed,
           nu.process,
           mu.process)

estimates = matrix(data=NA,nrow = length(truth),ncol = n.iterations)

for(iteration in 1:n.iterations){

  for(i in 1:n.pers)
  {
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))

  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  Vin[[i]] <- rep(1, n.obs[i])

  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
}
  mError_list <- list(Vs = Vin, noise = "NIG", sigma = sigma.e, nu = nu.error)
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = beta_random,
                            beta_fixed  = beta_fixed,
                            Sigma = Sigma,
                            noise = "NIG",
                            Sigma_epsilon=1,
                            nu = nu.mixed,
                            mu = mu.mixed)


  operator_list <- create_operator(locs, n, name = operator.type)
  if(operator.type == "matern"){
    operator_list$kappa <- kappa
  }

  operator_list$tau   <- tau

  processes_list = list(noise = "NIG",
                        nu  = nu.process,
                        mu  = mu.process)
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
                              silent = TRUE,
                              estimation.options = list(nIter.gauss = 10,nIter = nIter,
                                                        pSubsample = 0.1,
                                                        nPar_burnin = nPar_burnin))


  estimates[,iteration] =  c(res.est$mixedEffect_list$betaf_vec[nIter],
                            res.est$mixedEffect_list$betar_vec[nIter,],
                            res.est$mixedEffect_list$Sigma_vec[nIter,c(1,2,4)],
                            res.est$measurementError_list$sigma_vec[nIter],
                            res.est$operator_list$tauVec[nIter],
                            res.est$measurementError_list$nu_vec[nIter],
                            res.est$mixedEffect_list$nu_vec[nIter],
                            res.est$mixedEffect_list$mu_vec[nIter,],
                            res.est$processes_list$nu_vec[nIter],
                            res.est$processes_list$mu_vec[nIter])

}

res.est <- estimate.wrapper(Y = sim_res$Y,
                            locs = locs,
                            B_random= B_random,
                            B_fixed = B_fixed,
                            operator.type = "fd2",
                            n.process = NULL,
                            measurement.distribution = "NIG",
                            random.effect.distribution = "NIG",
                            process.distribution = "NIG",
                            silent = TRUE,
                            estimate_fisher = TRUE,
                            estimation.options = list(nIter.gauss = 10,nIter = nIter,
                                                      pSubsample = 0.1,
                                                      nPar_burnin = nPar_burnin))


estimate.mean <- apply(estimates,1,mean)
estimate.var <- apply(estimates,1,var)
fisher <- diag(solve(res.est$FisherMatrix))
result = data.frame(mean = estimate.mean,
                    var = estimate.var,
                    fisher = fisher,
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
