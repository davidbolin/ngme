graphics.off()
library(LDMod)

n.iterations <- 10
nIter <- 40
n.pers <- 100
nSim  <- 2
n.obs  <- 30 + 0*(1:n.pers)
n <- 30
nBurnin = 50
pSubsample = 0.1
nPar_burnin = 100
operator.type = "fd2"
meas.dist = "Normal"
mixe.dist = "Normal"
proc.dist = "Normal"

use.proc = TRUE

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

truth =  c(beta_fixed, beta_random)

estimates = matrix(data=NA,nrow = length(truth),ncol = n.iterations)

for(iteration in 1:n.iterations){
  cat("iteration = ", iteration, "\n")
  for(i in 1:n.pers)
  {
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))

  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  Vin[[i]] <- rep(1, n.obs[i])

  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
}
  mError_list <- list(Vs = Vin, noise = meas.dist, sigma = sigma.e, nu = nu.error)
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = beta_random,
                            beta_fixed  = beta_fixed,
                            Sigma = Sigma,
                            noise = mixe.dist,
                            Sigma_epsilon=1,
                            nu = nu.mixed,
                            mu = mu.mixed)

  if(use.proc){
    operator_list <- create_operator(locs, n, name = operator.type)
    if(operator.type == "matern"){
      operator_list$kappa <- kappa
    }

    operator_list$tau   <- tau

    processes_list = list(noise = proc.dist,
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
                                n.process = n,
                                measurement.distribution = meas.dist,
                                random.effect.distribution = mixe.dist,
                                process.distribution = proc.dist,
                                silent = TRUE,
                                estimation.options = list(nIter.gauss = 100,nIter = nIter,
                                                          pSubsample = pSubsample,
                                                          nPar_burnin = nPar_burnin))

  } else {
    sim_res <- simulateLongPrior( Y                 = Y,
                                  locs              = locs,
                                  mixedEffect_list  = mixedEffect_list,
                                  measurment_list   = mError_list)


    res.est <- estimate.wrapper(Y = sim_res$Y,
                                locs = locs,
                                B_random= B_random,
                                B_fixed = B_fixed,
                                use.process = FALSE,
                                measurement.distribution = meas.dist,
                                random.effect.distribution = mixe.dist,
                                silent = TRUE,
                                estimation.options = list(nIter.gauss = 10,nIter = nIter,
                                                          pSubsample = 0.5,
                                                          nPar_burnin = nPar_burnin))


  }


  estimates[,iteration] =  c(res.est$mixedEffect_list$betaf_vec[nIter],
                            res.est$mixedEffect_list$betar_vec[nIter,])

}

if(use.proc){
  res <- estimate.wrapper(Y = sim_res$Y,
                              locs = locs,
                              B_random= B_random,
                              B_fixed = B_fixed,
                              operator.type = "fd2",
                              n.process = n,
                              measurement.distribution = meas.dist,
                              random.effect.distribution = mixe.dist,
                              process.distribution = proc.dist,
                              silent = TRUE,
                              estimate_fisher = TRUE,
                              estimation.options = list(nIter.gauss = 10,
                                                        nIter = 20,
                                                        nSim = 2,
                                                        nSim.fisher = 20,
                                                        pSubsample = 0.1,
                                                        nPar_burnin = nPar_burnin))

} else {
  res <- estimate.wrapper(Y = sim_res$Y,
                              locs = locs,
                              B_random= B_random,
                              B_fixed = B_fixed,
                              use.process = FALSE,
                              measurement.distribution = meas.dist,
                              random.effect.distribution = mixe.dist,
                              silent = TRUE,
                              estimate_fisher = TRUE,
                              estimation.options = list(nIter.gauss = 10,
                                                        nIter = 1000,
                                                        nSim = 2,
                                                        nSim.fisher = 1000,
                                                        pSubsample = 0.1,
                                                        nPar_burnin = nPar_burnin))

}


estimate.mean <- apply(estimates,1,mean)
estimate.var <- apply(estimates,1,var)
fisher <- diag(solve(res$FisherMatrix))[1:3]
result = data.frame(mean = estimate.mean,
                    std = sqrt(estimate.var),
                    fisher = sqrt(fisher),
                    row.names = c("beta.fixed",
                                  "beta.random1","beta.random2"))
print(result)


K = operator_list$tau*operator_list$Q[[1]]
Sigma.Z = solve(t(K)%*%diag(1/operator_list$h[[1]])%*%K)

F_random <- F_fixed <- 0
for(i in 1:length(locs))
{
  A = spde.A(locs[[i]],operator_list$loc[[1]])
  Q =  B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + A%*%Sigma.Z%*%t(A) + sigma.e^2*diag(n.obs[i])
  F_fixed <-  F_fixed  + t(B_fixed[[i]] )%*%solve(Q,B_fixed[[i]])
  F_random <- F_random + t(B_random[[i]])%*%solve(Q,B_random[[i]])
}
diag(res$FisherMatrix[1:3,1:3])/c(F_fixed[1],diag(F_random))
