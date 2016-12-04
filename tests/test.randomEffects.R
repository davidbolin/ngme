graphics.off()
library(LDMod)

test.pred = TRUE
n.threads <- 1
nIter <- 1000
n.pers <- 1000
nSim  <- 2
nSim.pred <- 200
nBurnin.pred <- 100
n.obs  <- 40 + 0*(1:n.pers)
n <- 100
error.dist = "NIG"
mixed.dist = "NIG"
nBurnin = 50
pSubsample = 0.1
nPar_burnin = 100
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()
Vin <- list()


for(i in 1:n.pers)
{
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))

  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  Vin[[i]] <- rep(1, n.obs[i])

  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
}
mError_list <- list(Vs = Vin, noise = error.dist, sigma = 0.1, nu = 1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(0.1,0.2)),
                          beta_fixed  = as.matrix(c(0.1)),
                          Sigma = diag(c(0.1, 0.2)),
                          noise = mixed.dist,
                          Sigma_epsilon=1,
                          nu = 1,
                          mu = matrix(c(2,2),2,1))


sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = mError_list)


res.est <- estimate.wrapper(Y = sim_res$Y,
                              locs = locs,
                              B_random= B_random,
                              B_fixed = B_fixed,
                              use.process = FALSE,
                              measurement.distribution = error.dist,
                              random.effect.distribution = mixed.dist,
                              estimation.options = list(nIter.gauss = 10,nIter = nIter,
                                                        pSubsample = pSubsample,
                                                        nPar_burnin = nPar_burnin))



par(mfrow = c(2,3))
matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1)
matplot(res.est$mixedEffect_list$betar_vec,type="l",main="random effects",col=1)
matplot(res.est$mixedEffect_list$Sigma_vec,type="l",main="RE Sigma",col=1)
if(error.dist == "NIG"){
  plot(res.est$measurementError_list$nu_vec,type="l",main="error nu")
}
if(mixed.dist == "NIG"){
  plot(res.est$mixedEffect_list$nu_vec,type="l",main="mixed nu")
  matplot(res.est$mixedEffect_list$mu_vec,type="l",main="RE mu",col=1)
}

truth =  c(mixedEffect_list$beta_fixed,
           mixedEffect_list$beta_random,
           mixedEffect_list$Sigma[c(1,2,4)],
           mError_list$sigma)

start.values  = c(res.est$mixedEffect_list$betaf_vec[1],
                  res.est$mixedEffect_list$betar_vec[1,],
                  res.est$mixedEffect_list$Sigma_vec[1,c(1,2,4)],
                  res.est$measurementError_list$sigma_vec[1])

estimates =  c(res.est$mixedEffect_list$betaf_vec[nIter],
               res.est$mixedEffect_list$betar_vec[nIter,],
               res.est$mixedEffect_list$Sigma_vec[nIter,c(1,2,4)],
               res.est$measurementError_list$sigma_vec[nIter])

row.names = c("beta.fixed",
              "beta.random1","beta.random2",
              "Sigma1", "Sigma2","Sigma3",
              "sigma")

if(error.dist == "NIG"){
  truth = c(truth,mError_list$nu)
  start.values = c(start.values,res.est$measurementError_list$nu_vec[1])
  estimates = c(estimates, res.est$measurementError_list$nu_vec[nIter])
  row.names = c(row.names,"nu.error")
}
if(mixed.dist == "NIG"){
  truth = c(truth,mixedEffect_list$nu, mixedEffect_list$mu)
  start.values= c(start.values, res.est$mixedEffect_list$nu_vec[1], res.est$mixedEffect_list$mu_vec[1,])
  estimates = c(estimates, res.est$mixedEffect_list$nu_vec[nIter], res.est$mixedEffect_list$mu_vec[nIter,])
  row.names = c(row.names,"nu.mixed", "mu.mixed1","mu.mixed2")
}
result = data.frame(start = start.values,
                    estimate = estimates,
                    true = truth,
                    row.names = row.names)
print(result)

if(test.pred){
  res <- predictLong( Y = sim_res$Y,
                      pInd = c(1,2),
                      locs.pred = locs,
                      Brandom.pred = B_random,
                      Bfixed.pred = B_fixed,
                      type = "Filter",
                      nSim             = nSim.pred,
                      nBurnin = nBurnin.pred,
                      locs             = locs,
                      mixedEffect_list = mixedEffect_list,
                      measurment_list  = mError_list,
                      return.samples = TRUE,
                      quantiles = c(0.05,0.95),
                      max.num.threads = n.threads)


  par(mfrow = c(1,2))
  for(k in c(1,2)){
    plot(locs[[k]],sim_res$Y[[k]],ylim=c(min(res$X.summary[[k]]$quantiles[[1]]$field),
                                         max(res$X.summary[[k]]$quantiles[[2]]$field)))
    lines(res$locs[[k]],res$X.summary[[k]]$Mean)
    lines(res$locs[[k]],res$Y.summary[[k]]$quantiles[[1]]$field,col=2)
    lines(res$locs[[k]],res$Y.summary[[k]]$quantiles[[2]]$field,col=2)
  }


}
