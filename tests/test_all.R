graphics.off()
library(ngme)

test.pred = FALSE
test.est = FALSE
test.fisher = TRUE

#data options
n.pers <- 20
n.obs  <- rep(5,n.pers)#10 + (1:n.pers)
cutoff = 0.1
max.dist = 1
operator.type = "matern"
process.noise = "NIG"
process.nu = 0.1
process.mu = 0.1

#estimation options
nIter <- 1000
nSim <- 5
nBurnin = 10
pSubsample = 0.5
subsample.type = 1

#prediction options
n.pred <- n.obs[[1]]#10
pred.type <- "Filter"
nSim.pred  <- 100 #simulations in prediction

#Fisher options
nIter.fisher = 1
nSim.fisher = 1000

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
  locs[[i]] <- sort(1 + 9*runif(n.obs[i]))
  #locs.pred[[i]] <-c(locs[[i]][1])
  locs.pred[[i]] <- unique(sort(c(locs[[i]],seq(from = 0.8*locs[[i]][1], to = 1.2*locs[[i]][n.obs[i]], length.out = n.pred))))
  #locs.pred[[i]] <- locs[[i]]
  #locs.pred[[i]] <-seq(from = locs[[i]][1], to = locs[[i]][n.obs[i]], length.out = n.pred)
  n.pred <- length(locs.pred[[i]])
  Vin[[i]] <- rep(1, n.obs[i])

  #random effects, 1 and t
  B_random[[i]] <- cbind(rep(1, n.obs[i]), locs[[i]])
  B_random.pred[[i]] <- cbind(rep(1, n.pred), locs.pred[[i]])

  #fixed effects, sqrt(t) and 1/t
  B_fixed[[i]]  <- cbind(sqrt(locs[[i]]),1/locs[[i]])
  B_fixed.pred[[i]]  <- cbind(sqrt(locs.pred[[i]]),1/locs.pred[[i]])

}
mError_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1, nu = 1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(0.1,0.2)),
                          beta_fixed  = as.matrix(c(0.1,0.2)),
                          Sigma = diag(c(0.1, 0.2)),
                          noise = "Normal",
                          Sigma_epsilon=1)


operator_list <- create_operator(locs, max.dist=max.dist,cutoff = cutoff, name = operator.type,extend=0.5)
if(operator.type == "matern"){
  operator_list$kappa <- 2
}

operator_list$tau   <- 1



processes_list = list(noise = process.noise,
                      nu  = process.nu,
                      mu  = process.mu)
processes_list$V <- list()
for(i in 1:length(locs))
{
  processes_list$V[[i]] <- operator_list$h[[i]]
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

if(test.est){
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
                          subsample.type = subsample.type,
                          silent = FALSE)


  n.plots <- 3
  if(operator.type == "matern")
    n.plots = 4
  if(process.noise == "NIG")
    n.plots <- n.plots + 2
  if(n.plots <= 3) {
    par(mfrow = c(1,3))
  } else if(n.plots <=4){
    par(mfrow = c(2,2))
  } else {
    par(mfrow = c(2,3))
  }
  matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1)
  matplot(t(matrix(rep(mixedEffect_list$beta_fixed,nIter),nrow=length(mixedEffect_list$beta_fixed),ncol = nIter)),add=TRUE,col=2,type="l")
  matplot(res.est$mixedEffect_list$betar_vec,type="l",main="random effects",col=1)
  matplot(t(matrix(rep(mixedEffect_list$beta_random,nIter),nrow=length(mixedEffect_list$beta_random),ncol = nIter)),add=TRUE,col=2,type="l")

  plot(res.est$operator_list$tauVec,type="l",main="process tau")
  lines(rep(operator_list$tau,nIter),col=2)
  if(operator.type == "matern"){
    plot(res.est$operator_list$kappaVec,type="l",main="process kappa")
    lines(rep(operator_list$kappa,nIter),col=2)
  }
  if(process.noise == "NIG"){
    plot(1:nIter,res.est$processes_list$nu_vec,type="l",main="process nu")
    lines(rep(processes_list$nu,nIter),col=2)
    plot(1:nIter,res.est$processes_list$mu_vec,type="l",main="process mu")
    lines(rep(processes_list$mu,nIter),col=2)
  }
}


if(test.pred){
  res <- predictLong( Y = sim_res$Y,
                      pInd = c(1,2),
                      locs.pred = locs.pred,
                      Brandom.pred = B_random.pred,
                      Bfixed.pred = B_fixed.pred,
                      type = pred.type,
                      nSim             = nSim.pred,
                      locs             = locs,
                      mixedEffect_list = mixedEffect_list,
                      measurment_list  = mError_list,
                      processes_list   = processes_list,
                      operator_list    = operator_list,
                      return.samples = FALSE,
                      quantiles = c(0.05,0.95))


  #x11()
  par(mfrow=c(1,2))
  for(k in 1:2){
    plot(locs[[k]],sim_res$Y[[k]],
         ylim=c(min(res$X.summary[[k]]$quantiles[[1]]$field),max(res$X.summary[[k]]$quantiles[[2]]$field)),
         xlim = c(min(min(locs[[k]]),min(res$locs[[k]])),max(max(locs[[k]]),max(res$locs[[k]])))
    )
    lines(res$locs[[k]],res$X.summary[[k]]$Mean)
    points(res$locs[[k]],res$X.summary[[k]]$Mean,pch=4)
    lines(locs[[k]],sim_res$Y[[k]]-sim_res$E[[k]],col=3)
    lines(res$locs[[k]],res$X.summary[[k]]$quantiles[[1]]$field,col=2)
    lines(res$locs[[k]],res$X.summary[[k]]$quantiles[[2]]$field,col=2)
  }
}

if(test.fisher){

  res.est <- estimateLong(Y                = sim_res$Y,
                          nIter            = nIter.fisher,
                          nSim             = nSim.fisher,
                          locs             = locs,
                          nBurnin           = nBurnin,
                          mixedEffect_list = mixedEffect_list,
                          nBurnin_learningrate = 0,
                          measurment_list  = mError_list,
                          processes_list   = processes_list,
                          operator_list    = operator_list,
                          pSubsample = 1,
                          subsample.type = 1,
                          silent = FALSE,
                          estimate_fisher = 2)


  F_random <- F_fixed <- 0
  for(i in 1:length(locs))
  {
    K = 0.5*operator_list$tau*((operator_list$kappa)^(-1.5)*operator_list$G[[i]] + operator_list$C[[i]]*(operator_list$kappa)^(0.5))
    Sigma.Z = solve(t(K)%*%operator_list$Ci[[i]]%*%K)

    A = spde.A(locs[[i]],operator_list$loc[[i]])
    Q =  B_random[[i]]%*%mixedEffect_list$Sigma%*%t(B_random[[i]]) + A%*%Sigma.Z%*%t(A) + mError_list$sigma^2*diag(n.obs[i])
    F_fixed <-  F_fixed  + t(B_fixed[[i]] )%*%solve(Q,B_fixed[[i]])
    F_random <- F_random + t(B_random[[i]])%*%solve(Q,B_random[[i]])
  }
  F_fixed.est <- res.est$FisherMatrix[1:2,1:2]
  F_random.est <- res.est$FisherMatrix[3:4,3:4]

  cat("random:\n")
  print(F_random.est/F_random)
  cat("fixed:\n")
  print(F_fixed.est/F_fixed)


}
