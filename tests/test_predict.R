graphics.off()
library(ngme)

n.threads <- 1
nIter <- 100
n.pers <- 2
nSim  <- 100
use.random.effect = FALSE
n.obs  <- 3 + 0*(1:n.pers)
n.proc = n.obs + 2
grid.extend = c(0,0.1)
n <- 10
n.pred <- 15
nBurnin = 10
pred.type <- "Filter"
pSubsample = 0.99
#subsample.type = 2
test.pred = TRUE
Y <- list()
locs <- list()
B_random <- B_fixed  <- locs.pred <- locs.proc <- list()

list()
B_random.pred <- list()
B_fixed.pred <- list()

delta = 0.25
B_random.pred1 <- B_fixed.pred1 <- list()

Vin <- list()
for(i in 1:n.pers)
{
  B_random[[i]] <- cbind(rep(1, n.obs[i]), (1:n.obs[i]) / n.obs[i] )
  B_random.pred[[i]] <- cbind(rep(1, n.pred), (1:n.pred) / n.pred )
  B_random.pred1[[i]] <- cbind(rep(1, n.pred), delta + (1:n.pred) / n.pred )
  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  locs.proc[[i]] <- 0:n.obs[i]
  locs.pred[[i]] <- seq(from = 1, to = n.obs[i], length.out = n.pred)
  Vin[[i]] <- rep(1, n.obs[i])

  B_fixed[[i]]  <- as.matrix(locs[[i]])
  B_fixed.pred[[i]]  <- as.matrix(locs.pred[[i]])
  B_fixed.pred1[[i]]  <- as.matrix(delta + locs.pred[[i]])
}

if(use.random.effect){
  derivative_list = list(Bfixed = B_fixed.pred1,
                         Brandom = B_random.pred1,
                         delta = delta)
} else {
  derivative_list = list(Bfixed = B_fixed.pred1,
                         delta = delta)
}

mError_list <- list(Vs = Vin, noise = "NIG", sigma = 0.01, nu = 1)
if(use.random.effect){
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = as.matrix(c(2,-1)),
                            beta_fixed  = as.matrix(c(.1)),
                            Sigma = diag(c(0.1, 0.2)),
                            noise = "Normal",
                            Sigma_epsilon=1)

} else {
  mixedEffect_list  <- list(B_fixed  = B_fixed,
                            beta_fixed  = as.matrix(c(.1)),
                            noise = "Normal")

}

operator_list <- create_operator(locs.proc, n, name = "fd2",extend = grid.extend)
operator_list$type  <- "fd2"
operator_list$tau   <- 5


processes_list = list(noise = "Normal",
                      nu  = 0.,
                      mu  = 0.)
processes_list$V <- list()
for(i in 1:length(locs))
{
  processes_list$V[[i]] <- operator_list$h[[1]]
}


sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = mError_list,
                              processes_list    = processes_list,
                              operator_list     = operator_list)

processes_list$X <- sim_res$X
prediction.indices= c(1)
if(use.random.effect){
  res.pre0 <- predictLong( Y = sim_res$Y,
                           pInd = prediction.indices,
                           locs.pred = locs.pred,
                           Brandom.pred = B_random.pred,
                           Bfixed.pred = B_fixed.pred,
                           type = pred.type,
                           nSim             = nSim,
                           locs             = locs,
                           mixedEffect_list = mixedEffect_list,
                           predict.derivatives = derivative_list,
                           measurment_list  = mError_list,
                           processes_list   = processes_list,
                           operator_list    = operator_list,
                           excursions       = list(list(type = '<', level = -0.05, process = 'Xderivative')),
                           return.samples = TRUE,
                           quantiles = c(0.05,0.95),
                           max.num.threads = n.threads)

} else {
  res.pre0 <- predictLong( Y = sim_res$Y,
                           pInd = prediction.indices,
                           locs.pred = locs.pred,
                           Bfixed.pred = B_fixed.pred,
                           type = pred.type,
                           nSim             = nSim,
                           locs             = locs,
                           mixedEffect_list = mixedEffect_list,
                           measurment_list  = mError_list,
                           processes_list   = processes_list,
                           operator_list    = operator_list,
                           predict.derivatives = derivative_list,
                           excursions       = list(list(type = '<', level = -0.05, process = 'Xderivative')),
                           return.samples = TRUE,
                           quantiles = c(0.05,0.95),
                           max.num.threads = n.threads,
                           seed = 1)

}

k = 1
locs.pred <-  list(seq(from = 1, to = n.obs[k], length.out = n.pred))
B_random.pred  <- list(cbind(rep(1, n.pred), (1:n.pred) / n.pred ))
B_fixed.pred  <- list(as.matrix(locs.pred[[1]]))
B_random.pred1  <- list(cbind(rep(1, n.pred), delta + (1:n.pred) / n.pred ))
B_fixed.pred1  <- list(as.matrix(delta + locs.pred[[1]]))

derivative_list = list(Bfixed = B_fixed.pred1,
                       Brandom = B_random.pred1,
                       delta = delta)

if(use.random.effect){
  pred.lists <- updateLists(mixedEffect_list = mixedEffect_list,
                            processes_list = processes_list,
                            operator_list = operator_list,
                            measurement_list = mError_list,
                            Bfixed = B_fixed[prediction.indices],
                            Brandom = B_random[prediction.indices],
                            locs = locs[prediction.indices])

    res.pre <- predictLong(Y = sim_res$Y[k],
                         locs.pred        = locs.pred,
                         locs             = locs[k],
                         Brandom.pred = B_random.pred,
                         Bfixed.pred = B_fixed.pred,
                         type             = "Filter",
                         excursions       = list(list(type = '<', level = -0.05, process = 'Xderivative')),
                         mixedEffect_list = pred.lists$mixedEffect_list,
                         measurment_list  = pred.lists$measurement_list,
                         processes_list   = pred.lists$processes_list,
                         predict.derivatives = derivative_list,
                         operator_list    = operator_list,
                         nSim             = nSim,
                         max.num.threads = 1,
                         quantiles = c(0.025,0.975),
                         return.samples = TRUE)

} else {
  pred.lists <- updateLists(mixedEffect_list = mixedEffect_list,
                            processes_list = processes_list,
                            operator_list = operator_list,
                            measurement_list = mError_list,
                            Bfixed = B_fixed[prediction.indices],
                            locs = locs[prediction.indices])


  res.pre <- predictLong(Y = sim_res$Y[k],
                         locs.pred        = locs.pred,
                         locs             = locs[k],
                         Bfixed.pred = B_fixed.pred,
                         type             = "Filter",
                         excursions       = list(list(type = '<', level = -0.05, process = 'Xderivative')),
                         mixedEffect_list = pred.lists$mixedEffect_list,
                         measurment_list  = pred.lists$measurement_list,
                         predict.derivatives = derivative_list,
                         processes_list   = pred.lists$processes_list,
                         operator_list    = operator_list,
                         nSim             = nSim,
                         max.num.threads = 1,
                         quantiles = c(0.025,0.975),
                         return.samples = TRUE)

}


n.pred = 2*n.pred-1
k = 1
locs.pred <- list(seq(from = 1, to = n.obs[k], length.out = n.pred))
B_random.pred  <- list(cbind(rep(1, n.pred), (1:n.pred) / n.pred ))
B_fixed.pred  <- list(as.matrix(locs.pred[[1]]))

B_random.pred1  <- list(cbind(rep(1, n.pred), delta + (1:n.pred) / n.pred ))
B_fixed.pred1  <- list(as.matrix(delta + locs.pred[[1]]))

derivative_list = list(Bfixed = B_fixed.pred1,
                       Brandom = B_random.pred1,
                       delta = delta)

if(use.random.effect){
  pred.lists <- updateLists(mixedEffect_list = mixedEffect_list,
                            processes_list = processes_list,
                            operator_list = operator_list,
                            measurement_list = mError_list,
                            Bfixed = B_fixed[prediction.indices],
                            Brandom = B_random[prediction.indices],
                            locs = locs[prediction.indices])

  res.pre2 <- predictLong(Y = sim_res$Y[k],
                         locs.pred        = locs.pred,
                         locs             = locs[k],
                         Brandom.pred = B_random.pred,
                         Bfixed.pred = B_fixed.pred,
                         type             = "Filter",
                         excursions       = list(list(type = '<', level = -0.05, process = 'Xderivative')),
                         mixedEffect_list = pred.lists$mixedEffect_list,
                         measurment_list  = pred.lists$measurement_list,
                         predict.derivatives = derivative_list,
                         processes_list   = pred.lists$processes_list,
                         operator_list    = operator_list,
                         nSim             = nSim,
                         max.num.threads = 1,
                         quantiles = c(0.025,0.975),
                         return.samples = TRUE)

} else {
  pred.lists <- updateLists(mixedEffect_list = mixedEffect_list,
                            processes_list = processes_list,
                            operator_list = operator_list,
                            measurement_list = mError_list,
                            Bfixed = B_fixed[prediction.indices],
                            locs = locs[prediction.indices])


  res.pre2 <- predictLong(Y = sim_res$Y[k],
                         locs.pred        = locs.pred,
                         locs             = locs[k],
                         Bfixed.pred = B_fixed.pred,
                         type             = "Filter",
                         excursions       = list(list(type = '<', level = -0.05, process = 'Xderivative')),
                         mixedEffect_list = pred.lists$mixedEffect_list,
                         measurment_list  = pred.lists$measurement_list,
                         processes_list   = pred.lists$processes_list,
                         predict.derivatives = derivative_list,
                         operator_list    = operator_list,
                         nSim             = nSim,
                         max.num.threads = 1,
                         quantiles = c(0.025,0.975),
                         return.samples = TRUE)

}


par(mfrow = c(1,3))
pr <- c(min(min(res.pre$X.summary[[k]]$quantiles[[1]]$field),min(Y[[prediction.indices[k]]])),
        max(max(res.pre$X.summary[[k]]$quantiles[[2]]$field),max(Y[[prediction.indices[k]]])))
plot(res.pre$locs[[k]],res.pre$X.summary[[k]]$Mean,type="l",ylim=pr,
     xlab = "Follow-up time (in years)",ylab="log(eGFR)")
points(res.pre$locs[[k]],res.pre$X.summary[[k]]$Mean,pch = 4)
#lines(res.pre$locs[[k]],res.pre$X.summary[[k]]$quantiles[[1]]$field,col=2)
#lines(res.pre$locs[[k]],res.pre$X.summary[[k]]$quantiles[[2]]$field,col=2)
points(locs[[prediction.indices[k]]],sim_res$Y[[prediction.indices[k]]],col=2)


lines(res.pre2$locs[[k]],res.pre2$X.summary[[k]]$Mean,col=1,lty=2)
points(res.pre2$locs[[k]],res.pre2$X.summary[[k]]$Mean,pch = 5)
#lines(res.pre2$locs[[k]],res.pre2$X.summary[[k]]$quantiles[[1]]$field,col=2,lty=2)
#lines(res.pre2$locs[[k]],res.pre2$X.summary[[k]]$quantiles[[2]]$field,col=2,lty=2)


lines(res.pre0$locs[[k]],res.pre0$X.summary[[k]]$Mean,col=1,lty=3)
points(res.pre0$locs[[k]],res.pre0$X.summary[[k]]$Mean,pch = 6)
#lines(res.pre0$locs[[k]],res.pre0$X.summary[[k]]$quantiles[[1]]$field,col=2,lty=3)
#lines(res.pre0$locs[[k]],res.pre0$X.summary[[k]]$quantiles[[2]]$field,col=2,lty=3)

pr <- c(min(res.pre$Xderivative.summary[[k]]$Mean),
        max(res.pre$Xderivative.summary[[k]]$Mean))
plot(res.pre$locs[[k]],res.pre$Xderivative.summary[[k]]$Mean,type="l",ylim=pr,
     xlab = "Follow-up time (in years)",ylab="log(eGFR)")
points(res.pre$locs[[k]],res.pre$Xderivative.summary[[k]]$Mean,pch = 4)
#lines(res.pre$locs[[k]],res.pre$X.summary[[k]]$quantiles[[1]]$field,col=2)
#lines(res.pre$locs[[k]],res.pre$X.summary[[k]]$quantiles[[2]]$field,col=2)

lines(res.pre2$locs[[k]],res.pre2$Xderivative.summary[[k]]$Mean,col=1,lty=2)
points(res.pre2$locs[[k]],res.pre2$Xderivative.summary[[k]]$Mean,pch = 5)
#lines(res.pre2$locs[[k]],res.pre2$X.summary[[k]]$quantiles[[1]]$field,col=2,lty=2)
#lines(res.pre2$locs[[k]],res.pre2$X.summary[[k]]$quantiles[[2]]$field,col=2,lty=2)


lines(res.pre0$locs[[k]],res.pre0$Xderivative.summary[[k]]$Mean,col=1,lty=3)
points(res.pre0$locs[[k]],res.pre0$Xderivative.summary[[k]]$Mean,pch = 6)


cr = c(0,max(max(res.pre0$Xderivative.summary[[k]]$excursions$P),max(res.pre2$Xderivative.summary[[k]]$excursions$P)))
plot(res.pre0$locs[[k]],res.pre0$Xderivative.summary[[k]]$excursions$P,type="l",
     xlab = "time",ylab="Probability",
     ylim = cr)

lines(res.pre$locs[[k]],res.pre$Xderivative.summary[[k]]$excursions$P,type="l",
      xlab = "time",ylab="Probability",lty=2)

lines(res.pre2$locs[[k]],res.pre2$Xderivative.summary[[k]]$excursions$P,type="l",
     xlab = "time",ylab="Probability",lty=2)



