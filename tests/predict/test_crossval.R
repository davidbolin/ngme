graphics.off()
library(ngme)

n.threads <- 1
n.pers <- 2
nSim  <- 10
use.random.effect = FALSE
n.obs  <- 10 + 0*(1:n.pers)
n.proc = n.obs + 2
grid.extend = c(0,0.1)
n <- 10
n.pred <- 15
nBurnin = 10
pred.type <- "LOOCV"
pSubsample = 0.99
#subsample.type = 2
test.pred = TRUE
Y <- list()
locs <- list()
B_random <- B_fixed  <- list()


locs.proc <- Vin <- list()
for(i in 1:n.pers)
{
  B_random[[i]] <- cbind(rep(1, n.obs[i]), (1:n.obs[i]) / n.obs[i] )
  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  locs.proc[[i]] <- 0:n.obs[i]
  Vin[[i]] <- rep(1, n.obs[i])
  
  B_fixed[[i]]  <- as.matrix(locs[[i]])
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

operator_list <- create_operator(locs.proc, name = "exponential",extend = grid.extend,max.dist = 1)
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
  res.pre <- predictLong( Y = sim_res$Y,
                          pInd = prediction.indices,
                          type = pred.type,
                          Bfixed.pred = B_fixed,
                          Brandom.pred = B_random,
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
  res.pre <- predictLong( Y = sim_res$Y,
                           pInd = prediction.indices,
                           Bfixed.pred = B_fixed,
                           type = pred.type,
                           nSim             = nSim,
                           locs             = locs,
                           mixedEffect_list = mixedEffect_list,
                           measurment_list  = mError_list,
                           processes_list   = processes_list,
                           operator_list    = operator_list,
                           return.samples = TRUE,
                           quantiles = c(0.05,0.95),
                           max.num.threads = n.threads,
                           seed = 1)
  
}

k=1
pr <- c(min(min(res.pre$X.summary[[k]]$Mean),min(sim_res$Y[[prediction.indices[k]]])),
        1.1*max(max(res.pre$X.summary[[k]]$Mean),max(sim_res$Y[[prediction.indices[k]]])))
plot(res.pre$locs[[k]],res.pre$X.summary[[k]]$Mean,type="l",ylim=pr)
points(res.pre$locs[[k]],res.pre$X.summary[[k]]$Mean)
points(res.pre$locs[[k]],sim_res$Ystar[[k]],col=3)
points(locs[[prediction.indices[k]]],sim_res$Y[[prediction.indices[k]]],col=2)

