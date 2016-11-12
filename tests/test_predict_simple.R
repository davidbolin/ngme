graphics.off()
library(LDMod)

n.threads <- 1
nIter <- 1000
n.pers <- 2
nSim  <- 1000
use.random.effect = FALSE
n.obs  <- 3 + 0*(1:n.pers)
grid.extend = c(0,0.1)
n <- 10
n.pred <- 5
nBurnin = 10
pred.type <- "Filter"
pSubsample = 0.99
#subsample.type = 2
test.pred = TRUE
Y <- list()
locs <- list()
B_random <- B_fixed  <- locs.pred <- list()

list()
B_random.pred <- list()
B_fixed.pred <- list()

delta = 0.1
B_random.pred1 <- B_fixed.pred1 <- list()

Vin <- list()
for(i in 1:n.pers)
{
  Y[[i]] <- rep(1,n.obs[i])
  locs[[i]] <- 1:n.obs[i]
  locs.pred[[i]] <- seq(from = 1, to = n.obs[i], length.out = n.pred)
  Vin[[i]] <- rep(1, n.obs[i])
  B_fixed[[i]]  <- matrix(rep(1,n.obs[i]))
  B_fixed.pred[[i]]  <- matrix(rep(1,n.pred))
  B_fixed.pred1[[i]]  <- matrix(rep(1,n.pred))
}

derivative_list = list(Bfixed = B_fixed.pred1,
                         delta = delta)

mError_list <- list(Vs = Vin, noise = "Normal", sigma = 0.01, nu = 1)

mixedEffect_list  <- list(B_fixed  = B_fixed,
                            beta_fixed  = as.matrix(c(1)),
                            noise = "Normal")

operator_list <- create_operator(locs, n, name = "fd2",extend = grid.extend)
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
res.pre0 <- predictLong( Y = sim_res$Y,
                           pInd = prediction.indices,
                           locs.pred = locs.pred,
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


k = 1
pr <- c(min(min(res.pre0$X.summary[[k]]$quantiles[[1]]$field),min(Y[[prediction.indices[k]]])),
        max(max(res.pre0$X.summary[[k]]$quantiles[[2]]$field),max(Y[[prediction.indices[k]]])))
plot(res.pre0$locs[[k]],res.pre0$X.summary[[k]]$Mean,type="l",ylim=pr,
     xlab = "Follow-up time (in years)",ylab="log(eGFR)")
points(res.pre0$locs[[k]],res.pre0$X.summary[[k]]$Mean,pch = 4)
points(locs[[prediction.indices[k]]],sim_res$Y[[prediction.indices[k]]],col=2)
