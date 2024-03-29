n.obs  <- 5
n <- 10
mix.dist = "NIG"
error.dist = "NIG"
pred.type <- "Filter"
nSim <- 1000
Y <- locs <- B_random <- B_fixed  <- X <- V <- Vin <- list()
for(i in 1:5)
{
  B_random[[i]] <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  Y[[i]] <- rep(0,n.obs)
  Vin[[i]] <- rep(1,n.obs)
  Y[[i]][i] = 1
  locs[[i]] <- (1:n.obs)*i
  B_fixed[[i]]  <- as.matrix(locs[[i]])
  X[[i]] <- rep(0, n)
  V[[i]] <- rep(1, n)
}

mError_list <- list(noise = error.dist, sigma = 0.1, nu = 1,common_V = FALSE,Vs = Vin)
mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = as.matrix(c(2,-1)),
                            beta_fixed  = as.matrix(c(.1)),
                            Sigma = diag(c(0.1, 0.2)),
                            noise = mix.dist,
                            Sigma_epsilon=1)


if(mix.dist == "NIG"){
  mixedEffect_list$nu = as.matrix(c(1))
  mixedEffect_list$mu = as.matrix(c(1,1))
}


operator_list <- create_operator(locs, n, name = "Matern")
operator_list$kappa <- 2
operator_list$tau   <- 3


processes_list = list(noise = "Normal")
processes_list$X = X
processes_list$V = V



res <- predictLong( Y = Y,
                    locs = locs,
                    Brandom.pred = B_random,
                    Bfixed.pred = B_fixed,
                    type = pred.type,
                    nSim = nSim,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list  = mError_list,
                    processes_list = processes_list,
                    operator_list = operator_list,
                    max.num.threads = 1)

res <- predictLong( Y = Y,
                    locs = locs,
                    Brandom.pred = B_random,
                    Bfixed.pred = B_fixed,
                    type = pred.type,
                    nSim = nSim,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list  = mError_list,
                    processes_list = processes_list,
                    operator_list = operator_list,
                    max.num.threads = 2)

res2 <- predictLong(pInd = 2, Y = Y,locs = locs, Brandom.pred = B_random, Bfixed.pred = B_fixed, type = pred.type,
                    nSim = nSim, mixedEffect_list = mixedEffect_list, measurment_list  = mError_list,
                    processes_list = processes_list, operator_list = operator_list)

res$X.summary[[2]]$Mean  - res2$X.summary[[1]]$Mean
