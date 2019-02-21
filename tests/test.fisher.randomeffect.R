###
# testing estimating fisher for Gaussian random effect model
#
# D: 2018-03-13
###
graphics.off()
library(ngme)
library(doParallel)
test.fisher = TRUE

#data options
n.pers <- 20
n.obs  <- rep(10,n.pers)#10 + (1:n.pers)
cutoff = 0.1
max.dist = 1

#Fisher options
nIter.fisher = 500
nSim.fisher = 300
nBurnin = 10

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
  locs[[i]] <- 1:n.obs[i] #sort(1 + 9*runif(n.obs[i]))
  #random effects, 1 and t
  B_random[[i]] <- cbind(rep(1, n.obs[i]), locs[[i]])
  #fixed effects, sqrt(t) and 1/t
  B_fixed[[i]]  <- cbind(sqrt(locs[[i]]),1/locs[[i]])
}
mError_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(0.1,0.2)),
                          beta_fixed  = as.matrix(c(0.1,0.2)),
                          Sigma = diag(c(0.1, 0.2)),
                          noise = "Normal",
                          Sigma_epsilon=0)





###
# simulation
#
###
mixedEffect_list_in = mixedEffect_list

sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list_in,
                              measurment_list   = mError_list)




if(test.fisher){
  
  res.est <- estimateLong(Y                = sim_res$Y,
                          nIter            = nIter.fisher,
                          nSim             = nSim.fisher,
                          locs             = locs,
                          nBurnin           = nBurnin,
                          mixedEffect_list = mixedEffect_list,
                          nBurnin_learningrate = 0,
                          measurment_list  = mError_list,
                          pSubsample = 1,
                          subsample.type = 1,
                          silent = TRUE,
                          estimate_fisher = 0)
 
  res.fisher <- estimateLong(Y                = sim_res$Y,
                           locs             = locs,
                           mixedEffect_list = res.est$mixedEffect_list,
                           measurment_list  = res.est$measurementError_list,
                            nIter = nIter.fisher,
                           nSim             = nSim.fisher,
                            silent = FALSE,
                            nBurnin_base = nBurnin,
                            estimate_fisher = 2)
  F_random <- F_fixed <- 0
  Vgrad_F = matrix(0, 
                   nrow = dim(B_fixed[[1]])[2], 
                   ncol = dim(B_fixed[[1]])[2])
  Vgrad_R = matrix(0, 
                   nrow = dim(B_random[[1]])[2], 
                   ncol = dim(B_random[[1]])[2])
  gradFixed <- Vgrad_F
  for(i in 1:length(locs))
  {
    D       = B_random[[i]]
    B       = B_fixed[[i]]
    Sigma   = mixedEffect_list$Sigma
    SigmaE  = mError_list$sigma^2*diag(n.obs[i])
    SigmaEi = solve(SigmaE)
    VU      = solve(solve(Sigma) + t(D)%*%SigmaEi%*%D) #V[U|Y]
    Vgrad_F   = Vgrad_F + t(B)%*%SigmaEi%*%D%*%VU%*%t(D)%*%SigmaEi%*%B #V[\Delta_\beta\log(\pi)|Y]
    Vgrad_R   = Vgrad_R + t(D)%*%SigmaEi%*%D%*%VU%*%t(D)%*%SigmaEi%*%D #V[\Delta_\beta\log(\pi)|Y]
    S =  D%*%Sigma%*%t(D)  +  SigmaE
    F_fixed   <- F_fixed   + t(B )%*%solve(S, B)
    F_random  <- F_random  + t(D)%*%solve(S, D)
    gradFixed <- gradFixed +t(B )%*%SigmaEi%*%B
  }
  F_fixed.est <- res.fisher$FisherMatrix[1:2,1:2]
  F_random.est <- res.fisher$FisherMatrix[3:4,3:4]
  Vgrad_R.est <- res.fisher$GradientVariance[3:4,3:4]
  cat("random:\n")
  print(F_random.est/F_random)
  cat("inverse random:\n")
  print(solve(F_random.est)/solve(F_random))
  cat("fixed:\n")
  print(F_fixed.est/F_fixed)
  cat("inverse fixed:\n")
  print(solve(F_fixed.est)/solve(F_fixed))
  cat('variance estimate:\n')
  print(res.fisher$GradientVariance[1:2,1:2]/Vgrad_F)
  #(F_random-Vgrad_R)/length(locs) 
  #Vgrad_R/length(locs)
}
