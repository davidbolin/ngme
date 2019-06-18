###
# testing estimating fisher for Gaussian random effect model
#
# D: 2019-03-13
###
graphics.off()
library(ngme)
library(testthat)
library(doParallel)
library(Matrix)
library(numDeriv)
library(pracma)
set.seed(1)
use.process = T
estimate.parameters = FALSE
#data option
n.pers <- 10
n.obs  <- rep(10,n.pers)#10 + (1:n.pers)
cutoff = 0.1
max.dist = 1

#Fisher options
niter.est =10
nIter.fisher = 1
nSim.fisher = 1
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
  B_random[[i]] <- cbind(rep(1, n.obs[i]), locs[[i]]/max(locs[[i]]))
  #fixed effects, sqrt(t) and 1/t
  B_fixed[[i]]  <- cbind(sqrt(locs[[i]]),1/locs[[i]])
}
mError_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(0.1,0.2)),
                          beta_fixed  = as.matrix(c(-0.1,0.2)),
                          Sigma =  matrix(c(0.1,0.05,0.05,0.2),ncol=2),
                          noise = "Normal",
                          Sigma_epsilon=0)

mixedEffect_list_in = mixedEffect_list

if(use.process){
  operator_list <- create_operator(locs, name = "fd2",extend = 0,max,max.dist = 1)
  operator_list$type  <- "fd2"
  operator_list$tau   <- 5
  processes_list = list(noise = "Normal",
                        nu  = 0.,
                        mu  = 0.)
  processes_list$V <- list()
  for(i in 1:length(locs))
  {
    processes_list$V[[i]] <- operator_list$h[[1]]
    processes_list$X[[i]] <- 0 * operator_list$h[[1]]
  }
  
}

###
# simulation
#
###
if(use.process){
  sim_res <- simulateLongPrior( Y                 = Y,
                                locs              = locs,
                                mixedEffect_list  = mixedEffect_list,
                                measurment_list   = mError_list,
                                processes_list    = processes_list,
                                operator_list     = operator_list)
  processes_list$X <- sim_res$X
  
} else {
  sim_res <- simulateLongPrior( Y                 = Y,
                                locs              = locs,
                                mixedEffect_list  = mixedEffect_list_in,
                                measurment_list   = mError_list)
}

  
if(estimate.parameters){
  if(use.process){
    res.est <- estimateLong(Y                = sim_res$Y,
                            nIter            = niter.est,
                            nSim             = nSim.fisher,
                            locs             = locs,
                            nBurnin           = nBurnin,
                            mixedEffect_list = mixedEffect_list,
                            processes_list   = processes_list,
                            operator_list    = operator_list,
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
                               processes_list   = res.est$processes_list,
                               operator_list    = res.est$operator_list,
                               nIter = nIter.fisher,
                               nSim             = nSim.fisher,
                               silent = F,
                               nBurnin_base = nBurnin,
                               estimate_fisher = 2)
  } else {
    res.est <- estimateLong(Y                = sim_res$Y,
                            nIter            = niter.est,
                            nSim             = nSim.fisher,
                            locs             = locs,
                            nBurnin           = nBurnin,
                            mixedEffect_list = mixedEffect_list,
                            nBurnin_learningrate = 0,
                            measurment_list  = mError_list,
                            pSubsample = 1,
                            subsample.type = 1,
                            silent = FALSE,
                            estimate_fisher = 0)  
    res.fisher <- estimateLong(Y                = sim_res$Y,
                               locs             = locs,
                               mixedEffect_list = res.est$mixedEffect_list,
                               measurment_list  = res.est$measurementError_list,
                               nIter = nIter.fisher,
                               nSim             = nSim.fisher,
                               silent = T,
                               nBurnin_base = nBurnin,
                               estimate_fisher = 2)
  }
} else {
  if(use.process){
    res.fisher <- estimateLong(Y                = sim_res$Y,
                               locs             = locs,
                               mixedEffect_list = mixedEffect_list,
                               measurment_list  = mError_list,
                               processes_list   = processes_list,
                               operator_list    = operator_list,
                               nIter = nIter.fisher,
                               nSim             = nSim.fisher,
                               silent = T,
                               nBurnin_base = nBurnin,
                               estimate_fisher = 2)  
  } else {
    res.fisher <- estimateLong(Y                = sim_res$Y,
                               locs             = locs,
                               mixedEffect_list = mixedEffect_list,
                               measurment_list  = mError_list,
                               nIter = nIter.fisher,
                               nSim             = nSim.fisher,
                               silent = T,
                               nBurnin_base = nBurnin,
                               estimate_fisher = 2)  
  }
}
  
  F_random <- F_fixed <- 0
  Vgrad_F = matrix(0, 
                   nrow = dim(B_fixed[[1]])[2], 
                   ncol = dim(B_fixed[[1]])[2])
  Vgrad_R = matrix(0, 
                   nrow = dim(B_random[[1]])[2], 
                   ncol = dim(B_random[[1]])[2])
  gradFixed <- Vgrad_F
  F_total <- matrix(0, nrow = dim(B_fixed[[1]])[2] + dim(B_random[[1]])[2],
                       ncol = dim(B_fixed[[1]])[2] + dim(B_random[[1]])[2])
  grad_sigma <- rep(0,4)
  Lik_E <- rep(0,4)
  Lik_Em <- rep(0,4)
  F_numeric <- matrix(0, nrow=2+2+3,ncol=2+2+3)
  for(i in 1:length(locs))
  {
    SigmaE  = mError_list$sigma^2*diag(n.obs[i])
    SigmaEi = solve(SigmaE)
    if(use.process){
      Djoint <- cbind(B_random[[i]], res.fisher$A[[i]])
      D       = B_random[[i]]
      operator_list$Q[[i]] <- as.matrix(operator_list$Q[[i]])
      Sigma <- bdiag(mixedEffect_list$Sigma,solve(operator_list$tau^2*t(operator_list$Q[[i]])%*%diag(1/operator_list$h[[i]])%*%operator_list$Q[[i]]))
    } else {
      Djoint       = B_random[[i]]
      D       = B_random[[i]]
      Sigma   = mixedEffect_list$Sigma
    }
    B       = B_fixed[[i]]  
    Sigma<- as.matrix(Sigma)
    Djoint <- as.matrix(Djoint)
    VU      = solve(solve(Sigma) + t(Djoint)%*%SigmaEi%*%Djoint) #V[U|Y]
    Vgrad_F   = Vgrad_F + t(B)%*%SigmaEi%*%(Djoint%*%VU%*%t(Djoint))%*%SigmaEi%*%B #V[\Delta_\beta\log(\pi)|Y]
    Vgrad_R   = Vgrad_R + t(D)%*%SigmaEi%*%(Djoint%*%VU%*%t(Djoint))%*%SigmaEi%*%D #V[\Delta_\beta\log(\pi)|Y]
    S =  Djoint%*%Sigma%*%t(Djoint)  +  SigmaE
    F_fixed   <- F_fixed   + t(B )%*%solve(S, B)
    F_random  <- F_random  + t(D)%*%solve(S, D)
    gradFixed <- gradFixed +t(B )%*%SigmaEi%*%B
    
    Sigma_total  <- SigmaE + Djoint%*%Sigma%*%t(Djoint)
    iSigma_total <- solve(Sigma_total)
    # numerical differntation for cross gradient
    loglik <- function(betasigma){
     
      beta_fixed  <- c(betasigma[1],betasigma[2])
      beta_random <- c(betasigma[3],betasigma[4])
      Sigma[1,1] <- betasigma[5]
      Sigma[2,2] <- betasigma[6]
      Sigma[2,1] <- betasigma[7]
      Sigma[1,2] <- betasigma[7]
      print(Sigma[1:2,1:2])
      Sigma_total  <- SigmaE + Djoint%*%Sigma%*%t(Djoint)
      iSigma_total <- solve(Sigma_total)
      res <- sim_res$Y[[i]] - cbind(B,D)%*%c(beta_fixed,beta_random)
      Lik <- sum(log(diag(chol(iSigma_total)))) - 0.5*t(res)%*%iSigma_total%*%res
      return(Lik)
    }
   # vecH <- matrix(loglik, c(mixedEffect_list$beta_fixed,
    #                                 mixedEffect_list$beta_random,
    #                                 c(Sigma)[c(1,4,2)]),
    #                       method.args=list(eps=1e-3,d=0.001))
    #Hnum <- -matrix(vecH,nrow=sqrt(length(vecH)))
    Hnum <- -pracma::hessian(loglik,c(mixedEffect_list$beta_fixed,
                                                       mixedEffect_list$beta_random,
                                                       Sigma[1,1],Sigma[2,2],Sigma[1,2]),
                             h = 10^-4)
    F_numeric <-  F_numeric + Hnum
  }
  F_fixed.est <- res.fisher$FisherMatrix[1:2,1:2]
  F_random.est <- res.fisher$FisherMatrix[3:4,3:4]
  Vgrad_R.est <- res.fisher$GradientVariance[3:4,3:4]

test_that("Fisher Gaussian mixed effect", {  
    expect_equal(as.vector(as.matrix(F_random.est)),
                 as.vector(as.matrix(F_random)),
                 tolerance = 10^-4)
    expect_equal(as.vector(F_fixed.est),
                 as.vector(F_fixed),
                 tolerance = 10^-4)
    expect_equal(as.vector(res.fisher$FisherMatrix[1:7,1:7]),
                 as.vector(F_numeric),
                 tolerance = 10^-4)
  #(F_random-Vgrad_R)/length(locs) 
  #Vgrad_R/length(locs)
  })
