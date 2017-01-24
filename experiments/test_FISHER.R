#TODO: for processes, create idenity operator, then compute
#      the actual integral
rm(list=ls())
library(LDMod)
library(MASS)
seed     <- 9
silent   <- 1
plotflag <- 1

nIter <- 1
n.pers <- 1
nSim  <- 2
n.obs  <- 10 + 0*(1:n.pers)

nBurnin = 10
pSubsample = 0.1
nPar_burnin = 0
n = 20
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()

betaf <- as.matrix(c(2.1))
betar = as.matrix(c(0.8,0.5))
tau <- 75
Sigma <- matrix(c(0.02,0,0,0.01),2,2)

sd_Y = 0.2
set.seed(seed)


for(i in 1:n.pers)
{
  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))
  locs[[i]] <- 1:n.obs[i]
  Y[[i]] <- rep(1,n.obs[i])
}


operator_list <- create_operator(locs, n, name = "fd2")
operator_list$tau   <- tau
K = operator_list$tau*operator_list$Q[[1]]
Q_op = t(K)%*%diag(1/operator_list$h[[1]])%*%K
Sigma.Z = solve(Q_op)
#Sigma.Z =  Sigma.Z
processes_list = list(noise = "NIG", nu = 1, mu = 1)
processes_list$V <- list()

vvT <- 0
FSigma <- 0
dSigma_dsd <- 0
dSigma_dbeta <- 0
Fsbf = 0
Fsbr = 0 
Fs   = 0
Fbeta =  0
Ftau  = 0
for(i in 1:n.pers)
{
  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],rnorm(n.obs[i]))
  
  #Y[[i]] <- rnorm(n = n.obs[i], B_fixed[[i]]%*%betaf, sd = sd_Y) + B_random[[i]]%*%mvrnorm(n = 1,mu = betar, Sigma = Sigma)
  processes_list$V[[i]] <- operator_list$h
  A = spde.A(locs[[i]],operator_list$loc[[1]])
  Sigma_Y =  B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + sd_Y^2*diag(n.obs[i])
  Sigma_Y = Sigma_Y +  A%*%Sigma.Z%*%t(A) 
  Y[[i]] <- mvrnorm(n = 1, mu =  B_fixed[[i]]%*%betaf + B_random[[i]]%*%betar  , Sigma = Sigma_Y)
  v = Y[[i]] - B_fixed[[i]]%*%betaf - B_random[[i]]%*%betar  
  vvT <- v%*%t(v)
  Q = solve(Sigma_Y)
  #Sigma.Z
  D<-getDuplicateM(2)
  ddSigma <- 0.5*kronecker(t(B_random[[i]]),t(B_random[[i]]))%*%kronecker(Q,Q - 2*Q%*%vvT%*%Q)%*%kronecker(B_random[[i]],B_random[[i]])
  FSigma <- -t(D)%*%ddSigma%*%D + FSigma
  dSigma_dsd <- dSigma_dsd + sd_Y *  as.vector(Q%*%Q%*%vvT%*%Q + Q%*%vvT%*%Q%*%Q)%*%kronecker(B_random[[i]],B_random[[i]])%*%D
  dSigma_dsd <- dSigma_dsd - sd_Y * as.vector(Q%*%Q)%*%kronecker(B_random[[i]],B_random[[i]])%*%D
  #ddsigma_eps
  dl = (-sd_Y*sum(diag(Q)) + sd_Y*t(v)%*%Q%*%Q%*%v)
  Fs = Fs - (dl/sd_Y + 2*sd_Y^2*sum(diag(Q%*%Q)) - 4*(sd_Y^2)*t(v)%*%Q%*%Q%*%Q%*%v)
 }


mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = betar,
                          beta_fixed  = betaf,
                          Sigma = Sigma,
                          noise = "Normal")
mixedEffect_list2  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = betar,
                          beta_fixed  = betaf,
                          Sigma = Sigma,
                          noise = "NIG",
                          mu = as.matrix(c(0,1)),
                          nu = 2
                          )
sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = list(sigma = sd_Y, noise = "Normal"),
                              processes_list    = processes_list,
                              operator_list     = operator_list)

processes_list$X <- sim_res$X
processes_list$V <- sim_res$V

#measurment_list  = list(sigma = sd_Y, noise = "Normal")
measurment_list  = list(sigma = sd_Y, noise = "NIG", nu = 2, Vs = list(rep(1,n.obs[1])))
res <- estimateLong(Y                = Y,
                    locs             = locs,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list  = measurment_list,
                    processes_list   = processes_list,
                    operator_list    = operator_list,
                    nSim             = nSim,
                    nBurnin_learningrate = 0,
                    alpha = 0.3,
                    pSubsample = pSubsample,
                    step0 = 0.3,
                    nIter = nIter,
                    silent = silent,
                    polyak_rate = -1,
                    seed = seed,
                    nBurnin = nBurnin,
                    estimate_fisher = 2)

res$FisherMatrix-t(res$FisherMatrix)
