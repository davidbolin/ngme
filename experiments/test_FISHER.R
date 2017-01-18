#TODO: for processes, create idenity operator, then compute
#      the actual integral
rm(list=ls())
library(LDMod)
library(MASS)
seed     <- 9
silent   <- 1
plotflag <- 1

nIter <- 1
n.pers <- 20
nSim  <- 100
n.obs  <- 10 + 0*(1:n.pers)

nBurnin = 10
pSubsample = 1
nPar_burnin = 0
n = 10
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()

betaf <- as.matrix(c(2.1))
betar = as.matrix(c(0.8,0.5))
tau <- 1.1
kappa <- 2
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


operator_list <- create_operator(locs, n, name = "Matern")
operator_list$tau   <- tau
operator_list$kappa   <- kappa
K = 0.5  * operator_list$tau*(kappa^(-1.5) * operator_list$G[[1]] + kappa^(0.5) * operator_list$C[[1]])
Q_op = t(K)%*%diag(1/operator_list$h[[1]])%*%K
#operator_list <- create_operator(locs, n, name = "fd2")
#operator_list$tau   <- tau
#operator_list$Q[[1]] <- K/operator_list$tau 
Sigma.Z = solve(Q_op)
# debugging differntial matrices
A = spde.A(locs[[1]],operator_list$loc[[1]])
eps <- 10^-6
K_eps        = 0.5  * operator_list$tau*((kappa + eps)^(-1.5) * operator_list$G[[1]] + (kappa + eps)^(0.5) * operator_list$C[[1]])
K_meps       = 0.5  * operator_list$tau*((kappa - eps)^(-1.5) * operator_list$G[[1]] + (kappa - eps)^(0.5) * operator_list$C[[1]])
K_kappa_tau  = 0.5  * (operator_list$tau+eps)*((kappa + eps)^(-1.5) * operator_list$G[[1]] + (kappa + eps)^(0.5) * operator_list$C[[1]])
K_mkappa_tau = 0.5  * (operator_list$tau+eps)*((kappa - eps)^(-1.5) * operator_list$G[[1]] + (kappa - eps)^(0.5) * operator_list$C[[1]])
K_mkappa_mtau        = 0.5  * (operator_list$tau-eps)*((kappa - eps)^(-1.5) * operator_list$G[[1]] + (kappa - eps)^(0.5) * operator_list$C[[1]])
K_kappa_mtau        = 0.5  * (operator_list$tau-eps)*((kappa + eps)^(-1.5) * operator_list$G[[1]] + (kappa + eps)^(0.5) * operator_list$C[[1]])

Sigma.Z_eps = solve(t(K_eps)%*%diag(1/operator_list$h[[1]])%*%K_eps)
Sigma.Z_kappa_tau = solve(t(K_kappa_tau)%*%diag(1/operator_list$h[[1]])%*%K_kappa_tau)
Sigma.Z_meps = solve(t(K_meps)%*%diag(1/operator_list$h[[1]])%*%K_meps)
Sigma.Z_mkappa_tau = solve(t(K_mkappa_tau)%*%diag(1/operator_list$h[[1]])%*%K_mkappa_tau)
Sigma.Z_mkappa_mtau = solve(t(K_mkappa_mtau)%*%diag(1/operator_list$h[[1]])%*%K_mkappa_mtau)
Sigma.Z_kappa_mtau = solve(t(K_kappa_mtau)%*%diag(1/operator_list$h[[1]])%*%K_kappa_mtau)

#Sigma.Z =  Sigma.Z
processes_list = list(noise = "Normal")
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
Fkappa <- 0
Fkappa_tau <-  0
for(i in 1:n.pers)
{
  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],rnorm(n.obs[i]))
  
  #Y[[i]] <- rnorm(n = n.obs[i], B_fixed[[i]]%*%betaf, sd = sd_Y) + B_random[[i]]%*%mvrnorm(n = 1,mu = betar, Sigma = Sigma)
  processes_list$V[[i]] <- operator_list$h
  A = spde.A(locs[[i]],operator_list$loc[[1]])
  Sigma_Y =  B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + sd_Y^2*diag(n.obs[i])
  Sigma_Y_eps_kappa = Sigma_Y +  A%*%Sigma.Z_eps%*%t(A) 
  Sigma_Y_meps_kappa = Sigma_Y +  A%*%Sigma.Z_meps%*%t(A) 
  Sigma_Y_kappa_tau = Sigma_Y +  A%*%Sigma.Z_kappa_tau%*%t(A) 
  Sigma_Y_mkappa_tau       = Sigma_Y +  A%*%Sigma.Z_mkappa_tau%*%t(A) 
  Sigma_Y_mkappa_mtau       = Sigma_Y +  A%*%Sigma.Z_mkappa_mtau%*%t(A) 
  Sigma_Y_kappa_mtau       = Sigma_Y +  A%*%Sigma.Z_kappa_mtau%*%t(A) 
  
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
      #ddtau
      Q_A =  A%*%Sigma.Z%*%t(A) 
      Q_tilde = Q%*%Q_A%*%Q
      Z <-  v%*%t(v)
      Ftau = Ftau + 0.5 * (4/tau^2) * sum(diag( Q_A%*%Q_tilde ))
      Ftau = Ftau - 0.5 * (8/tau^2) * sum(diag( Sigma_Y%*%Q_tilde%*%Q_A%*%Q )) 
      Ftau = Ftau + 0.5 * (6/tau^2) * sum(diag( Sigma_Y%*%Q_tilde ))
      Ftau = Ftau + 0.5 * (8/tau^2) * sum(diag( Q_tilde%*%Q_A%*%Q%*%Z )) 
      Ftau = Ftau - 0.5 * (6/tau^2) * sum(diag( Q_tilde%*%Z ))
      
      
     Fsbf = Fsbf + (2*sd_Y)*t(v)%*%Q%*%Q%*%B_fixed[[i]]
      Fsbr = Fsbr + (2*sd_Y)*t(v)%*%Q%*%Q%*%B_random[[i]]
      B = cbind(B_fixed[[i]], B_random[[i]])
      Fbeta <- Fbeta + as.matrix(t(B)%*%Q%*%B)
  dSigma_dbeta <- dSigma_dbeta + kronecker(t(cbind(B_fixed[[i]],B_random[[i]]))%*%Q,t(v)%*%Q)%*%kronecker(B_random[[i]],B_random[[i]])%*%D

  #ddkappa 
  #
  
  lik = - 0.5 * log(det(Sigma_Y)) - sum(diag(solve(Sigma_Y,Z)))/2
  lik_eps = - 0.5 * log(det(Sigma_Y_eps_kappa)) - sum(diag(solve(Sigma_Y_eps_kappa,Z)))/2
  lik_meps = - 0.5 * log(det(Sigma_Y_meps_kappa)) - sum(diag(solve(Sigma_Y_meps_kappa,Z)))/2
  ddkappa =  (-2*lik + lik_eps + lik_meps)/eps^2
  Fkappa = Fkappa - ddkappa
  #dkappa dtau
  lik_tau_kappa = - 0.5 * log(det(Sigma_Y_kappa_tau)) - sum(diag(solve(Sigma_Y_kappa_tau,Z)))/2
  lik_mkappa_tau       = - 0.5 * log(det(Sigma_Y_mkappa_tau)) - sum(diag(solve(Sigma_Y_mkappa_tau,Z)))/2
  lik_mkappa_mtau       = - 0.5 * log(det(Sigma_Y_mkappa_mtau)) - sum(diag(solve(Sigma_Y_mkappa_mtau,Z)))/2
  lik_kappa_mtau       = - 0.5 * log(det(Sigma_Y_kappa_mtau)) - sum(diag(solve(Sigma_Y_kappa_mtau,Z)))/2
  dkappadtau = (lik_tau_kappa - lik_kappa_mtau - lik_mkappa_tau + lik_mkappa_mtau)/(4*eps^2)
  Fkappa_tau = Fkappa_tau - dkappadtau 
}


mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = betar,
                          beta_fixed  = betaf,
                          Sigma = Sigma,
                          noise = "Normal")

sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = list(sigma = sd_Y, noise = "Normal"),
                              processes_list    = processes_list,
                              operator_list     = operator_list)
processes_list$X <- sim_res$X


res <- estimateLong(Y                = Y,
                    locs             = locs,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list  = list(sigma = sd_Y, noise = "Normal"),
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

#expect_equal(max(abs(c(F_fixed[1],diag(F_random))/diag(res$FisherMatrix[1:3,1:3])-1)),0,tolerance=0.02)

Fish <- 0 * res$FisherMatrix
Fish[1:3,1:3  ]     <- Fbeta 
Fish[1:3, 4:6   ]   <- as.matrix(dSigma_dbeta)
Fish[4:6  , 1:3 ]   <- as.matrix(dSigma_dbeta)
Fish[4:6  , 4:6   ] <- as.matrix(FSigma)
Fish[7  , 7   ]     <-as.matrix(Fs)
Fish[1:3  , 7  ]    <- as.matrix(cbind(Fsbf,Fsbr))
Fish[7  ,   1:3 ]    <-t(as.matrix(cbind(Fsbf,Fsbr)))
Fish[7  , 4:6   ] <-dSigma_dsd
Fish[4:6  , 7   ] <-dSigma_dsd
Fish[8,8 ]        <- Ftau
Fish[8,9]         <- Fkappa_tau
Fish[9,8]         <- Fkappa_tau
Fish[9,9]         <- Fkappa
#print(res$FisherMatrix[1:7,1:7])
#print(Fish)
print(res$FisherMatrix[8:9,8:9]/Fish[8:9,8:9])