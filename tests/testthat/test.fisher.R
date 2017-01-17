require(testthat)
context("Fisher")
library(LDMod)
library(MASS)

test_that("Fisher, Gaussian fixed effects", {

  graphics.off()
  seed     <- 5
  silent   <- 1
  plotflag <- 1

  nIter <- 50
  pSubsample <- 1
  nSim <- 20
  n.pers <- 10 #number of patients
  n.obs  <- 10 #number of obs per patient

  sd_Y    <- 0.5 # error of the noise

  betaf <- c(1.,-1)
  nu <- 1
  betar_list <- list()
  Bf_list    <- list()
  V_list     <- list()
  Y_list     <- list()
  set.seed(seed)
  B_sum <- 0
  d2s = 0
  dsb = 0
  for(i in 1:n.pers)
  {
    Bf_list[[i]]    <- cbind(1:n.obs,rep(1, n.obs))
    Y_list[[i]]        <- rnorm(n = n.obs, Bf_list[[i]]%*%betaf, sd = sd_Y)
    B_sum <- B_sum + t(Bf_list[[i]])%*%Bf_list[[i]]/sd_Y^2
    #Sigma.hat = sd_Y^2*diag(n.obs)
    #Q.hat = solve(Sigma.hat)
    v = Y_list[[i]] - Bf_list[[i]]%*%betaf
    d2s = d2s - (n.obs/sd_Y^2 - (3/sd_Y^4)*t(v)%*%v)
    dsb = dsb + (2/sd_Y^3)*t(v)%*%Bf_list[[i]]
  }

  F = rBind(cBind(B_sum,t(dsb)),c(dsb,d2s))

  meas_list <- list(sigma = sd_Y, noise = "Normal")
  mixedEffect_list <- list(B_fixed  = Bf_list,
                           beta_fixed  = betaf,
                           noise = "Normal")



  res <- estimateME(Y = Y_list,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list = meas_list,
                    nSim = nSim,
                    alpha = 0.3,
                    pSubsample = pSubsample,
                    step0 = 0.3,
                    nIter = nIter,
                    silent = silent,
                    polyak_rate = -1,
                    seed = seed,
                    estimate_fisher = 2)

  #print(F)
  #print(res$FisherMatrix)

  expect_equal(max(abs(F/res$FisherMatrix-1)),0,tolerance=0.01)
})


test_that("Fisher, Gaussian random effects", {
  library(LDMod)
  library(MASS)
  seed     <- 4
  silent   <- 1
  plotflag <- 1
  
  nIter <- 1
  n.pers <- 2
  nSim  <- 20
  n.obs  <- 10 + 0*(1:n.pers)
  
  nBurnin = 40
  pSubsample = 1
  nPar_burnin = 100
  Y <- list()
  locs <- list()
  B_random <- list()
  B_fixed  <- list()
  Vin <- list()
  
  betaf <- as.matrix(c(2.1))
  betar = as.matrix(c(0.8,0.5))
  
  Sigma <- matrix(c(0.02,0,0,0.01),2,2)
  sd_Y = 0.1
  set.seed(seed)
  vvT <- 0
  FSigma <- 0
  dSigma_dsd <- 0
  dSigma_dbeta <- 0
  Fsbf = 0
  Fsbr = 0 
  Fs   = 0
  for(i in 1:n.pers)
  {
    B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
    B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],rnorm(n.obs[i]))
    
    Y[[i]] <- rnorm(n = n.obs[i], B_fixed[[i]]%*%betaf, sd = sd_Y) + B_random[[i]]%*%mvrnorm(n = 1,mu = betar, Sigma = Sigma)
    Q =  solve(B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + sd_Y^2*diag(n.obs[i]))
    v = Y[[i]] - B_fixed[[i]]%*%betaf - B_random[[i]]%*%betar
    vvT <- v%*%t(v)
    
    D<-getDuplicateM(2)
    ddSigma <- 0.5*kronecker(t(B_random[[i]]),t(B_random[[i]]))%*%kronecker(Q,Q - 2*Q%*%vvT%*%Q)%*%kronecker(B_random[[i]],B_random[[i]])
    FSigma <- -t(D)%*%ddSigma%*%D + FSigma
    dSigma_dsd <- dSigma_dsd + sd_Y *  as.vector(Q%*%Q%*%vvT%*%Q + Q%*%vvT%*%Q%*%Q)%*%kronecker(B_random[[i]],B_random[[i]])%*%D
    dSigma_dsd <- dSigma_dsd - sd_Y * as.vector(Q%*%Q)%*%kronecker(B_random[[i]],B_random[[i]])%*%D
     dl = (-sd_Y*sum(diag(Q)) + sd_Y*t(v)%*%Q%*%Q%*%v)
     Fs = Fs - (dl/sd_Y + 2*sd_Y^2*sum(diag(Q%*%Q)) - 4*(sd_Y^2)*t(v)%*%Q%*%Q%*%Q%*%v)
     Fsbf = Fsbf + (2*sd_Y)*t(v)%*%Q%*%Q%*%B_fixed[[i]]
     Fsbr = Fsbr + (2*sd_Y)*t(v)%*%Q%*%Q%*%B_random[[i]]
    dSigma_dbeta <- dSigma_dbeta + kronecker(t(cbind(B_fixed[[i]],B_random[[i]]))%*%Q,t(v)%*%Q)%*%kronecker(B_random[[i]],B_random[[i]])%*%D
  }
  
  
  
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = betar,
                            beta_fixed  = betaf,
                            Sigma = Sigma,
                            noise = "Normal")
  
  
  res <- estimateME(Y = Y,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list = list(sigma = sd_Y, noise = "Normal"),
                    nSim = nSim,
                    alpha = 0.3,
                    pSubsample = pSubsample,
                    step0 = 0.3,
                    nIter = nIter,
                    silent = silent,
                    polyak_rate = -1,
                    seed = seed,
                    nBurnin = nBurnin,
                    estimate_fisher = 2)
  expect_equal(mean(abs(FSigma - res$FisherMatrix[4:6,4:6])), 0, tolerance  = 100)
  expect_equal(mean(abs(dSigma_dsd-res$FisherMatrix[4:6,7])), 0, tolerance  = 100)
  expect_equal(mean(abs(res$FisherMatrix[1:3,4:6]-dSigma_dbeta)), 0, tolerance  = 100)
  
})

if(0){
test_that("Fisher, Gaussian process", {
library(LDMod)
library(MASS)
rm(list=ls())
seed     <- 5
silent   <- 1
plotflag <- 1

nIter <- 4000
n.pers <- 2
nSim  <- 25
n.obs  <- 10 + 0*(1:n.pers)

nBurnin = 40
pSubsample = 1
nPar_burnin = 100
n = 10
Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()

betaf <- as.matrix(c(2.1))
betar = as.matrix(c(0.8,0.5))

Sigma <- matrix(c(0.02,0,0,0.01),2,2)
sd_Y = 2
F_fixed <- 0
F_random <- 0
set.seed(seed)


for(i in 1:n.pers)
{
  B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
  B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))
  locs[[i]] <- 1:n.obs[i]
  Y[[i]] <- rep(1,n.obs[i])
}

operator_list <- create_operator(locs, n, name = "fd2")
operator_list$tau   <- 75
K = operator_list$tau*operator_list$Q[[1]]
Sigma.Z = solve(t(K)%*%diag(1/operator_list$h[[1]])%*%K)

processes_list = list(noise = "Normal")
processes_list$V <- list()
for(i in 1:length(locs))
{
  processes_list$V[[i]] <- operator_list$h
  A = spde.A(locs[[i]],operator_list$loc[[1]])
  Q =  B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + A%*%Sigma.Z%*%t(A) + sd_Y^2*diag(n.obs[i])
  F_fixed <-  F_fixed  + t(B_fixed[[i]] )%*%solve(Q,B_fixed[[i]])
  F_random <- F_random + t(B_random[[i]])%*%solve(Q,B_random[[i]])

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

res <- estimateLong(Y                = sim_res$Y,
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
                    estimate_fisher = TRUE)

expect_equal(max(abs(c(F_fixed[1],diag(F_random))/diag(res$FisherMatrix[1:3,1:3])-1)),0,tolerance=0.02)

})
}

lNIG <- function(U, sigma, nu)
{
  p = - 1
  b = U^2/sigma^2 + nu
  a =  nu
  logf =log(a)  - 0.5 * log(b)
  sqrt_ab = sqrt(a * b)
  K1 <- besselK(sqrt_ab, p, expon.scaled=T)
  
  logf = logf + log(K1) - sqrt_ab
  return(logf)
}


EiV_NIG <- function(U, sigma, nu)
{
  p <- -1
  b <- U^2/sigma^2 + nu
  a <- nu
  sqrt_ab = sqrt(a * b)
  K1 <- besselK(sqrt_ab, p, expon.scaled=T)
  K0 <- besselK(sqrt_ab, p+1, expon.scaled=T)
  
  sqrt_a_div_b <- sqrt(a/b)
  EiV = K0 / K1
  EiV = EiV * sqrt_a_div_b - (2 * p) * 1/b
  return(EiV)
}
EV_NIG <- function(U, sigma, nu)
{
  p <- -1
  b <- U^2/sigma^2 + nu
  a <- nu
  sqrt_ab = sqrt(a * b)
  K1 <- besselK(sqrt_ab, p, expon.scaled=T)
  K0 <- besselK(sqrt_ab, p+1, expon.scaled=T)
  
  sqrt_b_div_a <- sqrt(b/a)
  EV = K0 / K1
  EV = EV * sqrt_b_div_a
  return(EV)
}
dEiV_NIG <- function(U, sigma, nu)
{
  p <- -1
  b <- U^2/sigma^2 + nu
  a <- nu
  sqrt_ab = sqrt(a * b)
  K1 <- besselK(sqrt_ab, p, expon.scaled=T)
  K0 <- besselK(sqrt_ab, p+1, expon.scaled=T)
  
  sqrt_a_div_b <- sqrt(a/b)
  EiV = K0 / K1
  EiV = EiV * sqrt_a_div_b - (2 * p) * 1/b
  
  K0dK1 = K0 / K1
  dEiV = 0
  dEiV = -1 - (p+1) * K0dK1 / sqrt_ab
  dEiV = dEiV - (-K0^2 + (p/sqrt_ab) *K1 * K0)/K1^2
  dEiV = dEiV * 0.5 * sqrt_a_div_b
  dEiV = dEiV  * sqrt_a_div_b
  dEiV = dEiV - 0.5 * K0dK1 * sqrt_a_div_b / b
  dEiV = dEiV +  (2 * p) / b^2
  dEiV = c(dEiV * 2) * (U / sigma^2)
  return(dEiV)
}

ddU_NIG <- function(U, sigma, nu)
{
  # 
  EiV = EiV_NIG(U, sigma, nu)
  ddU = EiV /sigma^2
  dEiV <- dEiV_NIG(U, sigma, nu)
  d_ <- dEiV*U/sigma^2
  ddU  = ddU + d_
  return(ddU)
}


test_that("Fisher, NIG noise", {
silent   <- 1

nIter <- 1
pSubsample <- 1
nSim <- 200
n.pers <- 10 #number of patients
n.obs  <- 400 #number of obs per patient

sd_Y    <- 0.5 # error of the noise

betaf <- c(1.1,-2.1)
nu <- 1
betar_list <- list()
Bf_list    <- list()
V_list     <- list()
Y_list     <- list()
set.seed(seed)
B_sum <- 0
d2s = 0
dBeta <- 0
ddBeta <- 0
dsb = 0
resid <- c()
dBetadSigma <- 0
dBetadnu <- 0
ddsigma <- 0
ddnu <- 0
dnuds <- 0
for(i in 1:n.pers)
{
  Bf_list[[i]]    <- cbind(1:n.obs,rep(1, n.obs))
  V_list[[i]] <- LDMod::rGIG(rep(-0.5,n.obs), rep(nu, n.obs), rep(nu,n.obs), sample.int(10^6,1))
  Y_list[[i]]        <- Bf_list[[i]]%*%betaf  + sqrt(V_list[[i]]) * rnorm(n = n.obs, 0, sd = sd_Y)
  v = Y_list[[i]] - Bf_list[[i]]%*%betaf
  resid <- c(resid,v)
  eps <- 10^-6
  ddNIG_U_num <- (lNIG(v - eps,sd_Y, nu) + lNIG(v + eps,sd_Y, nu) -  2* lNIG(v,sd_Y, nu))/eps^2
  ddNIG <- ddU_NIG(v, sd_Y, nu) 
  ddBeta <-  ddBeta - t(Bf_list[[i]]  )%*%diag(as.vector(ddNIG_U_num), nrow=n.obs)%*%Bf_list[[i]] 
  
  EiV = EiV_NIG(v, sd_Y,  nu)
  EV = EV_NIG(v, sd_Y, nu)
  ddsigma <- ddsigma + 3 * t(v)%*% (EiV*v) / sd_Y^4 - n.obs/sd_Y^2
  dEiV = (EiV_NIG(v, sd_Y + eps ,  nu) - EiV )/eps
  EiV_eps <- EiV_NIG(v, sd_Y  ,  nu +  eps)
  dEiV_nu = (EiV_eps - EiV )/eps
  ddsigma <- ddsigma - t(v)%*% (dEiV*v) / sd_Y^3
  
  dnuds <- dnuds - t(v)%*% (dEiV_nu*v) / sd_Y^3
  EV_eps <- EV_NIG(v, sd_Y  ,  nu +  eps)
  
  dnu <-  (0.5 * ( n.obs / nu + 2 * n.obs -   (sum(EiV) + sum(EV)) ))
  dnu_eps <-  (0.5 * ( n.obs / (nu+eps) + 2 * n.obs -   (sum(EiV_eps) + sum(EV_eps)) ))
  ddnu <- ddnu - (dnu_eps - dnu) / eps
  dBeta = dBeta - (EiV * v) /sd_Y^2
  dBetadnu = dBetadnu -  t(Bf_list[[i]]  )%*%(dEiV_nu * v) /sd_Y^2
  dBetadSigma = dBetadSigma + 2 *t(Bf_list[[i]]  )%*% (EiV * v) /sd_Y^3  -t(Bf_list[[i]]  )%*%(dEiV * v) /sd_Y^2
}

dU_num =  (lNIG(v + eps,sd_Y, nu) -   lNIG(v,sd_Y, nu))/eps

meas_list <- list(sigma = sd_Y,  nu = nu, Vs = V_list, noise = "NIG")
mixedEffect_list <- list(B_fixed  = Bf_list,
                         beta_fixed  = betaf,
                         noise = "Normal")
#x_<-seq(min(resid),max(resid),length=1000)
#x11()
#par(mfrow=c(2,1))
#plot(density(resid))
#plot(x_,exp(lNIG(x_,sd_Y,nu)),type='l')


res <- estimateME(Y = Y_list,
                  mixedEffect_list = mixedEffect_list,
                  measurment_list = meas_list,
                  nSim = nSim,
                  alpha = 0.3,
                  pSubsample = pSubsample,
                  step0 = 0.3,
                  nIter = nIter,
                  silent = silent,
                  polyak_rate = -1,
                  seed = seed,
                  estimate_fisher = 2)

Fish <- 0 * res$FisherMatrix
Fish[1:2,1:2  ] <- ddBeta 
Fish[1:2, 3   ] <- dBetadSigma
Fish[3  , 1:2 ] <- dBetadSigma
Fish[1:2, 4   ] <- dBetadnu
Fish[4  , 1:2 ] <- dBetadnu
Fish[3  , 3   ] <-ddsigma
Fish[4  , 4   ] <-ddnu
Fish[3  , 4   ] <-dnuds
Fish[4  , 3   ] <-dnuds
expect_equal((max(abs(solve(Fish)-solve(res$FisherMatrix)))),0,tol=0.05)
})
