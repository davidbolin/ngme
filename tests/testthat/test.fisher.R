require(testthat)
context("Fisher")

test_that("Fisher, Gaussian fixed effects", {

  graphics.off()
  library(LDMod)
  library(MASS)
  seed     <- 5
  silent   <- 1
  plotflag <- 1

  nIter <- 50
  pSubsample <- 1
  nSim <- 2
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
                    estimate_fisher = TRUE)

  #print(F)
  #print(res$FisherMatrix)

  expect_equal(max(abs(F/res$FisherMatrix-1)),0,tolerance=0.01)
})


test_that("Fisher, Gaussian random effects", {
  library(LDMod)
  library(MASS)
  seed     <- 5
  silent   <- 1
  plotflag <- 1

  nIter <- 10
  n.pers <- 2
  nSim  <- 5
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
  sd_Y = 2
  Ff <- 0
  Fr <- matrix(0,2,2)
  Frf <- matrix(0,2,1)
  set.seed(seed)
  FS <- diag(c(0,0,0))
  FSs <- FSbf <- c(0,0,0)
  FSbr <- matrix(0,nrow = 3, ncol = 2)
  for(i in 1:n.pers)
  {
    B_fixed[[i]]  <- as.matrix(rep(1, n.obs[i]))
    B_random[[i]] <- cbind((1:n.obs[i])/n.obs[i],((1:n.obs[i])/n.obs[i])^(2))

    Y[[i]] <- rnorm(n = n.obs[i], B_fixed[[i]]%*%betaf, sd = sd_Y) + B_random[[i]]%*%mvrnorm(n = 1,mu = betar, Sigma = Sigma)
    Q =  solve(B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + sd_Y^2*diag(n.obs[i]))
    Ff  <- Ff  + t(B_fixed[[i]] )%*%Q%*%B_fixed[[i]]
    Fr  <- Fr  + t(B_random[[i]])%*%Q%*%B_random[[i]]
    Frf <- Frf + t(B_random[[i]])%*%Q%*%B_fixed[[i]]

    v = Y[[i]] - B_fixed[[i]]%*%betaf - B_random[[i]]%*%betar
    dl = (-sd_Y*sum(diag(Q)) + sd_Y*t(v)%*%Q%*%Q%*%v)
    Fs = Fs - (dl/sd_Y + 2*sd_Y^2*sum(diag(Q%*%Q)) - 4*sd_Y^2**t(v)%*%Q%*%Q%*%Q%*%v)
    Fsbf = Fsbf + (2*sd_Y)*t(v)%*%Q%*%Q%*%B_fixed[[i]]
    Fsbr = Fsbr + (2*sd_Y)*t(v)%*%Q%*%Q%*%B_random[[i]]

    vvv = matrix(c(1,0,0,0,1,0,0,0,0,0,0,1),3,4)
    for(k in 1:3){
      dSigma  <- matrix(vvv[k,],2,2)
      dQS = B_random[[i]]%*%dSigma%*%t(B_random[[i]]) + sd_Y^2*diag(n.obs[i])
      dQs = B_random[[i]]%*%Sigma%*%t(B_random[[i]]) + 2*sd_Y*diag(n.obs[i])
      FS[k,k] = FS[k,k] - (0.5*sum(diag(-Q%*%dQS%*%Q%*%dQS)) - t(v)%*%Q%*%dQS%*%Q%*%dQS%*%Q%*%v)
      FSs[k] = FSs[k] - (0.5*sum(diag(-Q%*%dQs%*%Q%*%dQS)) - t(v)%*%Q%*%dQs%*%Q%*%dQS%*%Q%*%v - t(v)%*%Q%*%dQS%*%Q%*%dQs%*%Q%*%v)
      FSbf[k] = FSbf[k] + 0.5*t(v)%*%Q%*%dQS%*%Q%*%B_fixed[[i]]
      FSbr[k,] = FSbr[k,] + 0.5*t(v)%*%Q%*%dQS%*%Q%*%B_random[[i]]
    }

  }

  F <- rBind(c(Ff,t(Frf),FSbf,Fsbf),
             cBind(Frf,Fr,t(FSbr),t(Fsbr)),
             cBind(FSbf,FSbr,FS,FSs),
             c(Fsbf,Fsbr,FSs,Fs))


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
                    estimate_fisher = TRUE)
  print(F)
  print(res$FisherMatrix)

  expect_equal(max(abs(F/res$FisherMatrix-1)),0,tolerance=0.01)
})


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
