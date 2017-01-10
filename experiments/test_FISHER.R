library(LDMod)
library(MASS)
seed     <- 2
silent   <- 1
plotflag <- 1

nIter <- 1
n.pers <- 2
nSim  <- 1000
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
Ff <- 0
Fr <- matrix(0,2,2)
Frf <- matrix(0,2,1)
set.seed(seed)
FS <- diag(c(0,0,0))
FSs <- FSbf <- c(0,0,0)
FSbr <- matrix(0,nrow = 3, ncol = 2)
Fs <- 0
Fsbf <- 0
Fsbr <- 0
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
                  estimate_fisher = 2)
print(F)
print(res$FisherMatrix)
