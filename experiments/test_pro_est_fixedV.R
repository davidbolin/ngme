#simulate V
#simulate X
#simulate Y
#use second order RW
rm(list=ls())
library(testthat)
library(LDMod)
library(methods)
graphics.off()

plot_flag <- TRUE
seed <- 4
set.seed(seed)
noises <- c("CH")
for(k in 1:length(noises)){
  pSubsample <- 0.5
  nBurnin <- 200
nobs  <- 50
nIter <- 1000 #100 is good enough
n     <- 100 #n grid points
learning_rate <- 0.99
nu_true <- 10
mu_true <- 0
nu_guess <- 20
mu_guess <- 20
tau_geuss <- 0.5
theta <- list()
theta$sigma <- 0.1 # meas error
theta$tau   <- 0.5
theta$nu    <- nu_true
theta$mu    <- mu_true
locs   <- list()
for(i in 1:nobs)
{ 
  locs[[i]]   <- seq(0, 1, length = nobs)
}

output_sim <- simulateLong.R(locs, 
               theta,
               noise = noises[k],
               operatorType = "fd2",
               n = n)
operator_list <- create_operator(locs, n, name = "fd2")
obs_list <- list()
X        <- list()
V        <- list()
for(i in 1:length(locs)){
  obs_list[[i]] <- list(A =  spde.A(x = operator_list$loc[[1]], loc = locs[[i]]), 
                        Y=output_sim$Y[[i]], 
                        locs = locs[[i]])
  X[[i]] <- rep(0, n) 
  V[[i]] <- operator_list$h
}


mError_list <- list(noise = "Normal",
                    sigma = theta$sigma)

mixedEffect_list  <- list(noise="Normal")
operator_list$tau <- tau_geuss
processes_list <- list(nu = nu_guess, 
                       mu = mu_guess, 
                       X = output_sim$X, 
                       V = output_sim$V, 
                       noise = noises[k])
input <- list( obs_list         = obs_list,
               operator_list    = operator_list,
               processes_list   = processes_list,
               nBurnin_base     = 3,
               nIter            = nIter,     # iterations to run the stochastic gradient
               nSim             = 4,
               nBurnin          = nBurnin,   # steps before starting gradient estimation
               silent           = 0, # print iteration info)
               step0            = 0.3,
               alpha            = 0.1,
               learning_rate    = learning_rate,
               pSubsample       = pSubsample,
               polyak_rate      = -1,
               subsample_type = 1,
               measurementError_list   = mError_list,
               mixedEffect_list = mixedEffect_list,
               sampleX = 1,
               sampleV = 0,
               seed   = seed
              )
output <- estimateLong_cpp(input)

if(plot_flag){
x11()
par(mfrow=c(3,2))
plot(locs[[1]],output_sim$Y[[5]])
lines(output$operator_list$loc[[1]], output_sim$X[[5]])
lines(output$operator_list$loc[[1]], output$Xs[[5]],col='red',lty='dashed')
n_ <- length(output$operator_list$tauVec)
if(noises[k] != "CH"){
  
  plot(output$processes_list$mu_vec)
  lines(c(1, n_), c( mu_true, mu_true), col='red' )
  plot(output$processes_list$nu_vec)
  lines(c(1, n_), c( nu_true, nu_true), col='red' )
  
}
plot(output$operator_list$tauVec)
lines(c(1, n_), c( theta$tau, theta$tau), col='red' )

}
if(noises[k] != "CH"){
  
  test_that(paste("nu with known X,V, noise = ",noises[k],sep=""),{
    expect_equal( theta$nu, output$processes_list$nu, tolerance  = 0.2)
  })
  test_that(paste("mu with known X,V, noise = ",noises[k],sep=""),{
    expect_equal( theta$mu, output$processes_list$mu, tolerance  = 0.2)
  } )
}
test_that(paste("tau with known X,V, noise = ",noises[k],sep=""),{
  expect_equal( theta$tau, output$operator_list$tau, tolerance  = 0.2)
})

}

#Eigen::VectorXd temp(Vs[i].size());
#temp.array() = Vs[i].array().log();
#dnu  +=  h_sum[i] * (1. + log(nu)) + h[i].dot(temp) - Vs[i].sum() - h_digamma[i];
#ddnu += h_sum[i]/ nu - h_trigamma[i];