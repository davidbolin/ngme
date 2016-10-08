##
# Test for general GAL processes
#
##

rm(list=ls())
#library(testthat)
library(LDMod)
#library(rGIG)
library(methods)
graphics.off()

plot_flag <- TRUE

nobs  <- 10
nIter <- 400
n     <- 200 #n grid points
operatorType <- "matern"
nu_true <- 30
mu_true <- 1
nu_guess <- 10
mu_guess <- 20
tau_geuss <- 0.5
theta <- list()
theta$sigma <- 0.1 # meas error
theta$tau   <- 0.5
theta$nu    <- nu_true
theta$mu    <- mu_true
theta$kappa <- 1
locs   <- list()
for(i in 1:nobs)
{ 
  locs[[i]]   <- seq(0, 1, length = nobs)
}

output_sim <- simulateLong.R(locs, 
                             theta,
                             noise = "GAL",
                             operatorType =operatorType,
                             n = n)
operator_list <- create_operator(locs, n, name = operatorType)

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
operator_list$kappa <- theta$kappa
processes_list <- list(nu = nu_guess, 
                       mu = mu_guess, 
                       X = output_sim$X, 
                       V = output_sim$V, 
                       noise = "GAL")
input <- list( obs_list         = obs_list,
               operator_list    = operator_list,
               processes_list   = processes_list,
               nIter            = nIter,     # iterations to run the stochastic gradient
               nSim             = 1,
               nBurnin          = 10,   # steps before starting gradient estimation
               silent           = 0, # print iteration info)
               step0            = 1,
               alpha            = 0.01,
               pSubsample       = 1,
               subsample_type   = 1,
               measurementError_list   = mError_list,
               mixedEffect_list = mixedEffect_list,
               sampleX = 1,
               sampleV = 1
)
output <- estimateLong_cpp(input)
