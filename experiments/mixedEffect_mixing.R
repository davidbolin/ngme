####
# testing the mixing properties of the mixed effect for case
# when Sigma close to singular!
#
#####
rm(list = ls())
graphics.off()
seed <- 2
set.seed(seed)
sim <- 5000
burnin = 200
library(rGIG)
library(mvtnorm)
nu.mixed <- 0.0910863
mu.mixed <- c(-0.1 , 0.1 )
beta_random <- c(4.739068 , -0.03971646 )
Sigma <-  10^-4*matrix(c(0.2,0.15,0.15,0.2), ncol = 2, nrow = 2)
sigma_eps =  0.01442434
source("Brandom.data")
source("mixedEffect_MALA.R")
n <- dim(B_random)[1]

# sampling, prior
V <- rGIG(-0.5, nu.mixed, nu.mixed)
mu <-   - mu.mixed +  mu.mixed * V
U <- t(rmvnorm(n = 1, mean = mu, sigma = sqrt(V)*Sigma))
Y <- B_random%*%(beta_random + U) + sigma_eps * rnorm(n,1)
##
# sampling posterior version 1
# sampling (X, V)
##


source("Gibbs_mixed.test.R")

loglik <- vector(mode = "numeric", length= sim)
for(i in 1:sim){
  #sample U
  U1[i, ] <- sampleU(V1[i])
  #sample V
  V1[i+1] <- sampleV(U1[i,])
  loglik[i] <-lNIG(U1[i, ], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma_inv, mu.mixed, nu.mixed)
}

V1 <- V1[1:sim]
x11()
par(mfrow=c(4,1))
acf(V1[burnin:sim],200)
plot(V1)
plot(U1[,1])
plot(loglik)
U1[1,] <- U1[sim,] 
U2 <- U1
loglik[1] <-  loglik[sim]
for(i in 2:sim)
{
  
  U2[i,] =  MALA(U2[i-1, ], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma, mu.mixed, nu.mixed,
                 scale = 1)
 loglik[i] <-lNIG(U2[i, ], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma_inv, mu.mixed, nu.mixed)
}
x11()
par(mfrow=c(3,1))
acf(U2[burnin:sim,1],200)
plot(U2[,1])
plot(loglik)
print(mean(abs(diff(U2[,1]))>0))
du_lik <- dU_NIG(U1[1,], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma, mu.mixed, nu.mixed)
loglik_1 <-lNIG(U1[1, ], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma_inv, mu.mixed, nu.mixed)
loglik_1e <-lNIG(U1[1, ] + c(0,10^-6), Y - B_random%*%beta_random, B_random, sigma_eps, Sigma_inv, mu.mixed, nu.mixed)
loglik_2e <-lNIG(U1[1, ] + c(10^-6,0), Y - B_random%*%beta_random, B_random, sigma_eps, Sigma_inv, mu.mixed, nu.mixed)
du_ <- c(loglik_1, loglik_1)/10^-6 - c(loglik_2e, loglik_1e)/10^-6
#dEiV_NIG
dEiV = dEiV_NIG(U, Sigma, mu.mixed, nu.mixed)
EiV = EiV_NIG(U, Sigma, mu.mixed, nu.mixed)
EiV_e = EiV_NIG(U+ c(0,10^-6), Sigma, mu.mixed, nu.mixed)
dEiV_num = (EiV - EiV_e)/10^-6
