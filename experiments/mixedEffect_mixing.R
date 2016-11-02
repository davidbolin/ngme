####
# testing the mixing properties of the mixed effect for case
# when Sigma close to singular!
#
#####
#rm(list = ls())
graphics.off()
sim <- 500
burnin = 100
library(rGIG)
library(mvtnorm)
nu.mixed <- 0.9910863
mu.mixed <- c(-0.2 , 0.2 )
beta_random <- c(4.739068 , -0.03971646 )
Sigma <-  10^-2*matrix(c(0.2,0.15,0.15,0.2), ncol = 2, nrow = 2)
sigma_eps =  0.1442434
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
V_0 <- rGIG(-0.5, nu.mixed, nu.mixed) 
V1    <- matrix(0, nrow = sim + 1, ncol = 1)
V1[1] <- V_0
U1    <- matrix(0, nrow= sim, ncol=  2)
p_GIG = 0.5 * (-1 - dim(B_random)[2])
eSigma <- eigen(Sigma) 
Sigma_inv <- eSigma$vectors%*%diag(1/eSigma$values)%*%t(eSigma$vectors)
R         <- eSigma$vectors%*%diag(sqrt(1/eSigma$values))%*%t(eSigma$vectors)
R_sigma   <- eSigma$vectors%*%diag(sqrt(eSigma$values))%*%t(eSigma$vectors)

b_U <- (t(B_random)%*%(Y - B_random%*%beta_random))/sigma_eps^2
Q_U <- (t(B_random)%*%B_random)/sigma_eps^2

sampleU <- function(V)
{
  Q   =  Sigma_inv/V
  #b_U <- (t(B_random)%*%Y)/sigma_eps^2
  #Q_U <- (t(B_random)%*%B_random)/sigma_eps^2
  b   = b_U + Q%*%(-mu.mixed + V*mu.mixed)
  return(rmvnorm(n = 1, mean = solve(Q_U + Q, b), sigma = solve(Q + Q_U)))
}

sampleV <- function(U)
{
  
  #sample V
  U_ <- U + mu.mixed #- beta_random
  b = t(U_)%*%Sigma_inv%*%(U_ )
  b = b + nu.mixed
  a_GIG = t(mu.mixed)%*%Sigma_inv%*%mu.mixed + nu.mixed
  return(rGIG(p_GIG, a_GIG, b))  
}
sampleV2 <- function(E)
{
  
  #sample V
  E_ <- E + R%*%(mu.mixed - beta_random)
  b = t(E_)%*%E_ 
  b = b + nu.mixed
  a_GIG = t(mu.mixed)%*%Sigma_inv%*%mu.mixed + nu.mixed
  return(rGIG(p_GIG, a_GIG, b))  
}

sampleU2 <- function(V)
{
  Q   =  diag(2)/V
  B_ <-B_random%*%R_sigma
  b_U <- (t(B_)%*%Y)/sigma_eps^2
  Q_U <- (t(B_)%*%B_)/sigma_eps^2
  b   = b_U + Q%*%(R%*%(  beta_random -mu.mixed + V*mu.mixed))
  return(rmvnorm(n = 1, mean = solve(Q_U + Q, b), sigma = solve(Q + Q_U)))
}



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
loglik[1] <-  loglik[sim]
for(i in 2:sim)
{
  
 U1[i,] =  MALA(U1[i-1, ], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma, mu.mixed, nu.mixed)
 loglik[i] <-lNIG(U1[i, ], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma_inv, mu.mixed, nu.mixed)
}
x11()
par(mfrow=c(3,1))
acf(U1[burnin:sim,1],200)
plot(U1[,1])
plot(loglik)
print(mean(abs(diff(U1[,1]))>0))
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
if(0){
for(i in 1:sim){
  #sample U
  U1[i, ] <- sampleU2(V1[i])
  #sample V
  V1[i+1] <- sampleV2(U1[i,])
}
x11()
par(mfrow=c(3,1))
acf(V1[burnin:sim],200)
plot(V1)
plot(U1[,1])

##
# sampling posterior version 2
# sampling (V | X-mu*V, X - mu * V |V)
##
V2    <- matrix(0, nrow = sim + 1, ncol = 1)
V2[1] <- V_0
U2    <- matrix(0, nrow= sim, ncol=  2)
p_GIG = 0.5 * (-1 - dim(B_random)[2])
Sigma_inv <- solve(Sigma)
a_GIG =  nu.mixed
b_U <- (t(B_random)%*%(Y - B_random%*%beta_random))/sigma_eps^2
Q_U <- (t(B_random)%*%B_random)/sigma_eps^2
for(i in 1:sim){
  #sample U
  U2[i, ] <- sampleU(V2[i])
  #sample V
  U_ <- U2[i, ] + mu.mixed - V2[i]*mu.mixed
  b = t(U_)%*%Sigma_inv%*%U_ + nu.mixed
  V2[i+1] <- rGIG(p_GIG, nu.mixed, b)
  print(V2[i+1])
  if(V2[i+1]>1000)
    break
}
V2 <- V2[1:sim]
x11()
par(mfrow=c(2,1))
acf(V2[burnin:sim],200)
plot(V2)
}

