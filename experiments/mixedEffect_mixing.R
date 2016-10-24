####
# testing the mixing properties of the mixed effect for case
# when Sigma close to singular!
#
#####
rm(list = ls())
graphics.off()
sim <- 50000
burnin = 10000
library(rGIG)
library(mvtnorm)
nu.mixed <- 0.9910863
mu.mixed <- c(-0.0906397 , -0.00919812 )
beta_random <- c(4.739068 , -0.03971646 )
Sigma <- matrix(c(0.04879072 , -0.0008268011 , -0.0008268011 , 1.413677e-05 ), nrow= 2 , ncol= 2 )
sigma_eps =  0.1442434
source("Brandom.data")
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
  b_U <- (t(B_random)%*%Y)/sigma_eps^2
  Q_U <- (t(B_random)%*%B_random)/sigma_eps^2
  b   = b_U + Q%*%(-mu.mixed + V*mu.mixed + beta_random)
  return(rmvnorm(n = 1, mean = solve(Q_U + Q, b), sigma = solve(Q + Q_U)))
}
sampleV <- function(U)
{
  
  #sample V
  U_ <- U + mu.mixed - beta_random
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




for(i in 1:sim){
  #sample U
  U1[i, ] <- sampleU(V1[i])
  #sample V
  V1[i+1] <- sampleV(U1[i,])
}
V1 <- V1[1:sim]
x11()
par(mfrow=c(3,1))
acf(V1[burnin:sim],200)
plot(V1)
plot(U1[,1])


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

if(0){
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