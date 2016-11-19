##
# write a Gibbs sampler source file, use it store simple data
# explore the posterior!
##
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
graphics.off()
seed <- 2
set.seed(seed)
sim <- 20000
burnin = 2000
library(rGIG)
library(mvtnorm)
nu.mixed <- 100.9910863
mu.mixed <- c(-0.2 , 0.2 )
beta_random <- c(4.739068 , -0.03971646 )
Sigma <-  10^-4*matrix(c(0.2,0.15,0.15,0.2), ncol = 2, nrow = 2)
sigma_eps =  0.01442434
source("Brandom.data")
source("mixedEffect_MALA.R")
n <- dim(B_random)[1]
# sampling, prior
V <- LDMod::rGIG(-1, nu.mixed, nu.mixed, seed)
mu <- - mu.mixed +  mu.mixed * V
U_true <- t(rmvnorm(n = 1, mean = mu, sigma = sqrt(V)*Sigma))
Y <- B_random%*%(beta_random + U_true) + sigma_eps * rnorm(n,1)


source("Gibbs_mixed.test.R")
loglik_v <- vector(mode = "numeric", length= sim)
for(i in 1:sim){
  #sample U
  U1[i, ] <- sampleU(V1[i])
  #sample V
  V1[i+1] <- sampleV(U1[i,])
  loglik_v[i] <-lNIG(U1[i, ], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma_inv, mu.mixed, nu.mixed)
}

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
data <- data.frame(U1 = U1[burnin:sim,1], U2 = U1[burnin:sim,2])
p <- ggplot(data, aes( U1, U2))
p <- p + stat_bin2d() + stat_bin2d(bins=20) + scale_fill_gradientn(colours=r)
x11()
plot(U1[,1])
x11()
print(p)

U_out =  MALA(U1[sim, ], Y - B_random%*%beta_random, B_random, sigma_eps, Sigma, mu.mixed, nu.mixed)


  U <- U1[sim,]
  sigma <- sigma_eps
  Sigma <- Sigma
  mu    <- mu.mixed
  nu <- nu.mixed
  B <- B_random
  res <- Y - B_random%*%beta_random
  sigma2_ <- 1
  
  #compusted gradient and second derivative and expectatin of second derivative?
  dU_list <- dU_NIG(U, res, B, sigma, Sigma, mu, nu)
  dU_list$ddU <- dU_list$ddU/sigma2_
  
  
  # sample the proposal
  R <- chol(dU_list$ddU)
  U_mean <- U + solve(dU_list$ddU, t(dU_list$dU))/2
  Ustar = U_mean + solve(R, rnorm(length(U)))
  Ustar = c(Ustar)
  
  # computes for the new location
  dU_star_list <- dU_NIG(Ustar, res, B, sigma, Sigma, mu, nu)
  dU_star_list$ddU <- dU_star_list$ddU/sigma2_
  U_star_mean <- Ustar + solve(dU_star_list$ddU, t(dU_star_list$dU))/2
  Rstar <- chol(dU_star_list$ddU)
  
  # the proposal distributions
  qstar <- sum( log( diag( R))) - 0.5* t(Ustar - U_mean)%*%dU_list$ddU%*%(Ustar - U_mean)
  q     <-sum( log( diag( Rstar)))  - 0.5* t(U - U_star_mean)%*%dU_star_list$ddU%*%(U - U_star_mean)
  
  # the loglikelihoods
  loglik <- lNIG(U, res, B, sigma, solve(Sigma), mu, nu)
  loglikstar <- lNIG(Ustar, res, B, sigma, solve(Sigma), mu, nu)
  
  d_num <- (loglik - lNIG(U + c(10^-6,0), res, B, sigma, solve(Sigma), mu, nu))/10^-6
  d_num <- c(d_num,(loglik - lNIG(U + c(0,10^-6), res, B, sigma, solve(Sigma), mu, nu))/10^-6)

  
  df <- expand.grid(x = seq(min(U1[burnin:sim,1]), max(U1[burnin:sim,1]), length=20),
                    y =  seq(min(U1[burnin:sim,2]), max(U1[burnin:sim,2]), length=20))
  
  df$z <- vector(mode = "numeric", length = length(df$x))
  for(i in 1:length(df$x)) 
    df$z[i] <- lNIG(c(df$x[i], df$y[i]), res, B, sigma, solve(Sigma), mu, nu)
  # default is compatible with geom_tile()
  df$z <- exp(df$z - max(df$z))
  x11()
  print(ggplot(df, aes(x, y, fill = z)) + geom_raster())
  VarU <- var(U1[burnin:sim,])
