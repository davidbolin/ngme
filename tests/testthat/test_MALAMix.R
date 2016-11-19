####
# Comparing the Gibbs sampler vs MALA for NIG mixed effect
#
#####
rm(list = ls())
library(rGIG)
library(mvtnorm)
graphics.off()

sim <- 5000
burnin = 200
seed <- 2
n_B <- 4 # number of coeff
n_Y <- 20
Sigma <-  10^-4*diag(n_B)
sigma_eps =  0.01442434
mu.mixed <- seq(0,1, length = n_B)
nu.mixed <- 0.0910863


B_random <- cbind(rep(1, n_Y),
           seq(0, 1, length=n_Y),
           matrix(rnorm((n_B-2) * n_Y),nrow= n_Y, ncol=n_B-2))


set.seed(seed)

beta_random <- rep(0, n_B)

# sampling, prior
V <- rGIG(-0.5, nu.mixed, nu.mixed)
mu <-   - mu.mixed +  mu.mixed * V
U <- t(rmvnorm(n = 1, mean = mu, sigma = sqrt(V)*Sigma))
Y <- B_random%*%(beta_random + U) + sigma_eps * rnorm(n_Y,1)