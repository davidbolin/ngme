##
#  test mixed effect p measurement error
#  D:2019-01-15
##
library(testthat)
library(ngme)
set.seed(4)

n_iter <- 100000

nindv <- 100
n     <- 110

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5


B_fixed  <- list()
B_random <- list()
B_sigma  <- list()
Y        <- list()
for(indv in 1:nindv){
  B_fixed[[indv]]  <- cbind(rep(1, n), rnorm(n))
  B_random[[indv]] <- as.matrix(1:n)
  Y[[indv]]        <- B_fixed[[indv]]%*%beta_fixed +
                   B_random[[indv]]%*%(beta_random  + sigma_random*rnorm(1)) +
                   sigma * rnorm(n)
  B_sigma[[indv]]  <- as.matrix(rep(1, n))
  
}

mixed_list <- list(B_fixed    = B_fixed, 
                   B_random   = B_random,
                   beta_random = beta_random,
                   beta_fixed  = beta_fixed,
                   Sigma       = as.matrix(sigma_random),
                   name       = 'Normal')
error_list <- list(name       = 'Normal',
                   B          = B_sigma,
                   sigma      = sigma)
res <- test_mixed_Fisher(n_iter, 
                  Y, 
                  mixed_list,
                  error_list)
res_est <- test_mixed(n_iter, 
                         Y, 
                         mixed_list,
                         error_list)
Fish <- res$d2Given-res$VVt + outer(res$V,res$V)
print(Fish)