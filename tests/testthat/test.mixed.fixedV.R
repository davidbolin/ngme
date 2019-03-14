##
#  test mixed effect p measurement error
#  D:2019-01-15
##
library(testthat)
library(ngme)
graphics.off()
set.seed(130)

n_iter <- 4000

nindv <- 200
n     <- 110

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
mu_mixed     <- 1
nu_mixed     <- 1#0.5


B_fixed  <- list()
B_random <- list()
B_sigma  <- list()
Y        <- list()
V_mixed  <- rep(nindv)
for(indv in 1:nindv){
  B_fixed[[indv]]  <- cbind(rep(1, n), rnorm(n))
  B_random[[indv]] <- as.matrix(1:n)
  V_mixed[indv] <-rGIG(rep(-0.5,1),
                 rep( nu_mixed, 1),
                 rep( nu_mixed, 1),
                 as.integer(1000 * runif(1) ))
  V_   <- V_mixed[indv]
  Y[[indv]]        <- B_fixed[[indv]]%*%beta_fixed +
                   B_random[[indv]]%*%(beta_random + mu_mixed*(-1+V_) + sqrt(V_)*sigma_random*rnorm(1)) +
                   sigma * rnorm(n)
  
  B_sigma[[indv]]  <- as.matrix(rep(1, n))
  
}

context("mixed with fixed V")
mixed_list <- list(B_fixed    = B_fixed, 
                   B_random   = B_random,
                   V          = V_mixed,
                   fixedV     = 1,
                   name       = 'NIG')
error_list <- list(name       = 'Normal',
                   B          = B_sigma)
res2 <- test_mixed(n_iter, 
                   Y, 
                   mixed_list,
                   error_list)
mixed_list <- list(B_fixed    = B_fixed, 
                   B_random   = B_random,
                   V          = V_mixed,
                   fixedV     = 0,
                   name       = 'NIG')
res3 <- test_mixed(n_iter, 
                   Y, 
                   mixed_list,
                   error_list)
print(res2$mixedEffect_list$mu)
print(res3$mixedEffect_list$mu)
x11()
par(mfrow=c(2,1))
plot(sqrt(res2$mixedEffect_list$Sigma_vec))
plot(sqrt(res3$mixedEffect_list$Sigma_vec))
cat('sigma (V known)= ',sqrt(res2$mixedEffect_list$Sigma),'\n')
cat('sigma = ',sqrt(res3$mixedEffect_list$Sigma),'\n')
x11()
x <-seq(-4,10,length.out = 2000)
dens <-dnig(x, -as.vector(res3$mixedEffect_list$mu), as.vector(res3$mixedEffect_list$mu), res3$mixedEffect_list$nu,  sqrt(res3$mixedEffect_list$Sigma[1,1]))
dens2 <-dnig(x, -mu_mixed,mu_mixed, nu_mixed, sigma_random)
plot(x,dens2,type='l')
lines(x, dens, col='red')
x11()
par(mfrow=c(3,1))
hist(res2$mixedEffect_list$U,20)
hist(res3$mixedEffect_list$U,20)
hist(res3$mixedEffect_list$U-res2$mixedEffect_list$U)


