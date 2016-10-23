###
# example for testing estimation of random effect.
# using Normal Mixed effect NIG measurement error
###
graphics.off()
library(testthat)
library(LDMod)
library(rGIG)
set.seed(1)
nu <- 3
n.pers <- 10
n.obs  <- 200

sd_beta <- 0.01
sd_Y    <- 0.1

B_list <- list()
beta <- c(0.9,0.4)
beta_list <- list()
Y_list <- list()
Vin        <- list()
for(i in 1:n.pers)
{
  B_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  beta_list[[i]] <-beta+  rnorm(n = length( beta), 0, sd = sd_beta)
  V <- rGIG(rep(-0.5,n.obs), rep(nu, n.obs), rep(nu, n.obs))
  Y_list[[i]]        <-  B_list[[i]]%*%beta_list[[i]] + sqrt(V)*rnorm(n = n.obs, 0, sd = sd_Y)
  Vin[[i]] <- rep(1, n.obs)
}
meas_list <- list(Vs = Vin, sigma.eps = 1., nu = 3., noise = "NIG")

mixedEffect_list <- list(B_random = B_list, 
                         Sigma = sd_beta*diag(2), 
                         beta_random = c(0.,0.),
                         noise = "Normal")
input <- list(Y = Y_list, 
              mixedEffect_list = mixedEffect_list,
              measurementError_list = meas_list,
              nSim = 2,
              alpha = 0.3,
              step0 = 1,
              nIter = 2000,
              silent = 0)
res <- estimateME(Y = Y_list, 
                  mixedEffect_list = mixedEffect_list,
                  measurment_list = meas_list,
                  nSim = 2,
                  alpha = 0.3,
                  step0 = 1,
                  nIter = 2000,
                  silent = 0)


test_that("simple Gaussian-IG random effect",
{
  expect_equal(c(res$mixedEffect_list$beta_random), beta, tolerance  = 0.1)
})
test_that("simple Gaussian-IG measurement sigma",
{
  expect_equal(res$measurementError_list$sigma, sd_Y, tolerance  = 0.1)
})
if(1){
x11()
par(mfrow=c(3,1))
plot(res$measurementError_list$sigma_vec, type='l', col='red')
n_ <- length(res$measurementError_list$sigma_vec)
lines(c(1, n_), c(sd_Y[1],sd_Y[1]))
#lines(c(1, n_), c(sd_Y[2],sd_Y[2]))
plot(res$measurementError_list$nu_vec, type='l', col='red')
lines(c(1, n_), c(nu[1],nu[1]))

plot(res$mixedEffect_list$betar_vec[,2], type='l', col='red')
lines(c(1, n_), c(beta[2],beta[2]))
}