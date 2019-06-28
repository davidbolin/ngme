##
#  test mixed effect p measurement error
#  D:2019-03-15
##
rm(list=ls())
graphics.off()
library(testthat)
library(ngme)
set.seed(4)

n_iter <- 500

nindv <-  200
n     <-  50

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
nu           <- 0.6

data <-c()
for(indv in 1:nindv){
  B_fixed  <- cbind(rep(1, n), rnorm(n))
  B_random <- as.matrix(1:n)
  V      = 1
  E      = sigma_random * sqrt(V) * rnorm(n)
  Y         <- B_fixed%*%beta_fixed +
               B_random%*%(beta_random  + sigma*rnorm(1)) + E
  Ya <- Y 
  id <- rep(indv, n)
  data <- rbind(data, cbind(B_fixed,
                            B_random,
                            V,
                            Y,
                            Ya,
                            id))
}
dimnames(data)[[2]] <- c('B1','B2','B3','V','Y','Ya','id')
data <- as.data.frame(data)
beta_f <- c(1.13732, 3.06337)
beta_r <- 0.959745
sigma_random = 0.448103
sigma        = sqrt( 0.00212802)
grad_sum <- rep(0,3)
for(i in 1:2){
  data1 <- data[data$id==i,]
  res  <- data1$Y - data1$B1*beta_f[1] - data1$B2*beta_f[2] - data1$B3*beta_r
  
  Qhat <- (1/sigma_random^2)* t(data1$B3)%*%data1$B3  + 1/sigma^2
  mu_hat <- solve(Qhat, t(data1$B3)%*%res/sigma_random^2)
  #print(mu_hat)
  #cat('res = ',res,'\n')
  #cat('(res - Ajoint * mu_hat) = ',(res- data1$B3%*%mu_hat),'\n')
  grad <- t(data1[,c("B1","B2","B3")])%*%(res- data1$B3%*%mu_hat)/sigma_random^2
  grad_sum <- grad_sum + grad
}
GAUSMIX <- ngme( fixed       = Ya ~ -1 + B1 + B2,
                random      = ~ -1+B3|id,
                data        = as.data.frame(data),
                error       = 'Normal',
                error_assymetric=T,
                use.process = F,
                silent      = T,
                nIter       = n_iter,
                controls    = list(estimate.fisher = FALSE,
                                   subsample.type  = 0,
                                   nSim  =1,
                                   nBurnin = 2,
                                   alpha = 0.1,
                                   step0 = 0.9))



test_that("Gaussian mixed effect", {
theta_est <- c(
               GAUSMIX$measurementError_list$sigma,
               GAUSMIX$mixedEffect_list$beta_fixed,
               GAUSMIX$mixedEffect_list$beta_random,
               GAUSMIX$mixedEffect_list$Sigma)
theta <- c( sigma_random, beta_fixed,beta_random, sigma)
expect_equal(theta,
             theta_est,
             tolerance = 0.1)
})