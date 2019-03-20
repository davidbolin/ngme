##
#  test mixed effect p measurement error
#  D:2019-01-15
##
rm(list=ls())
graphics.off()
library(testthat)
library(ngme)
set.seed(7)

n_iter <- 200

nindv <- 200
n     <- 30

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
nu_mixed     <- 1
mu_mixed     <- 2
data <-c()
for(indv in 1:nindv){
  B_fixed  <- cbind(as.matrix(1:n), rnorm(n))
  B_random <- as.matrix(rep(1, n))
  V_mixed <-rGIG(rep(-0.5,1),
                 rep( nu_mixed, 1),
                 rep( nu_mixed, 1),
                 as.integer(1000 * runif(1) ))
  E      = sigma_random * rnorm(n)
  Y         <- B_fixed%*%beta_fixed +
               B_random%*%(beta_random  + (-1+V_mixed)*mu_mixed+ sqrt(V_mixed)*sigma*rnorm(1)) + E
  Ya <- Y 
  id <- rep(indv, n)
  data <- rbind(data, cbind(B_fixed,
                            B_random,
                            V_mixed,
                            Y,
                            Ya,
                            id))
}
dimnames(data)[[2]] <- c('B1','B2','B3','V','Y','Ya','id')


NIGMVD_ass <- ngme( fixed       = Ya ~ B1 + B2,
                random      = ~ -1+B3|id,
                data        = as.data.frame(data),
                reffects    = 'NIG',
                use.process = F,
                silent      = T,
                controls.init = list(nIter.init=100),
                nIter  = n_iter,
                controls    = list(estimate.fisher = FALSE,
                                   subsample.type  = 0,
                                   nSim  =2,
                                   nBurnin = 2,
                                   alpha = 0.1,
                                   step0 = 0.5))




theta_est <- c(NIGMVD_ass$measurementError_list$sigma,
               NIGMVD_ass$mixedEffect_list$nu,
               NIGMVD_ass$mixedEffect_list$mu,
               NIGMVD_ass$mixedEffect_list$beta_random,
               NIGMVD_ass$mixedEffect_list$beta_fixed,
               sqrt(NIGMVD_ass$mixedEffect_list$Sigma))
theta <- c(sigma_random,nu_mixed, mu_mixed, beta_random,beta_fixed, sigma)



###
# debug plots
###
if(0){
  x11()
  par(mfrow=c(2,4))
  plot(NIGMVD_ass$mixedEffect_list$betaf_vec[,1],type='l',ylab='betaf 1')
  plot(NIGMVD_ass$mixedEffect_list$betaf_vec[,2],type='l',ylab='betaf 2')
  plot(NIGMVD_ass$mixedEffect_list$betar_vec,type='l',ylab='betar')
  plot(sqrt(NIGMVD_ass$mixedEffect_list$Sigma_vec),type='l',ylab='mixed sigma')
  plot(NIGMVD_ass$mixedEffect_list$nu_vec,type='l',ylab='mixed nu')
  plot(NIGMVD_ass$mixedEffect_list$mu_vec,type='l',ylab='mixed mu')
}
test_that("NIG mixed effect", {
expect_equal(theta,
             theta_est,
             tolerance = 0.1)
})