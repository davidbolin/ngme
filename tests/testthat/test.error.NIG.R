##
#  test mixed effect p measurement error
#  D:2019-01-15
##
library(testthat)
library(ngme)
set.seed(4)

n_iter <- 1000

nindv <- 20
n     <- 200

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
nu           <- 0.6
mu           <- 1.2


data <-c()
for(indv in 1:nindv){
  B_fixed  <- cbind(rep(1, n), rnorm(n))
  B_random <- as.matrix(1:n)
  V      = rGIG(rep(-0.5,n),
                rep( nu, n),
                rep( nu, n),
                as.integer(1000 * runif(1) ))
  E      = sigma_random * sqrt(V) * rnorm(n)
  Y         <- B_fixed%*%beta_fixed +
               B_random%*%(beta_random  + sigma*rnorm(1)) + E
  Ya <- Y + (V-1)*mu
  id <- rep(indv, n)
  data <- rbind(data, cbind(B_fixed,
                            B_random,
                            V,
                            Y,
                            Ya,
                            id))
}
dimnames(data)[[2]] <- c('B1','B2','B3','V','Y','Ya','id')


NIGMVD_ass <- ngme( fixed       = Ya ~ B1 + B2,
                random      = ~ -1+B3|id,
                data        = as.data.frame(data),
                error       = 'NIG',
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
if(0){
fish.est      <- ngme.fisher(NIGMVD_ass,
                                            nSim = 20,
                                            nIter = 200,
                                            nBurnin=20,
                                            n.cores = 1,
                                            observed = TRUE,
                                            silent = TRUE)}


context("error")
test_that("assymetric NIG error ", {
theta_est <- c(NIGMVD_ass$measurementError_list$nu,
               NIGMVD_ass$measurementError_list$sigma,
               NIGMVD_ass$mixedEffect_list$beta_random,
               NIGMVD_ass$mixedEffect_list$beta_fixed,
               NIGMVD_ass$mixedEffect_list$Sigma,
               NIGMVD_ass$measurementError_list$mu)
theta <- c(nu, sigma, beta_random, beta_fixed, sigma_random^2,mu)
expect_equal(theta,
             theta_est,
             tolerance = 0.1)
})

