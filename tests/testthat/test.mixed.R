##
#  test mixed effect vs lmm
#  D:2019-05-31
##
rm(list=ls())
graphics.off()
library(testthat)
library(ngme)
library(lme4)
set.seed(4) #4

n_iter <- 50

nindv <- 50
n     <- 30

beta_random  <- c(0.8, 0.6)
beta_fixed   <- c(1.1, 2.2)
Sigma        <- diag(c(2,1))
sigma_random <- 0.5
data <-c()
for(indv in 1:nindv){
  B_fixed  <- cbind(as.matrix(1:n)/n, rnorm(n))
  B_random <- cbind(rep(1, n), rnorm(n))
  U = mvrnorm(n=1, mu = c(0,0),Sigma = Sigma)
  E      = sigma_random * rnorm(n)
  Y         <- B_fixed%*%beta_fixed +
    B_random%*%(beta_random  + U) + E
  Ya <- Y 
  id <- rep(indv, n)
  data <- rbind(data, cbind(B_fixed,
                            B_random,
                            Y,
                            Ya,
                            id))
}
dimnames(data)[[2]] <- c('B1','B2','B3','B4','Y','Ya','id')


res.ngme <- ngme( fixed       = Ya ~ -1 + B1 + B2,
                    random      = ~ -1 + B3+B4|id,
                    data        = as.data.frame(data),
                    reffects    = 'Normal',
                    use.process = F,
                    silent      = T,
                    controls.init = list(nIter.init=100),
                    nIter  = n_iter,
                    controls    = list(estimate.fisher = FALSE,
                                       subsample.type  = 0,
                                       pSubsample = 1,
                                       nSim  =2,
                                       nBurnin = 2,
                                       alpha = 0.3,
                                       step0 = 1))




fisher.ngme <- ngme.fisher(res.ngme,
                         nSim = 2,
                         nIter = 1,
                         nBurnin = 1,
                         n.cores = 1,
                         n.rep = 1,
                         std.threshold = 2,
                         observed = TRUE,
                         only.effects=T,
                         silent = T)

lmm <- lmer(Y ~-1 + B1 + B2 +B3 + B4+ (-1 +  B3 + B4| id), data = as.data.frame(data),
            REML = FALSE)
slmm <- summary(lmm)
beta_lmm <- as.vector(slmm$coefficients[,1])
beta_est<- c(res.ngme$mixedEffect_list$beta_fixed,
             res.ngme$mixedEffect_list$beta_random)
###
# debug plots
###

test_that("Normal mixed effect", {
  
  expect_equal(beta_lmm,
               beta_est,
               tolerance = 10^-4)
  expect_equal(slmm$sigma,
               res.ngme$measurementError_list$sigma,
               tol=10^-2)
  expect_equal( c(as.matrix(slmm$vcov)),
                c(as.matrix(solve(fisher.ngme$fisher_est))),
               tol=10^-3)
  expect_equal( c(res.ngme$mixedEffect_list$Sigma),
               c(as.matrix(as.data.frame(slmm$varcor$id))),
                tol = 10^-3)
})