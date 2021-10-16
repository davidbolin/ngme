##
#  test mixed effect p measurement error
##
rm(list=ls())
graphics.off()
library(ngme)
set.seed(8)
sim <-  1
n_iter <- 15000

nindv <- 500 
n     <- 30

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
nu_mixed     <- 17
mu_mixed     <- 2
trace_random <- trace_random2 <- matrix(0,nrow = sim, ncol = n_iter)
est_param <- est_param2 <- matrix(0, nrow=sim, ncol = 3)
sd_param  <- sd_param2 <- sd_param_mc <- matrix(0, nrow=sim, ncol = 3)
true_param <- c(beta_fixed, beta_random)

data <- c()
for(indv in 1:nindv){
  B_fixed  <- cbind(as.matrix(1:n)/n, rnorm(n))
  B_random <- as.matrix(rep(1, n))
  V_mixed <-rGIG(rep(-0.5,1),
                 rep( nu_mixed, 1),
                 rep( nu_mixed, 1),
                 as.integer(1000 * runif(1) ))
  E      = sigma_random * rnorm(n)
  Y         <- B_fixed%*%beta_fixed + B_random%*%(beta_random  + (-1+V_mixed)*mu_mixed+ sqrt(V_mixed)*sigma*rnorm(1)) + E
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
data <- as.data.frame(data)
if(nindv >= 100){
  control = list(estimate.fisher = FALSE,
                 pSubsample = 0.2,
                 alpha = 0.2,
                 subsample.type  = 1,
                 nSim  =2,
                 nBurnin = 2)
}else{
  control = list(estimate.fisher = FALSE,
                 pSubsample = 0.2,
                 subsample.type  = 0,
                 nSim  =2,
                 alpha = 0.2,
                 nBurnin = 2)
}

res <- ngme(fixed       = Ya ~ B1 + B2,
                random      = ~ -1+B3|id,
                data        = data,
                reffects    = 'NIG',
                use.process = FALSE,
                silent      = FALSE,
                controls.init = list(nIter.init=500),
                nIter  = 1000,
                controls    = control)

res <- ngme.par(n.cores = 5, std.lim = Inf,max.rep = 10,
                fixed       = Ya ~ B1 + B2,
                random      = ~ -1+B3|id,
                data        = data,
                reffects    = 'NIG',
                use.process = FALSE,
                silent      = FALSE,
                controls.init = list(nIter.init=500),
                nIter  = 1000,
                controls    = control)

fiher_NIG <- ngme.fisher(res,
                         nSim = 20,
                         nIter = 10,
                         nBurnin=1,
                         n.cores = 1,
                         n.rep = 1,
                         std.threshold = 2,
                         observed = TRUE,
                         only.effects=F,
                         silent = TRUE)
V1 <- sqrt(diag(solve(fiher_NIG$fisher_est))[3] + res$fixed_est_var[10,1])
V2 <- sqrt(diag(solve(fiher_NIG$fisher_est))[3]) 
cat(round(fiher_NIG$mixedEffect_list$beta_random -1.644854 * V1,2),
    round(fiher_NIG$mixedEffect_list$beta_random +1.644854 * V1,2),'\n')