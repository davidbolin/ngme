##
#  test mixed effect p measurement error
##
rm(list=ls())
graphics.off()
library(testthat)
library(ngme)
set.seed(3)
sim <-  200       
n_iter <- 5000

nindv <- 1200
n     <- 30

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
nu_mixed     <- 1
mu_mixed     <- 2
trace_random <- trace_random2 <- matrix(0,nrow = sim, ncol = n_iter)
est_param <- est_param2 <- matrix(0, nrow=sim, ncol = 3)
sd_param  <- sd_param2 <- matrix(0, nrow=sim, ncol = 3)
true_param <- c(beta_fixed, beta_random)
sum_res <- function(est_param, se_param, true_param){
  lower <- matrix(0, nrow=dim(est_param)[1], ncol = dim(est_param)[2])
  upper <- lower
  in_CI <- lower
  for(i in 1:dim(est_param)[1]){
  lower[i,]   <- est_param[i,] -  1.644854*se_param[i,]
  upper[i,]   <- est_param[i,] +  1.644854*se_param[i,]
  in_CI[i,]   <-  lower[i,] <= true_param & true_param <= upper[i,]
  }
  return(list(est   = est_param,
              se    = se_param,
              in_CI = in_CI,
              lower = lower,
              upper = upper ))
}


for(i in 1:sim){
  cat('iter = ',i,'of ',sim,'\n')
  data <-c()
  for(indv in 1:nindv){
    B_fixed  <- cbind(as.matrix(1:n)/n, rnorm(n))
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
                                     pSubsample = 0.4,
                                     subsample.type  = 0,
                                     nSim  =2,
                                     nBurnin = 2,
                                     alpha = 0.3,
                                     step0 = 0.99))
  
  fiher_NIG <- ngme.fisher(NIGMVD_ass,
                           nSim = 20,
                           nIter = 10,
                           nBurnin=1,
                           n.cores = 1,
                           n.rep = 1,
                           std.threshold = 2,
                           observed = TRUE,
                           only.effects=F,
                           silent = TRUE)
  est_param[i,] <- c(NIGMVD_ass$mixedEffect_list$beta_fixed, NIGMVD_ass$mixedEffect_list$beta_random)
  sd_param[i,]  <-  sqrt(diag(solve(fiher_NIG$fisher_est)))[1:length(est_param[i,])]
  trace_random[i,]<-NIGMVD_ass$mixedEffect_list$betar_vec
  if(0){
  NIGMVD_ass <- ngme( fixed       = Ya ~ B1 + B2,
                      random      = ~ -1+B3|id,
                      data        = as.data.frame(data),
                      reffects    = 'NIG',
                      use.process = F,
                      silent      = T,
                      nIter  = n_iter,
                      controls    = list(estimate.fisher = FALSE,
                                         pSubsample = 0.4,
                                         subsample.type  = 0,
                                         nSim  =2,
                                         nBurnin = 2,
                                         alpha = 0.3,
                                         step0 = 0.99),
                      init.fit = NIGMVD_ass)
  fiher_NIG <- ngme.fisher(NIGMVD_ass,
                            nSim = 20,
                            nIter = 10,
                            nBurnin=1,
                            n.cores = 1,
                            n.rep = 1,
                            std.threshold = 2,
                            observed = TRUE,
                            only.effects=F,
                            silent = TRUE)
  
  est_param2[i,] <- c(NIGMVD_ass$mixedEffect_list$beta_fixed, NIGMVD_ass$mixedEffect_list$beta_random)
  sd_param2[i,]  <-  sqrt(diag(solve(fiher_NIG$fisher_est)))[1:length(est_param[i,])]
  trace_random2[i,]<-NIGMVD_ass$mixedEffect_list$betar_vec
  
  cat('r = ',est_param2[i,3],'\n')
  cat('r = ',est_param2[i,3],'\n')
  }
  cat('low = ',est_param[i,3] - 1.6*sd_param[i,3]-true_param[3],'\n')
  cat('upp = ',est_param[i,3] + 1.6*sd_param[i,3]-true_param[3],'\n')
}
res <- sum_res(est_param,sd_param, true_param)
cat(colMeans(res$in_CI),'\n')
x11()
par(mfrow=c(3,1))
for(i in 1:3){
  plot(res$lower[,i] - true_param[i], type='l',col='red', ylim = c(min(res$lower[,i]-true_param[i]),max(res$upper[,i]-true_param[i])) )
  lines(res$upper[,i] - true_param[i],col='red')
  abline(h=0,col='blue')
}
if(0){
res <- sum_res(est_param2,sd_param2, true_param)
cat(colMeans(res$in_CI))
x11()
par(mfrow=c(3,1))
for(i in 1:3){
  plot(res$lower[,i] - true_param[i], type='l',col='red', ylim = c(min(res$lower[,i]-true_param[i]),max(res$upper[,i]-true_param[i])) )
  lines(res$upper[,i] - true_param[i],col='red')
  abline(h=0,col='blue')
}
}
x11()
par(mfrow=c(2,1))
plot(trace_random[1,])
#plot(trace_random2[1,])
