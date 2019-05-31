##
#  test mixed effect p measurement error
##
rm(list=ls())
graphics.off()
library(ngme)
set.seed(7)
save.file=0
sim <-  100

nindv <- 500
n     <- 100

beta_random  <- as.vector(0.8)
beta_fixed   <- c(1.1, 2.2)
sigma        <- 0.5
sigma_random <- 0.5
nu_mixed     <- 1
mu_mixed     <- 2
est_param <- est_param2 <- matrix(0, nrow=sim, ncol = 3)
sd_param  <- sd_param2 <- matrix(0, nrow=sim, ncol = 3)
sd_param2  <- sd_param2 <- matrix(0, nrow=sim, ncol = 3)
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

m_ = 0
m_2 = 0
for(i in 1:2){
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
  
  if(nindv >= 100){
    control = list(estimate.fisher = FALSE,
         pSubsample = 0.2,
         subsample.type  = 1,
         alpha = 0.3,
         nSim  =2,
         nBurnin = 2,
         step0 = 1,
         nBurnin.learningrate=1000)
  }else{
    control = list(estimate.fisher = FALSE,
         pSubsample =1,
         subsample.type  = 1,
         step0 = 1,
         alpha = 0.3,
         nSim  =2,
         nBurnin = 2,
         nBurnin.learningrate=1000)
  }
  NIGMVD_ass <-  ngme.par(n.cores = 4, std.lim = 2, max.rep = 6,
                                fixed       = Ya ~ B1 + B2,
                                random      = ~ -1+B3|id,
                                data        = as.data.frame(data),
                                reffects    = 'Normal',
                                use.process = FALSE,
                                silent      = FALSE,
                                controls.init = list(nIter.init=500),
                                nIter  = 1000,
                                controls    = control)
  
  fiher_NIG <- ngme.fisher(NIGMVD_ass,
                           nSim = 30,
                           nIter = 10,
                           nBurnin = 2,
                           n.cores = 1,
                           n.rep = 1,
                           std.threshold = 2,
                           observed = TRUE,
                           only.effects=F,
                           silent = F)
  est_param[i,] <- c(NIGMVD_ass$mixedEffect_list$beta_fixed, NIGMVD_ass$mixedEffect_list$beta_random)
  V <- diag(solve(fiher_NIG$fisher_est[c(1:3), c(1:3)]))[1:3]
  sd_param[i,]  <-  sqrt(V + NIGMVD_ass$fixed_est_var[dim(NIGMVD_ass$fixed_est_var)[1],c(2,3,1)])
  sd_param2[i,]  <-  sqrt(V)
  
 
 # cat('mu  = ',round(NIGMVD_ass$mixedEffect_list$mu,2),'\n')
  cat('b_r = ',round(NIGMVD_ass$mixedEffect_list$beta_random,2),'\n')
  cat('low = ',est_param[i,3] - 1.64*sd_param[i,3]-true_param[3],'\n')
  cat('upp = ',est_param[i,3] + 1.6*sd_param[i,3]-true_param[3],'\n')
  m_ = m_ + ((est_param[i,3] - 1.64*sd_param[i,3]-true_param[3]) * (est_param[i,3] + 1.64*sd_param[i,3]-true_param[3])<0)
  m_2 = m_2 + ((est_param[i,3] - 1.64*sd_param2[i,3]-true_param[3]) * (est_param[i,3] + 1.64*sd_param2[i,3]-true_param[3])<0)
  
  cat('CI  = ',100*round(m_/i,2),'%\n')
  cat('CI_v2  = ',100*round(m_2/i,2),'%\n')
}
res <- sum_res(est_param,sd_param, true_param)
cat(colMeans(res$in_CI),'\n')
if(save.file){
if(mu_mixed ==0){
  file.name = paste('res_',nindv,'_mu0.RData',sep="")
}else{
  file.name = paste('res_',nindv,'_mu1.RData',sep="")
}
save(res, file = file.name)
}
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
