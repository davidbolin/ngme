rm(list=ls())
graphics.off()
library(LDMod)
library(testthat)
library(LDMod)
library(MASS)
seed     <- 3
nsamples <- 5000
burnin <- ceiling(nsamples*0.1)
silent   <- 1
plotflag <- 0
d <- 5
n.pers <- 2 #number of patients
n.obs  <- 50 #number of obs per patient
A <- matrix(rnorm(d*d), ncol= d, nrow=d)
COV_ <-  t(A)%*%A
eig <- eigen(COV_)
eig$values[5] <- 10^-2.5  *eig$values[5]
COV_beta <- eig$vectors %*% diag(eig$values) %*%t(eig$vectors)
sd_Y    <- 0.1 # error of the noise

Br_list <- list()
betar <- rep(0,d)
betaf <- c(1.)
mu   <- c(0.6, -0.6, rnorm(d-2))
nu <- 10
betar_list <- list()
Bf_list    <- list()
V_list     <- list()
Y_list     <- list()
set.seed(seed)
for(i in 1:n.pers)
{
  Bf_list[[i]]    <- as.matrix(runif(n = n.obs))
  Br_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs , matrix(rnorm(n.obs*(d-2)),nrow=n.obs, ncol=d-2 ))
  V <- LDMod::rGIG(-0.5, nu, nu, sample.int(10^6,1))
  V_list[[i]] <- V
  betar_list[[i]] <- betar - mu * 1  + V * mu +
    sqrt(V) * mvrnorm(n = 1, mu  =rep(0,d), Sigma = COV_beta)
  Y_list[[i]]        <- rnorm(n = n.obs,
                              Br_list[[i]]%*%(betar_list[[i]] - betar), sd = sd_Y)
  
}


meas_list <- list(Y = Y_list, sigma_eps = sd_Y, noise = "Normal")
mixedEffect_list <- list(B_random = Br_list,
                         B_fixed  = Bf_list,
                         Sigma = COV_beta,
                         beta_random = rep(0,d),
                         beta_fixed  = c(0.),
                         mu          = as.matrix(mu),
                         nu          = as.matrix(1.) + nu,
                         noise = "NIG")

out = test_sampling_NIG(mixedEffect_list,
                  meas_list,
                  nsamples)
x11()
par(mfrow=c(2, 3))
plot(out$U_MALA[burnin:nsamples ,1], xlab="i", ylab="U", main = "MALA")
acf(out$U_MALA[burnin:nsamples  ,1],  main = "MALA")
plot(density(out$U_MALA[burnin:nsamples ,1]), main="MALA")
plot(out$U_Gibbs[burnin:nsamples  ,1], xlab="i", ylab="U", main = "Gibbs")
acf(out$U_Gibbs[burnin:nsamples  ,1], main = "Gibbs")
plot(density(out$U_Gibbs[burnin:nsamples ,1]), main="Gibbs")
cat("acc = ", out$acc_MALA, '\n')
x11()
X_MALA  <- colSums(diag(eig$vectors[,5])%*%t(out$U_MALA[burnin:nsamples,]))
X_Gibbs <- colSums(diag(eig$vectors[,5])%*%t(out$U_Gibbs[burnin:nsamples,]))
par(mfrow=c(2, 3))
plot(X_MALA, xlab="i", ylab="U", main = "MALA")
acf(X_MALA,  main = "MALA")
plot(density(X_MALA), main="MALA")
plot(X_Gibbs, xlab="i", ylab="U", main = "Gibbs")
acf(X_Gibbs, main = "Gibbs")
plot(density(X_Gibbs), main="Gibbs")
