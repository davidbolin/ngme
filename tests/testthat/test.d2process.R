rm(list=ls())
require(testthat)
context("Fisher")
library(LDMod)
library(MASS)

test_that("Fisher, d2 latent process", {
library(LDMod)
library(MASS)

n.obs <- 20
n = n.obs
nu = 2.12
mu = -0.1
sd_Y = 1
locs <- list()
locs[[1]] <- 1:n.obs

mixedEffect_list  <- list(noise = "Normal")
operator_list <- create_operator(locs, n, name = "fd2")
Ylist <- list(c(rnorm(n.obs))) 


processes_list = list(nu = nu, mu = mu, noise = "NIG")

sim_res <- simulateLongPrior( Y                 = Ylist,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = list(sigma = sd_Y, noise = "Normal"),
                              processes_list    = processes_list,
                              operator_list     = operator_list)
processes_list$V <- sim_res$V
processes_list$X <- sim_res$X

Y <- sim_res$Y
Y$Y <- Y[[1]]
Y$A <- sparseMatrix(1:n.obs,1:n.obs, x=rep(1,n.obs))
res <- test_d2_process(Y,
                processes_list,
                operator_list)

resid <- as.vector(operator_list$Q[[1]]%*%sim_res$X[[1]])
ddnu <-  0.5 * n / nu^2 
ddmu <- t((sim_res$V[[1]] - operator_list$h[[1]])/sim_res$V[[1]])%*%(+sim_res$V[[1]] - operator_list$h[[1]])
expect_equal(res$d2,matrix(c(ddmu,0,0,ddnu),nrow=2) , tolerance  = 100)
processes_list$noise = "GAL"
res <- test_d2_process(Y,
                       processes_list,
                       operator_list)


})