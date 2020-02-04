library(testthat)
library(ngme)
library(INLA)
source('CharFunc.R')
graphics.off()
n <- 100
locY = cbind(runif(n), runif(n))
mesh = inla.mesh.2d(loc=locY,
                    # the inner edge and outer edge
                    max.edge = c(0.5,1), 
                    # offset extension distance inner and outer extenstion
                    offset = c(1.5, 1.5)
)

operator_list <- create_operator_matern2Dbivariate(mesh)  
operator_list$kappa1 <- 2.2
operator_list$kappa2 <- 1.3
operator_list$tau1   <- 2
operator_list$tau2   <- 0.4
operator_list$rho   <- -0.9

operator_list$theta   <- -0.4
process_list = list(noise = "MultiGH",
                    nu  = as.matrix(c(1,log(2))),
                    mu  = as.matrix(c(2,-2)))
sim <- 4000
sim.res <- simulate.process(sim, 
                            operator_list = operator_list, 
                            process_list = process_list)

obs_loc = cbind(c(0,0.5,0,0.8),c(0,0.5,0.5,0.8))
A = build.A.matrix(operator_list, obs_loc)

X1s <- matrix(0, nrow=sim, ncol=dim(obs_loc)[1])
X2s <- matrix(0, nrow=sim, ncol=dim(obs_loc)[1])
for(i in 1:sim){
  X1s[i,] <- as.vector(A%*%sim.res$X[[i]][1:mesh$n])
  X2s[i,] <- as.vector(A%*%sim.res$X[[i]][mesh$n + (1:mesh$n)])
}
maternParam <- list()
maternParam[[1]] <- c(2, operator_list$tau1, operator_list$kappa1)
maternParam[[2]] <- c(2, operator_list$tau2, operator_list$kappa2)
maternParam[[3]] <- c(operator_list$rho,operator_list$theta)
res<- density_2d_nig_multivariate(c(-3,3),
                                     c(-5,5),
                                     7,
                                     6,
                                     c(0,process_list$mu[1], 1, exp(process_list$nu[1])),
                                     c(0,process_list$mu[2], 1, exp(process_list$nu[2])), 
                                     maternParam)


Moment_emp1 <- calc_moment_emperical(X1s[,1] )
Moment_emp2 <- calc_moment_emperical(X2s[,1] )
cov_emp    <- cov(cbind(X1s[,1],X2s[,1]))[1,2]
moment2 <- calc_moment_density_2d(res)

test_that("2d multi nig matern", {
  expect_equal(cov_emp,moment2[[3]],tolerance=0.1)
  expect_equal(moment2[[1]][1] - Moment_emp1[1],0 ,tolerance=0.1)
  expect_equal(moment2[[1]][2] - Moment_emp1[2],0 ,tolerance=0.2)
  expect_equal(moment2[[1]][3] - Moment_emp1[3],0 ,tolerance=0.5)
  expect_equal(moment2[[1]][4] - Moment_emp1[4],0 ,tolerance=2)
  expect_equal(moment2[[2]][1] - Moment_emp2[1],0 ,tolerance=0.1)
  expect_equal(moment2[[2]][2] - Moment_emp2[2],0 ,tolerance=0.2)
  expect_equal(moment2[[2]][3] - Moment_emp2[3],0 ,tolerance=0.5)
  expect_equal(moment2[[2]][4] - Moment_emp2[4],0 ,tolerance=1)
  f_1 = function(x){2*pi*x*maternkernelMulti(x,2,operator_list$tau1,operator_list$kappa1,2)^2}
  f_2 = function(x){2*pi*x*maternkernelMulti(x,2,operator_list$tau2,operator_list$kappa2,2)^2}
})
if(0){
  x11()
  par(mfrow=c(3,1))
  hist(X1s[,1],prob=T,breaks=50)
  lines(res$x,colSums(res$density)*(res$y[2]-res$y[1]))
  hist(X2s[,1],prob=T,breaks=50)
  lines(res$y,rowSums(res$density)*(res$x[2]-res$x[1]))
  
  grid <- meshgrid(res$x,res$y)
  image(res$y,res$x,res$density,xlim = c(-2,2), ylim = c(-0.4,0.4))
}