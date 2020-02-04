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
operator_list$rho   <- 0.6

operator_list$theta   <- 0.6
process_list = list(noise = "Normal",
                    nu  = 1,
                    mu  = 0)
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
rho = operator_list$rho 
theta  =  operator_list$theta
B = c(cos(theta)+rho*sin(theta) , -sin(theta)*sqrt(1+rho^2)  ,
      sin(theta)-rho*cos(theta),cos(theta)*sqrt(1+rho^2))
B = t(matrix(B, nrow=2,ncol=2))
maternParam <- list()
maternParam[[1]] <- c(2, operator_list$tau1, operator_list$kappa1)
maternParam[[2]] <- c(2, operator_list$tau2, operator_list$kappa2)
maternParam[[3]] <- c(operator_list$rho,operator_list$theta)
res<- density_2d_normal_multivariate(c(-3,3),
                                     c(-5,5),
                                     7,
                                     6,
                                     c(0,1),
                                     c(0,1), 
                                     maternParam)



Moment_emp1 <- calc_moment_emperical(X1s[,1] )
Moment_emp2 <- calc_moment_emperical(X2s[,1] )
cov_emp    <- cov(cbind(X1s[,1],X2s[,1]))[1,2]
moment2 <- calc_moment_density_2d(res)

test_that("2d multi normal matern", {
  expect_equal(cov_emp,moment2[[3]],tolerance=0.01)
  f_1 = function(x){maternkernelMulti(x,2,operator_list$tau1,operator_list$kappa1,2)}
  f_2 = function(x){maternkernelMulti(x,2,operator_list$tau2,operator_list$kappa2,2)}
  VX1 <- integrate(f_1,0,20)$value*solve(B)[1,1]^2+integrate(f_1,0,20)$value*solve(B)[1,2]^2
  VX2 <- integrate(f_2,0,20)$value*solve(B)[2,2]^2+integrate(f_2,0,20)$value*solve(B)[2,1]^2
  expect_equal(moment2[[1]][2] -VX1,0 ,tolerance=10^-5)
  expect_equal(moment2[[2]][2] -VX2,0 ,tolerance=10^-5)
  expect_equal(moment2[[1]][1] - Moment_emp1[1],0 ,tolerance=10^-1)
  expect_equal(moment2[[1]][2] - Moment_emp1[2],0 ,tolerance=10^-1)
  expect_equal(moment2[[1]][3] - Moment_emp1[3],0 ,tolerance=2*10^-1)
  expect_equal(moment2[[1]][4] - Moment_emp1[4],0 ,tolerance=5*10^-1)
  expect_equal(moment2[[2]][1] - Moment_emp2[1],0 ,tolerance=10^-1)
  expect_equal(moment2[[2]][2] - Moment_emp2[2],0 ,tolerance=10^-1)
  expect_equal(moment2[[2]][3] - Moment_emp2[3],0 ,tolerance=2*10^-1)
  expect_equal(moment2[[2]][4] - Moment_emp2[4],0 ,tolerance=5*10^-1)
})

if(0){
  res<- density_2d_normal_multivariate(c(-3,3),
                                       c(-5,5),
                                       7,
                                       6,
                                       c(0,1),
                                       c(0,1), 
                                       maternParam)
x11()
par(mfrow=c(3,1))
hist(X1s[,1],prob=T,breaks=20)
lines(res$x,colSums(res$density)*(res$y[2]-res$y[1]))
hist(X2s[,1],prob=T,breaks=20)
lines(res$y,rowSums(res$density)*(res$x[2]-res$x[1]))
grid <- meshgrid(res$x,res$y)
image(res$y,res$x,res$density,xlim = c(-2,2), ylim = c(-0.4,0.4))
}

