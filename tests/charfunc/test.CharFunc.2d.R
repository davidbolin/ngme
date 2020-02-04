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
operator_list <- create_operator_matern2D(mesh)  
operator_list$kappa <- 2
operator_list$tau   <- 1
process_list = list(noise = "Normal",
                    nu  = 1,
                    mu  = 0)
sim <- 4000
sim.res <- simulate.process(sim, 
                            operator_list = operator_list, 
                            process_list = process_list)
obs_loc = cbind(c(0,0.5,0,0.8),c(0,0.5,0.5,0.8))
A = build.A.matrix(operator_list, obs_loc)

Xs <- matrix(0, nrow=sim, ncol=dim(obs_loc)[1])
for(i in 1:sim){
  Xs[i,] <- as.vector(A%*%sim.res$X[[i]])
}
#eucldian distnace
x <- sqrt(rowSums((obs_loc)^2))
Ecorr <- materncorr(x, 2, operator_list$kappa, 2)

MaternParameter = c(2, operator_list$tau , operator_list$kappa)
dens <- density_2d_normal(c(-16,16),
                          11,
                          param = c(0,1), 
                          maternParam= MaternParameter)
Moment <- calc_moment_density(dens)
Moment_emp <- calc_moment_emperical(Xs[,1] )
test_that("2d normal matern", {
  expect_equal(cor(Xs)[1,] - Ecorr,rep(0,length(Ecorr)),tolerance=0.05)
  expect_equal(Moment[1] - Moment_emp[1],0 ,tolerance=10^-1)
  expect_equal(Moment[2] - Moment_emp[2],0 ,tolerance=10^-1)
  expect_equal(Moment[3] - Moment_emp[3],0 ,tolerance=2*10^-1)
  expect_equal(Moment[4] - Moment_emp[4],0 ,tolerance=0.5)
  
})