library(testthat)
source('charfuncUtil.R')
library(mvtnorm)
d <- characteristic_function_to_density(
  function(t,mu=1,sigma=.5) 
     {1i*t*mu - sigma^2/2*t^2 },
  8,
  c(-4,4)
)

test_that("univarate normal phi", {
  
  expect_equal( sum(d$density-dnorm(d$x,1,.5)),0 ,tolerance=10^-4)
})


mu <- c(0,2)
Sigma <-  matrix(c(1,1,1,3), nrow=2)
phi_multi <- function(t_x, 
                      t_y, 
                      mu, 
                      Sigma
                      ){
  lphi <- matrix(0, nrow=length(t_x), ncol=length(t_y))
  for(i_x in 1:length(t_x) ){
    for(i_y in 1:length(t_y) ){
      t_ <- c(t_x[i_x], t_y[i_y])
      lphi[i_x, i_y] <- 1i*t(mu)%*%t_ - 0.5  * t(t_)%*%Sigma%*%t_
    }
  }
  return(lphi)
}
lphi <- function(t_x,t_y){phi_multi(t_x,t_y, mu, Sigma)}
dens <- characteristic_function_to_density2d(lphi, 
                                             7,
                                             7,
                                             c(-10,10),
                                             c(-10,10))



grid <- meshgrid(dens$x,dens$y)

dens_true <- dmvnorm(cbind(c(grid$X),c(grid$Y)), mean=mu, sigma=Sigma, log=FALSE)

test_that("multivariate normal phi", {
  
  expect_equal(sum(abs(dens_true-c(t(dens$density)))),0 ,tolerance=10^-4)
})