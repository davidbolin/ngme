library(testthat)
graphics.off()
library(ngme)
source('CharFunc.R')
source('charfuncUtil.R')

kappa = .8
alpha    = 3.
d     = 1
tau   = 1
MaternParameter = c(alpha, tau, kappa)

delta = 1
mu    = 1
sigma = 0.5
nu    = 3
NIGParameter = c(delta, mu, sigma, nu)

x = seq(-20,20,length.out = 2000)

f <-maternkernel(x, MaternParameter[1],  MaternParameter[2] , MaternParameter[3], d)


sim <- 5000

#x11()
#hist(X,probability = T)
#lines(dens$x, dens$density,type='l')
test_that("1d normal matern", {
  
  dens <- density_1d_normal(c(-16,16),
                            11,
                            param = c(delta, sigma), 
                            maternParam= MaternParameter)
  X <- rep(0,sim)
  for(i in 1:sim){
    X[i] <- sqrt(x[2]-x[1]) * sigma * sum( f * (sqrt(x[2]-x[1])*delta/sigma + rnorm(length(x))))
  }
  Moment <- calc_moment_density(dens)
  Moment_emp <- calc_moment_emperical(X)
  
    
    expect_equal(Moment[1] - Moment_emp[1],0 ,tolerance=10^-1)
    expect_equal(Moment[2] - Moment_emp[2],0 ,tolerance=10^-1)
    expect_equal(Moment[3] - Moment_emp[3],0 ,tolerance=10^-1)
    expect_equal(Moment[4] - Moment_emp[4],0 ,tolerance=0.5)
})


#cat('emperical -  moment[X] = ',round(Moment_emp-Moment,2),'\n')
#x11()
#corr <- materncorr(x, MaternParameter[1],  MaternParameter[2],  1)
#plot(x, corr,type='l')

test_that("1d nig matern", {
X <- rep(0,sim)
for(i in 1:sim){
  W <- ngme::rNIG(1, 
                  NIGParameter[1],
                  NIGParameter[2], 
                  NIGParameter[3], 
                  NIGParameter[4],
                  rep(x[2]-x[1], length(f)))
  X[i] <-  sum( f *W)
}

dens <- density_1d_nig(c(-20,20),
                          12,
                          param =NIGParameter, 
                          maternParam= MaternParameter)


Moment <- calc_moment_density(dens)
Moment_emp <- calc_moment_emperical(X)

expect_equal(Moment[1] - Moment_emp[1],0 ,tolerance=0.1)
expect_equal(Moment[2] - Moment_emp[2],0 ,tolerance=0.1)
expect_equal(Moment[3] - Moment_emp[3],0 ,tolerance=0.2)
expect_equal(Moment[4] - Moment_emp[4],0 ,tolerance=0.5)
})


#x11()
#hist(X,probability = T)
#lines(dens$x, dens$density,type='l')
#cat('emperical - moment[X] = ',round(Moment_emp-Moment,2),'\n')