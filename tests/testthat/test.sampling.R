#seem to work

rm(list = ls())
chisquare_test <- function(m, Femp, EY)
{
  n <- length(EY)
  chi2 <- sum((EY-m*(Fn(0:(n-1))-c(0,Fn(0:(n-2)))))^2/EY)
  return(chi2)
}


set.seed(1)
library(testthat)
library(ngme)
n <- 1000
sim <- 1000

m <- 100
w <- rgamma(n, 10, 1)
w <- w/sum(w)
Y <- c()
for(i in 1:sim)
  Y <- c(Y, sampleR(m, w))
res<-hist(Y,breaks = 0:999, plot = F)
Fn <- ecdf(Y)
#nice brownian bridge
#plot(Fn(0:999)-cumsum(w))

test_that("max, min", {
  expect_equal(max(Y), n-1)
  expect_equal(min(Y), 0)
})
test_that("sample chisquare", {
  expect_lt(chisquare_test(length(Y),Fn,length(Y)*w), qchisq(0.95,n))
})
selected_in <- rep(0,n)
weight <- rep(0,n)
ans <- sample_internalR(m, 
                       w,
                       selected_in,
                       weight)
test_that("sample weighted, check weight is correct", {
  ans$w_in[ans$w_in>0]-1/(m*w[ans$w_in>0])
})
Y <- rep(0,n)
for(i in 1:sim)
{
  ans <- sample_internalR(m, 
                          w,
                          selected_in,
                          weight)
  
 Y[ans$ans + 1] = Y[ans$ans + 1] + 1 
}
# have not figured out good test yet
test_that("test stupid", {
expect_equal(mean(w*m*sim-Y) <1 ,TRUE)
})

ans2 <- sample_internalR(6*m, 
                        w,
                        selected_in,
                        weight)

test_that("sample weighted 1", {
  expect_equal(ans2$w_in[which.max((w))], 1)
})
index = order(w,decreasing = T)
w_mod <- 6*m*w
for(i in 1:n){
  if(w_mod[index[i]] < 1)
    break
  
  p_ <- w_mod[index[i]]
  w_mod[index[(i+1):n]] <- w_mod[index[(i+1):n]] +  + (p_-1) * w[index[(i+1):n]]
  w_mod[index[i]] <- 1
                            
}
w_mod[index[i:n]] <- 1/w_mod[index[i:n]]
test_that("test weights", {
  expect_equal(max(abs(w_mod[ans2$w_in>0] - ans2$w_in[ans2$w_in>0])), 0)
})

Y <- rep(0,n)
selected_in[sample(1:n,10)] <- 1
w_selected <- rep(0, 10)
for(i in 1:sim)
{
  ans <- sample_internalR(m, 
                          w,
                          selected_in,
                          weight)
  
  Y[ans$ans + 1] = Y[ans$ans + 1] + 1 
  w_selected <- w_selected + ans$w_in[selected_in==1]
}
test_that("selected should not be redrawn", {
  expect_equal(sum(Y[selected_in==1]),0)
})
test_that("selected should have weights", {
  expect_gt(sum(w_selected),0)
})

selected_in <- rep(0,n)
Y <- sampleR(m, w)
selected_in[Y+1] = 1
ans <- sample_internalR(m, 
                        w,
                        selected_in,
                        weight)
test_that("no interssection allowed", {
  expect_equal(length(intersect(Y, ans$ans)),0)
})
Y_emp <- rep(0,n)
for(i in 1:sim)
{
Y <- sampleR(m, w)
selected_in <- rep(0,n)
selected_in[Y+1] = 1
ans <- sample_internalR(m, 
                        w,
                        selected_in,
                        weight)
Y_emp[Y+1] = Y_emp[Y+1] + 1
Y_emp[ans$ans+1] = Y_emp[ans$ans+1] + 1
}