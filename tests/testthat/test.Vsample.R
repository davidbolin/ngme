####
#  testing GIG (nig) sampling
# D: 2019-04-02
####
rm(list=ls())
graphics.off()
library(ngme)

n<- 500
a <- 1.1
b <- 0.3
V      <-rGIG(rep(-0.5,n),
              rep(a, n),
              rep(b, n),
              as.integer(1000 * runif(1) ))
x11()
hist(V,100,probability = T)
x <- seq(min(V),max(V),length.out = 1000)
lines(x, dig(x, a, b), col='red')