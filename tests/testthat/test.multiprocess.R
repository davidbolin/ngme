####
#  testing estimation (and sampling of multiprocess)
#  simulating to types of variance processes
# D: 2019-04-02
####
library(testthat)
library(ngme)
context("process nig")
set.seed(12)
niter = 200
sim <- 5
nu <- c(0.5, 5)
mu <- c(-2,5) #c(-2,  10)
n     <- 80
nindv <- 20
h        <- list()
Vprocess <- list()
U        <- list()

Bnu <- matrix(0, nrow=2*n,ncol=2 )
Bnu[1:n         ,1] <- 1
Bnu[(n+1):(2*n) ,2] <- 1 
B <- list()
for(i in 1:nindv){
  h[[i]] <- rep(1,2*n)#runif(2*n)+0.1
  B[[i]] <- Bnu
  Vprocess[[i]]      <-rGIG(rep(-0.5,2*n),
                            B[[i]]%*%nu,
                            (B[[i]]%*%nu) * h[[i]]^2,
                            as.integer(1000 * runif(1) ))
  U[[i]] <-  (-h[[i]] + Vprocess[[i]])*(B[[i]]%*%mu)+ sqrt(Vprocess[[i]] ) * rnorm(2*n)
}
process_list_1 <- list(V = Vprocess, X = U, Bmu = B, Bnu = B, h = h, noise= 'NIG', nu = nu[1])
for(i in 1:length(process_list_1$V)){
  process_list_1$V[[i]] = process_list_1$V[[i]][1:n]
  process_list_1$X[[i]] = process_list_1$X[[i]][1:n]
  process_list_1$h[[i]] = process_list_1$h[[i]][1:n]
}

process_list <- list(V = Vprocess, X = U, Bmu = B, Bnu = B, h = h, noise= 'NIG', sampleV =F)
out <- test_Mprocess(niter,
                     U,
                     process_list)

process_list <- list(V = Vprocess, X = U, Bmu = B, Bnu = B, h = h, noise= 'NIG', sampleV =T)
out2 <- test_Mprocess(niter,
                     U,
                     process_list)



test_that("NIG multi process nu", {
  expect_equal(nu/nu,
               exp(c(out$process_list$nu))/nu,
               tolerance = 0.1)
})
test_that("NIG multi process mu", {
  expect_equal(mu/abs(mu),
               c(out$process_list$mu)/abs(mu),
               tolerance = 0.01)
})
test_that("NIG multi process nu", {
  expect_equal(mu/abs(mu),
               c(out2$process_list$mu)/abs(mu),
               tolerance = 0.1)
})
