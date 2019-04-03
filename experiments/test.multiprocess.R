####
#  testing estimation (and sampling of multiprocess)
#  simulating to types of variance processes
# D: 2019-04-02
####
rm(list=ls())
graphics.off()
library(ngme)
set.seed(112)
niter = 30
nu <- c(.5, 5)
mu <- c(0,0) #c(-2,  10)
n     <- 500
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
  U[[i]] <- (Vprocess[[i]]  - h[[i]]) * (B[[i]]%*%mu) + sqrt(Vprocess[[i]] ) * rnorm(2*n,1)
}
process_list_1 <- list(V = Vprocess, X = U, Bmu = B, Bnu = B, h = h, noise= 'NIG')
for(i in 1:length(process_list_1$V)){
  process_list_1$V[[i]] = process_list_1$V[[i]][1:n]
  process_list_1$X[[i]] = process_list_1$X[[i]][1:n]
  process_list_1$h[[i]] = process_list_1$h[[i]][1:n]
}

##
##
nu_est <- 1
nu_vec <- c()
for(j in 1:niter){
  dnu = 0
  ddnu = 0
  for(i in 1:nindv){
    X = process_list_1$X[[i]]
    Vpost <- rGIG( rep(-1, n),
                   rep(nu_est,n),
                   X^2 + nu_est * process_list_1$h[[i]]^2,
                   as.integer(1000 * runif(1)))
    dnu = dnu   + 0.5 * (n/nu_est - sum(process_list_1$h[[i]] /Vpost) - sum(Vpost) + sum(process_list_1$h[[i]]))
    ddnu = ddnu - 0.5 * n/nu_est^2
    
  }
  nu_est = nu_est - 0.5*dnu/ddnu
  
  nu_vec <- c(nu_vec, nu_est)
}


##
##
process_list <- list(V = Vprocess, X = U, Bmu = B, Bnu = B, h = h, noise= 'NIG')
out <- test_Mprocess(niter,
                     U,
                     process_list)

out2 <- test_process(niter,
                      process_list_1$U,
                      process_list_1)
x11()
par(mfrow=c(2,1))
ymax <- max(max(exp(out$process_list$nu_vec[,1])),out2$process_list$nu_vec[,1])
plot(exp(out$process_list$nu_vec[,1]),ylim=c(0,ymax))
lines(out2$process_list$nu_vec[,1],col='red')
lines(nu_vec,col='green')

plot(exp(out$process_list$nu_vec[,2]))
x11()
plot(process_list_1$V[[1]],out2$process_list$V[[1]])

x11()
V <- unlist(out2$process_list$V)
hist(V,100,probability = T)
x <- seq(min(V),max(V),length.out = 1000)
lines(x, dgig(x,-1,0.717451,4.71745), col='red')