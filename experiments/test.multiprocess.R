####
#  testing estimation (and sampling of multiprocess)
#  simulating to types of variance processes
# D: 2019-04-02
####
rm(list=ls())
graphics.off()
library(ngme)
#set.seed(112)
niter = 200
sim <- 5
nu <- c(1, 1)
mu <- c(0,0) #c(-2,  10)
n     <- 8000
nindv <- 1
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
  nus <- B[[i]]%*%nu
  for(ii in 1:length(nus))
    Vprocess[[i]][ii] <-rgig(n=1,
                     lambda = -0.5, 
                     chi   = h[[i]][ii]^2*nus[ii],
                     psi   = nus[ii])
  #U[[i]] <- (Vprocess[[i]]  - h[[i]]) * (B[[i]]%*%mu) + sqrt(Vprocess[[i]] ) * rnorm(2*n,1)
  U[[i]] <-  sqrt(Vprocess[[i]] ) * rnorm(2*n)
}
process_list_1 <- list(V = Vprocess, X = U, Bmu = B, Bnu = B, h = h, noise= 'NIG', nu = nu[1])
for(i in 1:length(process_list_1$V)){
  process_list_1$V[[i]] = process_list_1$V[[i]][1:n]
  process_list_1$X[[i]] = process_list_1$X[[i]][1:n]
  process_list_1$h[[i]] = process_list_1$h[[i]][1:n]
}
####
## estimation viewing data
###
nu_est <- nu[1]
nu_vec2 <- c()
for(j in 1:niter){
  dnu = 0
  ddnu = 0
  for(i in 1:nindv){
    Vpost <-  process_list_1$V[[i]]
    dnu = dnu   + 0.5 * (n/nu_est - sum(process_list_1$h[[i]]^2 /Vpost) - sum(Vpost) + 2*sum(process_list_1$h[[i]]))
    ddnu = ddnu - 0.5 * n/nu_est^2
    
  }
  nu_est = - n/(- sum(process_list_1$h[[i]]^2 /Vpost) - sum(Vpost) + 2*sum(process_list_1$h[[i]]))
  #nu_est = nu_est - 0.2*dnu/ddnu
  
  nu_vec2 <- c(nu_vec2, nu_est)
}
x11()
V<-process_list_1$V[[1]]
hist(V,100,probability = T)
x <- seq(min(V),max(V),length.out = 1000)
lines(x, dig(x, nu[1], h[[1]][1]^2*nu[1]), col='red')
####
## estimation sampling
###
nu_est <- nu[1]
nu_vec <- c()
for(j in 1:niter){
  dnu = 0
  ddnu = 0
  for(i in 1:nindv){
    X = process_list_1$X[[i]]
    EiV <- rep(0,length(V))
    EV  <- rep(0,length(V))
    p <- -1
    a <- rep(nu_est[1],n)
    b <- X^2 + nu_est[1] * process_list_1$h[[i]]^2
    for(k in 1:sim){
    Vpost <-  rGIG( rep(p, n),
                    a,
                    b,
                    as.integer(1000 * runif(1)))
    EiV = EiV + 1/Vpost
    EV  = EV  + Vpost
    }
    EiV = sqrt(a/b) * besselK(sqrt(a*b), p + 1)/besselK(sqrt(a*b), p ) - 2*p/b
    EV  = sqrt(b/a) * besselK(sqrt(a*b), p + 1)/besselK(sqrt(a*b), p )
    #EV  = EV  / sim
    dnu = dnu   + 0.5 * (n/nu_est - sum(process_list_1$h[[i]]^2 * EiV) - sum(EV) + 2*sum(process_list_1$h[[i]]))
    ddnu = ddnu - 0.5 * n/nu_est^2
    
  }
  nu_est = -n/(- sum(process_list_1$h[[i]]^2 * EiV) - sum(EV) + 2*sum(process_list_1$h[[i]]))
  #nu_est = nu_est - 0.5*dnu/ddnu
  
  nu_vec <- c(nu_vec, nu_est)
}
nu_est <- nu[1]
nu_vec3 <- c()
for(j in 1:niter){
  dnu = 0
  ddnu = 0
  for(i in 1:nindv){
    X = process_list_1$X[[i]]
    Vpost <-  rGIG( rep(-1, n),
                    rep(nu_est,n),
                    X^2 + nu_est * process_list_1$h[[i]]^2,
                    as.integer(1000 * runif(1)))
    dnu = dnu   + 0.5 * (n/nu_est - sum(process_list_1$h[[i]]^2 /Vpost) - sum(Vpost) + 2*sum(process_list_1$h[[i]]))
    ddnu = ddnu - 0.5 * n/nu_est^2
    
  }
  nu_est = nu_est - 0.5*dnu/ddnu
  
  nu_vec3 <- c(nu_vec3, nu_est)
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
lines(nu_vec2,col='yellow')
lines(nu_vec3, col='blue')
plot(exp(out$process_list$nu_vec[,2]))
#x11()
#plot(process_list_1$V[[1]],out2$process_list$V[[1]])

#x11()
#V <- unlist(out2$process_list$V)
#hist(V,100,probability = T)
#x <- seq(min(V),max(V),length.out = 1000)
#lines(x, dgig(x,-1,0.717451,4.71745), col='red')