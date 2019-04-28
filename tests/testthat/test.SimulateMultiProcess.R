graphics.off()
library(ngme)
library(testthat)
library(INLA)

use.process = TRUE
estimate.parameters = FALSE
#data options
n.obs=200 #number of observations per replicate
n.rep = 3 #number of replicates

cutoff = 0.1
max.dist = 1
n.lattice <- 50
kappa1 = 1
kappa2 = 1
tau1 = 5
tau2 = 5
rho = 1
theta = 0
nu <- c(log(10), log(10))
mu <- c(0     , 3)
sigma.e = c(0.01,0.1)
beta.fixed = c(0,0)


#create fixed effects-list with one intercept per field
Bf <- kronecker(diag(2),matrix(rep(1, n.obs)))
B_fixed <- list()
for(i in 1:n.rep){
  B_fixed[[i]]  <- Bf
}
mixedEffect_list  <- list(B_fixed  = B_fixed,
                          beta_fixed  = as.matrix(beta.fixed),
                          noise = "Normal")

#create measurement error list, with one sigma per field
B.e <- list()
for(i in 1:n.rep){
  B.e[[i]] <- kronecker(diag(2),matrix(rep(1, n.obs)))
}

mError_list <- list(noise = "nsNormal", 
                    B = B.e,
                    theta = matrix(log(sigma.e)))



#create mesh 
x=seq(from=0,to=10,length.out=n.lattice)
lattice=inla.mesh.lattice(x=x,y=x)
mesh=inla.mesh.create(lattice=lattice, extend=FALSE, refine=FALSE)
operator_list <- create_operator_matern2Dbivariate(mesh)

#create observation locations, both fields are observed at the same locations
obs.loc <- list()
for(i in 1:n.rep){
  obs.loc[[i]] = cbind(runif(n.obs)*diff(range(x))+min(x),
                       runif(n.obs)*diff(range(x))+min(x))
  
}

###
# simulation
#
###


operator_list$kappa1 <- kappa1
operator_list$kappa2 <- kappa2
operator_list$tau1   <- tau1
operator_list$tau2   <- tau2
operator_list$rho    <- rho
operator_list$theta  <- theta



Bmu <- list()
Bnu <- list()
n.grid <- length(operator_list$h[[1]])/2
for(i in 1:n.rep){
  Bmu[[i]] <- kronecker(diag(2),matrix(rep(1, n.grid)))
  Bnu[[i]] <- kronecker(diag(2),matrix(rep(1, n.grid)))
}
#toSampleV
processes_list = list(noise = "MultiGH", 
                      Bmu = Bmu, mu = as.matrix(mu),
                      Bnu = Bnu, nu = as.matrix(nu))
sim_res <- simulateLongPrior( locs              = obs.loc,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = mError_list,
                              processes_list    = processes_list,
                              operator_list     = operator_list)
x11()
h <- operator_list$h[[1]][ceiling(n.grid/2)]
index <-  abs(operator_list$h[[1]] - h) < 10^-10 
index1 <- index
index2 <- index
index1[(n.grid+1):length(index)] <- F
index2[1:n.grid] <- F
E <- c(sim_res$Z[[1]][index1],sim_res$Z[[2]][index1],sim_res$Z[[3]][index1])
hist(E, probability = T,200)
x_min <- min(E)
x_max <- max(E)
x_g <- seq(x_min, x_max, length.out = 1000) 
f_x <- dnigMeasure(x_g,h,-mu[1]*h, mu[1], exp(nu[1]), 1,log=F)
lines(x_g,f_x,col='red')
x11()
E <- c(sim_res$Z[[1]][index2],sim_res$Z[[2]][index2],sim_res$Z[[3]][index2])
hist(E, probability = T,200)
x_min <- min(E)
x_max <- max(E)
x_g <- seq(x_min, x_max, length.out = 1000) 
f_x <- dnigMeasure(x_g,h,-mu[2]*h, mu[2], exp(nu[2]), 1,log=F)
lines(x_g,f_x,col='red')