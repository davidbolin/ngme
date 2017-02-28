rm(list = ls())
graphics.off()
library(LDMod)
library(INLA)
library(fields)
#First estimate stationary model:
noise="Gaussian"
kappa2 = 2
phi2 = 0.5
sigma2.e = 0.01

sigma = 1
lambda = 20
delta = 0
mu = 0
theta.mu = 1
n.lattice = 50
n.obs=100

x=seq(from=0,to=10,length.out=n.lattice)
lattice=inla.mesh.lattice(x=x,y=x)
mesh=inla.mesh.create(lattice=lattice, extend=FALSE, refine=FALSE)
obs.loc = cbind(runif(n.obs)*diff(range(x))+min(x),
                runif(n.obs)*diff(range(x))+min(x))


Vin <- list(rep(1, n.obs))
B_fixed  <- list(matrix(rep(1, n.obs)))

mError_list <- list(noise = "Normal", sigma = 0.1)
mixedEffect_list  <- list(B_fixed  = B_fixed,
                          beta_fixed  = as.matrix(c(1)),
                          noise = "Normal")

operator_list <- create_operator_matern2D(mesh)

operator_list$kappa <- 2
operator_list$tau   <- 15

processes_list = list(noise = "Normal",nu  = 0., mu  = 0., V <- list())
processes_list$V[[1]] <- operator_list$h

sim_res <- simulateLongPrior( locs              = list(obs.loc),
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = mError_list,
                              processes_list    = processes_list,
                              operator_list     = operator_list)

proj <- inla.mesh.projector(mesh,dims=c(80,80))
image.plot(proj$x,proj$y,inla.mesh.project(proj,sim_res$X[[1]]),xlab="",ylab="")



