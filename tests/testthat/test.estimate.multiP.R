graphics.off()
library(ngme)
library(testthat)
library(INLA)
nIter = 30000

use.process = TRUE
estimate.parameters = FALSE
#data options
n.obs=200 #number of observations per replicate
n.rep = 10 #number of replicates

cutoff = 0.1
max.dist = 1
n.lattice <- 10
kappa1 = 1
kappa2 = 1
tau1 = 5
tau2 = 5
rho = 0
theta = 0
nu <- c(log(2), log(2))
mu <- c(3     , 3)
sigma.e = c(0.01,0.01)
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
print(processes_list$noise)

#processes_list = list(noise = "MultiGH", 
#                      Bmu = Bmu, mu = 0*as.matrix(mu),
#                      Bnu = Bnu, nu = 0*as.matrix(nu)
#                      )
operator_list$kappa1 <- 2
operator_list$kappa2 <- 2
operator_list$tau1   <- 4
operator_list$tau2   <- 4
operator_list$rho    <- rho
operator_list$theta  <- theta
res.est <- estimateLong(Y                = sim_res$Y,
                        nIter            = nIter,
                        nSim             = 2,
                        step0            = 0.5,
                        alpha            = 0.01,
                        locs             = obs.loc,
                        mixedEffect_list = mixedEffect_list,
                        measurment_list  = mError_list,
                        processes_list   = processes_list,
                        operator_list    = operator_list,
                        nBurnin_learningrate = 50,
                        silent = FALSE)