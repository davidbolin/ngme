rm(list = ls())
graphics.off()
library(ngme)
library(INLA)
library(fields)
#First estimate stationary model:
test.pred = FALSE
test.est = TRUE
nIter = 1000

noise="Gaussian"
kappa1 = 1
kappa2 = 1
tau1 = 5
tau2 = 5
rho = 0
theta = 0
sigma.e = c(0.01,0.1)
beta.fixed = c(0,0)

n.lattice = 40
n.obs=1000 #number of observations per replicate
n.rep = 1 #number of replicates

#create mesh 
x=seq(from=0,to=10,length.out=n.lattice)
lattice=inla.mesh.lattice(x=x,y=x)
mesh=inla.mesh.create(lattice=lattice, extend=FALSE, refine=FALSE)

#create observation locations, both fields are observed at the same locations
obs.loc = cbind(runif(n.obs)*diff(range(x))+min(x),
                runif(n.obs)*diff(range(x))+min(x))


#create fixed effects-list with one intercept per field
Bf <- kronecker(diag(2),matrix(rep(1, n.obs)))
B_fixed  <- list(Bf)
mixedEffect_list  <- list(B_fixed  = B_fixed,
                          beta_fixed  = as.matrix(beta.fixed),
                          noise = "Normal")

#create measurement error list, with one sigma per field

B.e <- kronecker(diag(2),matrix(rep(1, n.obs)))
mError_list <- list(noise = "nsNormal", 
                    B = list(B.e),
                    theta = matrix(log(sigma.e)))


operator_list <- create_operator_matern2Dbivariate(mesh)

operator_list$kappa1 <- kappa1
operator_list$kappa2 <- kappa2
operator_list$tau1   <- tau1
operator_list$tau2   <- tau2
operator_list$rho   <- rho
operator_list$theta   <- theta

processes_list = list(noise = "Normal", V <- list())
processes_list$V[[1]] <- c(operator_list$h,operator_list$h)

cat("Simulate\n")
sim_res <- simulateLongPrior( locs              = list(obs.loc),
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = mError_list,
                              processes_list    = processes_list,
                              operator_list     = operator_list)

n.proc <- length(sim_res$X[[1]])
proj <- inla.mesh.projector(mesh,dims=c(80,80))
par(mfrow=c(1,2))
image.plot(proj$x,proj$y,inla.mesh.project(proj,sim_res$X[[1]][1:(n.proc/2)]),xlab="",ylab="")
image.plot(proj$x,proj$y,inla.mesh.project(proj,sim_res$X[[1]][(n.proc/2+1):n.proc]),xlab="",ylab="")

processes_list$X <- sim_res$X
#operator_list$kappa <- 1
#operator_list$tau   <- 10
operator_list$rho   <- 1.5
#mixedEffect_list$beta_fixed <- 2

if(test.est){
  cat("Estimate\n")
  res.est <- estimateLong(Y                = sim_res$Y,
                          nIter            = nIter,
                          nSim             = 10,
                          locs             = list(obs.loc),
                          mixedEffect_list = mixedEffect_list,
                          measurment_list  = mError_list,
                          processes_list   = processes_list,
                          operator_list    = operator_list,
                          learning_rate = 0.9,
                          nBurnin_learningrate = 50,
                          silent = FALSE)
  
  par(mfrow = c(1,3))
  matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1)
  matplot(t(matrix(rep(matrix(beta.fixed),nIter),length(beta.fixed),nIter)),add=TRUE,lty=2)
  plot(res.est$operator_list$tauVec,type="l",main="process tau")
  lines(rep(tau,nIter),col=2)
  plot(res.est$operator_list$kappaVec,type="l",main="process kappa")
  lines(rep(kappa,nIter),col=2)
  
}


if(test.pred){
  locs.pred <- proj$lattice$loc
  Bfixed.pred <- list(matrix(rep(1, dim(locs.pred)[1])))
  res <- predictLong( Y                = sim_res$Y,
                      locs.pred        = list(locs.pred),
                      Bfixed.pred      = Bfixed.pred,
                      type             = "Smoothing",
                      nSim             = 100,
                      locs             = list(obs.loc),
                      mixedEffect_list = mixedEffect_list,
                      measurment_list  = mError_list,
                      processes_list   = processes_list,
                      operator_list    = operator_list)
  
  par(mfrow=c(1,3))
  image.plot(proj$x,proj$y,inla.mesh.project(proj, sim_res$X[[1]]),xlab="",ylab="")
  image.plot(proj$x,proj$y,matrix(res$W.summary[[1]]$Mean,80,80),xlab="",ylab="")
  image.plot(proj$x,proj$y,matrix(res$W.summary[[1]]$Var,80,80),xlab="",ylab="")
}
