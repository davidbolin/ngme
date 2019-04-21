rm(list = ls())
graphics.off()
library(ngme)
library(INLA)
library(fields)
library(ggplot2)
library(fields)
library(gridExtra)
#First estimate stationary model:
test.pred = FALSE
test.est = TRUE
nIter = 1000

noise="Gaussian"
kappa1 = 1
kappa2 = 1
tau1 = 5
tau2 = 5
rho = 1
theta = 0
sigma.e = c(0.01,0.1)
beta.fixed = c(0,0)

n.lattice = 40
n.obs=200 #number of observations per replicate
n.rep = 3 #number of replicates

#create mesh 
x=seq(from=0,to=10,length.out=n.lattice)
lattice=inla.mesh.lattice(x=x,y=x)
mesh=inla.mesh.create(lattice=lattice, extend=FALSE, refine=FALSE)

#create observation locations, both fields are observed at the same locations
obs.loc <- list()
for(i in 1:n.rep){
  obs.loc[[i]] = cbind(runif(n.obs)*diff(range(x))+min(x),
                  runif(n.obs)*diff(range(x))+min(x))
  
}


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


operator_list <- create_operator_matern2Dbivariate(mesh)

operator_list$kappa1 <- kappa1
operator_list$kappa2 <- kappa2
operator_list$tau1   <- tau1
operator_list$tau2   <- tau2
operator_list$rho   <- rho
operator_list$theta   <- theta

processes_list = list(noise = "Normal", V <- list())
for(i in 1:n.rep){
  processes_list$V[[i]] <- c(operator_list$h[[1]],operator_list$h[[1]])  
  processes_list$X[[i]] <- 0*processes_list$V[[i]]
}


cat("Simulate\n")
sim_res <- simulateLongPrior( locs              = obs.loc,
                              mixedEffect_list  = mixedEffect_list,
                              measurment_list   = mError_list,
                              processes_list    = processes_list,
                              operator_list     = operator_list)

n.proc <- length(sim_res$X[[1]])
proj <- inla.mesh.projector(mesh,dims=c(80,80))

df <- expand.grid(x= proj$x, y = proj$y)
df$z <- c(inla.mesh.project(proj,sim_res$X[[1]][1:(n.proc/2)]))
df2 <- expand.grid(x= proj$x, y = proj$y)
df2$z <- c(inla.mesh.project(proj,sim_res$X[[1]][(n.proc/2+1):n.proc]))
p1 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
p2 <- ggplot(df2, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=sim_res$Y[[1]][,1])
df2 = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=sim_res$Y[[1]][,2])
p3 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
p4 <- ggplot(df2) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
grid.arrange(p1,p2,p3,p4,ncol=2)

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
                          locs             = obs.loc,
                          mixedEffect_list = mixedEffect_list,
                          measurment_list  = mError_list,
                          processes_list   = processes_list,
                          operator_list    = operator_list,
                          learning_rate = 0.9,
                          nBurnin_learningrate = 50,
                          silent = FALSE)
  
  par(mfrow = c(2,3))
  matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1,xlab="",ylab="")
  abline(beta.fixed[1],0,col=2)
  abline(beta.fixed[2],0,col=2)
  matplot(res.est$measurementError_list$theta_vec,type="l",main="noise",col=1,xlab="",ylab="")
  abline(log(sigma.e[1]),0,col=2)
  abline(log(sigma.e[2]),0,col=2)
  plot(res.est$operator_list$tau1Vec,type="l",main="process tau")
  lines(res.est$operator_list$tau2Vec)
  abline(tau1,0,col=2)
  abline(tau2,0,col=2)
  plot(res.est$operator_list$kappa1Vec,type="l",main="process kappa")
  lines(res.est$operator_list$kappa2Vec)
  abline(kappa1,0,col=2)
  abline(kappa2,0,col=2)
  plot(res.est$operator_list$rhoVec,type="l",main="process rho")
  abline(rho,0,col=2)
  
}


if(test.pred){
  locs.pred <- proj$lattice$loc
  Bf <- kronecker(diag(2),matrix(rep(1, dim(locs.pred)[1])))
  Bfixed.pred  <- list(Bf)
  
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
  
  df = data.frame(x = obs.loc[,1],y=obs.loc[,2],z=sim_res$Y[[1]][,1])
  df2 = data.frame(x = obs.loc[,1],y=obs.loc[,2],z=sim_res$Y.star[[1]][,2])
  df3 = data.frame(x = locs.pred[,1],y=locs.pred[,2],z=res$X.summary[[1]]$Mean[1:length(locs.pred[,1])])
  df4 = data.frame(x = locs.pred[,1],y=locs.pred[,2],z=res$X.summary[[1]]$Mean[(length(locs.pred[,1])+1):(2*length(locs.pred[,1]))])
  p1 <- ggplot(df, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  p2 <- ggplot(df2, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  p3 <- ggplot(df3,aes(x,y))+geom_raster(aes(fill=z))+ scale_fill_gradientn(colours=tim.colors(100)) 
  p4 <- ggplot(df4,aes(x,y))+geom_raster(aes(fill=z))+ scale_fill_gradientn(colours=tim.colors(100)) 
  grid.arrange(p1,p2,p3,p4,ncol=2)
}
