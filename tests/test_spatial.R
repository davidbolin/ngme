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
nIter = 100

noise="Normal"
kappa = 1
tau = 0.5

sigma.e = 0.01
beta.fixed = c(0)

n.lattice = 20
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
Bf <- matrix(rep(1, n.obs))
B_fixed <- list()
for(i in 1:n.rep){
  B_fixed[[i]]  <- Bf
}
mixedEffect_list  <- list(B_fixed  = B_fixed,
                          beta_fixed  = as.matrix(beta.fixed),
                          noise = "Normal")

mError_list <- list(noise = "Normal", 
                    theta = matrix(sigma.e))

operator_list <- create_operator_matern2D(mesh)

operator_list$kappa <- kappa
operator_list$tau   <- tau

processes_list = list(noise = noise, V <- list())
for(i in 1:n.rep){
  processes_list$V[[i]] <- c(operator_list$h[[1]])  
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
df$z <- c(inla.mesh.project(proj,sim_res$X[[1]]))
p1 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=sim_res$Y[[1]])
p2 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
grid.arrange(p1,p2,ncol=2)



processes_list$X <- sim_res$X

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
  
  par(mfrow = c(2,2))
  matplot(res.est$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1,xlab="",ylab="")
  abline(beta.fixed[1],0,col=2)
  matplot(res.est$measurementError_list$sigma_vec,type="l",main="noise",col=1,xlab="",ylab="")
  abline(sigma.e,0,col=2)
  plot(res.est$operator_list$tauVec,type="l",main="process tau")
  abline(tau,0,col=2)
  plot(res.est$operator_list$kappaVec,type="l",main="process kappa")
  abline(kappa,0,col=2)
}


if(test.pred){
  locs.pred <- Bfixed.pred <- list()
  for(i in 1:n.rep){
    locs.pred[[i]] <- proj$lattice$loc
    Bfixed.pred[[i]] <- matrix(rep(1, dim(locs.pred[[i]])[1]))
  }
  
  
  
  res <- predictLong( Y                = sim_res$Y,
                      locs.pred        = locs.pred,
                      Bfixed.pred      = Bfixed.pred,
                      type             = "Smoothing",
                      nSim             = 100,
                      locs             = obs.loc,
                      mixedEffect_list = mixedEffect_list,
                      measurment_list  = mError_list,
                      processes_list   = processes_list,
                      operator_list    = operator_list)
  
  df <- expand.grid(x= proj$x, y = proj$y)
  df$z <- c(inla.mesh.project(proj,sim_res$X[[1]]))
  p1 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=sim_res$Y[[1]])
  df2 = data.frame(x = locs.pred[[1]][,1],y=locs.pred[[1]][,2],z=res$X.summary[[1]]$Mean)
  p2 <- ggplot(df, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  p3 <- ggplot(df2,aes(x,y))+geom_raster(aes(fill=z))+ scale_fill_gradientn(colours=tim.colors(100)) 
  grid.arrange(p1,p2,p3,ncol=2)
}

