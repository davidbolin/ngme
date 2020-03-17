rm(list = ls())
graphics.off()
library(ngme)
library(INLA)
library(fields)
library(ggplot2)
library(fields)
library(gridExtra)

#First estimate stationary model:
test.pred = TRUE
test.est = FALSE
test.cv = FALSE
test.missing = TRUE
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
  
  Y <- sim_res$Y
  if(test.missing){
    Y[[1]][1:10] = NA
  }
  
  res <- predictLong( Y                = Y,
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

if(test.cv){
  #define indices for where to predict and where to observe in each fold
  pred.ind <- list()
  obs.ind <- list()
  for(i in 1:n.obs){
    pred.ind[[i]] <- i
    obs.ind[[i]] <- setdiff(1:n.obs,i)
  }
  
  #setup cluster
  n.cores = 8
  cl <- makeCluster(n.cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(pred.ind), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  parallel::clusterExport(cl, varlist = c('pred.ind','obs.loc', 'mixedEffect_list','sim_res','obs.ind',
                                          'mError_list','processes_list','operator_list'), envir = environment())
  
  #run CV
  preds.list <- foreach(j = 1:length(pred.ind), .options.snow = opts) %dopar%
  {
    locs.pred <- Bfixed.pred <- list()
    Y.i <- Y.val <- loc.i <- list()
    for(i in 1:n.rep){
      locs.pred[[i]] <- obs.loc[[i]][pred.ind[[j]],,drop=FALSE]
      Bfixed.pred[[i]] <- mixedEffect_list$B_fixed[[i]][pred.ind[[j]],,drop=FALSE]
      Y.i[[i]] <- sim_res$Y[[i]][obs.ind[[j]]]
      Y.val[[i]] <- sim_res$Y[[i]][pred.ind[[j]]]
      loc.i[[i]] <- obs.loc[[i]][obs.ind[[j]],,drop=FALSE]
    }
    res <- ngme::predictLong( Y          = Y.i,
                        locs.pred        = locs.pred,
                        Bfixed.pred      = Bfixed.pred,
                        type             = "Smoothing",
                        nSim             = 1000,
                        locs             = loc.i,
                        mixedEffect_list = mixedEffect_list,
                        measurment_list  = mError_list,
                        processes_list   = processes_list,
                        operator_list    = operator_list,
                        crps             = TRUE,
                        Y.val            = Y.val,
                        silent           = TRUE)  
    Y.pred <- Y.crps <- Y.mse <- list()
    for(i in 1:n.rep){
      Y.pred[[i]] = res$Y.summary[[i]]$Mean 
      Y.crps[[i]] = res$Y.summary[[i]]$crps
      Y.mse[[i]] = sum((res$Y.summary[[i]]$Mean - sim_res$Y[[i]][pred.ind[[j]]])^2)
    }
    return(list(pred = Y.pred,
                crps = Y.crps,
                mse = Y.mse))
  }
  close(pb)
  stopCluster(cl)
  #collect results
  Y.pred <- Y.crps <- Y.mse <- list()
  for(i in 1:n.rep){
    Y.pred[[i]] <- 0*sim_res$Y[[i]]
    Y.crps[[i]]<- 0*sim_res$Y[[i]]
    Y.mse[[i]]<- 0*sim_res$Y[[i]]
  }
  for(j in 1:length(pred.ind)){
    for(i in 1:n.rep){
      Y.pred[[i]][pred.ind[[j]]] <- preds.list[[j]]$pred[[i]]
      Y.crps[[i]][pred.ind[[j]]] <- preds.list[[j]]$crps[[i]]
      Y.mse[[i]][pred.ind[[j]]] <- preds.list[[j]]$mse[[i]]
    }  
  }
  
  #plot somre results
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=Y.pred[[1]]-sim_res$Y[[1]])
  p1 <- ggplot(df, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=Y.crps[[1]])
  p2 <- ggplot(df, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  grid.arrange(p1,p2,ncol=2)
  
  #compare with direct CV using predictLong
  res2 <- predictLong(Y                = sim_res$Y,
                      type             = "LOOCV",
                      nSim             = 1000,
                      locs             = obs.loc,
                      mixedEffect_list = mixedEffect_list,
                      measurment_list  = mError_list,
                      processes_list   = processes_list,
                      operator_list    = operator_list,
                      crps             = TRUE)
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=res2$Y.summary[[1]]$Mean-sim_res$Y[[1]])
  p3 <- ggplot(df, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=res2$Y.summary[[1]]$crps)
  p4 <- ggplot(df, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  grid.arrange(p1,p2,p3,p4,ncol=2)
}
