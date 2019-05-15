rm(list = ls())
graphics.off()
library(ngme)
library(INLA)
library(fields)
library(ggplot2)
library(fields)
library(gridExtra)
library(RandomFields)
data(weather)

loc <- weather[,3:4]
pres <- weather[,1]
temp <- weather[,2]
n.obs <- length(pres)

data <- data.frame(weather)

df = data.frame(x = loc[,1],y=loc[,2],z=pres)
p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df = data.frame(x = loc[,1],y=loc[,2],z=temp)
p2 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
grid.arrange(p1,p2,ncol=2)

########################################
# Create mesh for process
########################################
mesh <- inla.mesh.create.helper( loc,
                                 cutoff=0,
                                 max.edge=c(1,1),
                                 offset=c(-0.1,-0.1),
                                 min.angle=20)

#define locations where we want predict the process
proj <- inla.mesh.projector(mesh,dims=c(80,80))
data.pred = data.frame(lon = proj$lattice$loc[,1],
                       lat = proj$lattice$loc[,2])


########################################################
# estimate univariate Gaussian model for preassure
#######################################################

res.est.pres <- ngme.spatial(pres ~ 1,
                             data = data,
                             location.names = c("lon","lat"),
                             silent = FALSE,
                             nIter = 10000,
                             mesh = mesh,
                             controls = list(learning.rate = 0.9,
                                             polyak.rate = 0.1,
                                             nBurnin = 100,
                                             nSim = 4,
                                             step0 = 1,
                                             alpha = 0.3))

cat("beta = ", res.est.pres$fixed_est, "kappa = ", res.est.pres$operator_kappa,
    "tau = ", res.est.pres$operator_tau, "sigma.e = ", res.est.pres$meas_error_sigma)

#compute prediction
res.pred.pres <-predict(res.est.pres, data = data.pred)

df = data.frame(x = loc[,1],y=loc[,2],z=pres)
p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df <- expand.grid(x= proj$x, y = proj$y)
df$z <- res.pred.pres$predictions$X.summary[[1]]$Mean
p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
grid.arrange(p1,p2,ncol=2)


#compute accuracy measure from leave-one-out crossvalidation
res.cv.pres <-predict(res.est.pres, type = "LOOCV")

cat("mae = ", res.cv.pres$median.mae.mean.predictor,"crps =",res.cv.pres$median.crps)

########################################
#same for temp
########################################
  
res.est.temp <- ngme.spatial(temp ~ 1,
                             process = c("Normal"),
                             error = "Normal",
                             data = data,
                             location.names = c("lon","lat"),
                             silent = FALSE,
                             nIter = 10000,
                             mesh = mesh,
                             controls = list(learning.rate = 0.9,
                                             polyak.rate = 0.1,
                                             nBurnin = 100,
                                             nSim = 4,
                                             step0 = 1,
                                             alpha = 0.3))

cat("beta = ", res.est.temp$fixed_est, "kappa = ", res.est.temp$operator_kappa,
    "tau = ", res.est.temp$operator_tau, "sigma.e = ", res.est.temp$meas_error_sigma)

#compute prediction
res.pred.temp <-predict(res.est.temp, data = data.pred)

df = data.frame(x = loc[,1],y=loc[,2],z=temp)
p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df <- expand.grid(x= proj$x, y = proj$y)
df$z <- res.pred.temp$predictions$X.summary[[1]]$Mean
p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
grid.arrange(p1,p2,ncol=2)


#compute accuracy measure from leave-one-out crossvalidation
res.cv.temp <-predict(res.est.pres, type = "LOOCV")

cat("mae = ", res.cv.temp$median.mae.mean.predictor,"crps =",res.cv.temp$median.crps)

####################################################
# Multivariate gaussian
#################################################
  
#Fit gaussian model to pressure data
res.est.gaus <- ngme.spatial(fixed = pres ~ 1,
                             fixed2 = temp ~ 1,
                             process = c("Normal","matern"),
                             error = "Normal",
                             data = data,
                             location.names = c("lon","lat"),
                             silent = FALSE,
                             nIter = 10000,
                             mesh = mesh,
                             controls = list(learning.rate = 0.9,
                                             polyak.rate = 0.1,
                                             nBurnin = 100,
                                             nSim = 4,
                                             step0 = 1,
                                             alpha = 0.3))

cat("beta = ", res.est.gaus$fixed_est, "kappa = ", res.est.gaus$operator_kappa,
    "tau = ", res.est.gaus$operator_tau, "sigma.e = ", res.est.gaus$meas_error_sigma,"\n")

#compute accuracy measure from leave-one-out crossvalidation
res.cv.gaus <-predict(res.est.gaus, type = "LOOCV")

cat("mae = ", res.cv.gaus$median.mae.mean.predictor,"crps =",res.cv.gaus$median.crps)

#define locations where we want predict the process
proj <- inla.mesh.projector(mesh,dims=c(80,80))
data.pred = data.frame(lon = proj$lattice$loc[,1],
                       lat = proj$lattice$loc[,2])

#compute prediction                         
res.pred.gaus <-predict(res.est.gaus, data = data.pred)

df = data.frame(x = loc[,1],y=loc[,2],z=pres)
p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df <- expand.grid(x= proj$x, y = proj$y)
df$z <- res.pred.gaus$predictions$X.summary[[1]]$Mean[,1]
p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
df = data.frame(x = loc[,1],y=loc[,2],z=temp)
p3 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df <- expand.grid(x= proj$x, y = proj$y)
df$z <- res.pred.gaus$predictions$X.summary[[1]]$Mean[,2]
p4 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 

grid.arrange(p1,p2,p3,p4,ncol=2)


########################################################
#univariate NIG for preassure
#######################################################

res.est.nig.pres <- ngme.spatial(pres ~ 1,
                             data = data,
                             location.names = c("lon","lat"),
                             process = c("NIG","matern"),
                             silent = FALSE,
                             nIter = 10000,
                             mesh = mesh,
                             controls = list(learning.rate = 0.9,
                                             polyak.rate = 0.1,
                                             nBurnin = 100,
                                             nSim = 4,
                                             step0 = 1,
                                             alpha = 0.3),
                             init.fit = res.est.pres)

cat("beta = ", res.est.nig.pres$mixedEffect_list$beta_fixed, 
    "kappa = ", res.est.nig.pres$operator_list$kappa,
    "sigma.e = ", res.est.nig.pres$measurementError_list$sigma,
    "nu = ", res.est.nig.pres$processes_list$nu,
    "mu = ", res.est.nig.pres$processes_list$mu)

par(mfrow = c(2,3))
matplot(res.est.nig.pres$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1,xlab="",ylab="")
matplot(res.est.nig.pres$measurementError_list$sigma_vec,type="l",main="noise",col=1,xlab="",ylab="")
plot(res.est.nig.pres$operator_list$tauVec,type="l",main="process tau")
plot(res.est.nig.pres$operator_list$kappaVec,type="l",main="process kappa")
plot(res.est.nig.pres$processes_list$nu_vec,type="l",main="process nu")
plot(res.est.nig.pres$processes_list$mu_vec,type="l",main="process mu")

#compute prediction
res.pred.nig.pres <-predict(res.est.nig.pres, data = data.pred)

df = data.frame(x = loc[,1],y=loc[,2],z=pres)
p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df <- expand.grid(x= proj$x, y = proj$y)
df$z <- res.pred.nig.pres$predictions$X.summary[[1]]$Mean
p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
grid.arrange(p1,p2,ncol=2)


#compute accuracy measure from leave-one-out crossvalidation
res.cv.nig.pres <-predict(res.est.nig.pres, type = "LOOCV")

cat("mae = ", res.cv.nig.pres$median.mae.mean.predictor,"crps =",res.cv.nig.pres$median.crps)
########################################################
#univariate NIG for temperature
#######################################################
  
res.est.nig.temp <- ngme.spatial(temp ~ 1,
                                 data = data,
                                 process = c("NIG","matern"),
                                 location.names = c("lon","lat"),
                                 silent = FALSE,
                                 nIter = 10000,
                                 mesh = mesh,
                                 controls = list(learning.rate = 0.9,
                                                 polyak.rate = 0.1,
                                                 nBurnin = 100,
                                                 nSim = 4,
                                                 step0 = 1,
                                                 alpha = 0.3),
                                 init.fit = res.est.temp)

cat("beta = ", res.est.nig.temp$mixedEffect_list$beta_fixed, 
    "kappa = ", res.est.nig.temp$operator_list$kappa,
    "sigma.e = ", res.est.nig.temp$measurementError_list$sigma,
    "nu = ", res.est.nig.temp$processes_list$nu,
    "mu = ", res.est.nig.temp$processes_list$mu)

par(mfrow = c(2,3))
matplot(res.est.nig.temp$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1,xlab="",ylab="")
matplot(res.est.nig.temp$measurementError_list$sigma_vec,type="l",main="noise",col=1,xlab="",ylab="")
plot(res.est.nig.temp$operator_list$tauVec,type="l",main="process tau")
plot(res.est.nig.temp$operator_list$kappaVec,type="l",main="process kappa")
plot(res.est.nig.temp$processes_list$nu_vec,type="l",main="process nu")
plot(res.est.nig.temp$processes_list$mu_vec,type="l",main="process mu")

#compute prediction
res.pred.nig.temp <-predict(res.est.nig.temp, data = data.pred)

df = data.frame(x = loc[,1],y=loc[,2],z=temp)
p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df <- expand.grid(x= proj$x, y = proj$y)
df$z <- res.pred.nig.temp$predictions$X.summary[[1]]$Mean
p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
grid.arrange(p1,p2,ncol=2)


#compute accuracy measure from leave-one-out crossvalidation
res.cv.nig.temp <-predict(res.est.nig.temp, type = "LOOCV")

####################################################
# Multivariate NIG general
#################################################

#Fit gaussian model to pressure data
res.est.nig <- ngme.spatial(fixed = pres ~ 1,
                            fixed2 = temp ~ 1,
                            process = c("NIG","matern"),
                            error = "Normal",
                            data = data,
                            location.names = c("lon","lat"),
                            silent = FALSE,
                            nIter = 10000,
                            mesh = mesh,
                            controls = list(learning.rate = 0.9,
                                           polyak.rate = 0.1,
                                           nBurnin = 100,
                                           nSim = 4,
                                           step0 = 1,
                                           alpha = 0.3),
                            init.fit = res.est.gaus)

cat("beta = ", res.est.nig$fixed_est, "kappa = ", res.est.nig$operator_kappa,
    "tau = ", res.est.nig$operator_tau, "sigma.e = ", res.est.nig$meas_error_sigma,"\n")


par(mfrow = c(3,3))
matplot(res.est.nig$mixedEffect_list$betaf_vec,type="l",main="fixed effects",xlab="",ylab="")
matplot(res.est.nig$measurementError_list$theta_vec,type="l",main="noise",xlab="",ylab="")
matplot(cbind(res.est.nig$operator_list$tau1Vec,res.est.nig$operator_list$tau2Vec),
        type="l",main="process tau",xlab="",ylab="")
matplot(cbind(res.est.nig$operator_list$kappa1Vec,res.est.nig$operator_list$kappa2Vec),
        type="l",main="process kappa",xlab="",ylab="")
plot(res.est.nig$operator_list$rhoVec,type="l",main="process rho",xlab="",ylab="")
plot(res.est.nig$operator_list$thetaVec,type="l",main="process theta",xlab="",ylab="")  
matplot(res.est.nig$processes_list$nu_vec, type="l",main="process nu",xlab="",ylab="")
matplot(res.est.nig$processes_list$mu_vec, type="l",main="process mu",xlab="",ylab="")

#compute accuracy measure from leave-one-out crossvalidation
res.cv.nig <-predict(res.est.nig, type = "LOOCV")

cat("mae = ", res.cv.nig$median.mae.mean.predictor,"crps =",res.cv.nig$median.crps)

#compute prediction                         
res.pred.nig <-predict(res.est.nig, data = data.pred)

df = data.frame(x = loc[,1],y=loc[,2],z=pres)
p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df <- expand.grid(x= proj$x, y = proj$y)
df$z <- res.pred.gaus$predictions$X.summary[[1]]$Mean[,1]
p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 

df = data.frame(x = loc[,1],y=loc[,2],z=temp)
p3 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df <- expand.grid(x= proj$x, y = proj$y)
df$z <- res.pred.gaus$predictions$X.summary[[1]]$Mean[,2]
p4 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 

grid.arrange(p1,p2,p3,p4,ncol=2)

cat(c(median(abs(res.cv.nig$Y.summary[[1]]$Mean[1:n.obs] - pres)),
      median(res.cv.nig$Y.summary[[1]]$crps[1:n.obs]),
      median(abs(res.cv.nig$Y.summary[[1]]$Mean[(n.obs+1):(2*n.obs)] - temp)),
      median(res.cv.nig$Y.summary[[1]]$crps[(n.obs+1):(2*n.obs)])))
  
  
GGi = c(median(abs(res.cv.pres$Y.summary[[1]]$Mean - pres)),
        median(res.cv.pres$Y.summary[[1]]$crps),
        median(abs(res.cv.temp$Y.summary[[1]]$Mean - temp)),
        median(res.cv.temp$Y.summary[[1]]$crps))
GGl = c(median(abs(res.cv.gauss$Y.summary[[1]]$Mean[1:n.obs] - pres)),
        median(res.cv.gauss$Y.summary[[1]]$crps[1:n.obs]), 
        median(abs(res.cv.gauss$Y.summary[[1]]$Mean[(n.obs+1):(2*n.obs)] - temp)),
        median(res.cv.gauss$Y.summary[[1]]$crps[(n.obs+1):(2*n.obs)]))
NNi = c(median(abs(res.cv.nig.pres$Y.summary[[1]]$Mean[1:n.obs] - pres)),
        median(res.cv.nig.pres$Y.summary[[1]]$crps),
        median(abs(res.cv.nig.temp$Y.summary[[1]]$Mean - temp)),
        median(res.cv.nig.temp$Y.summary[[1]]$crps))
NNg = c(median(abs(res.cv.nig.full$Y.summary[[1]]$Mean[1:n.obs] - pres)),
        median(res.cv.nig.full$Y.summary[[1]]$crps[1:n.obs]),
        median(abs(res.cv.nig.full$Y.summary[[1]]$Mean[(n.obs+1):(2*n.obs)] - temp)),
        median(res.cv.nig.full$Y.summary[[1]]$crps[(n.obs+1):(2*n.obs)]))

results <- data.frame(mae.pres = c(GGi[1],GGl[1],NNi[1],NNg[1]),
                      crps.pres = c(GGi[2],GGl[2],NNi[2],NNg[2]),
                      mae.temp = c(GGi[3],GGl[3],NNi[3],NNg[3]),
                      crps.temp = c(GGi[4],GGl[4],NNi[4],NNg[4]),
                      row.names = c("Gaus indep", "Gauss lower", "NIG indep", "NIG general"))
print(results)
  