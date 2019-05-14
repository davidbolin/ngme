rm(list = ls())
graphics.off()
library(ngme)
library(INLA)
library(fields)
library(ggplot2)
library(fields)
library(gridExtra)
library(RandomFields)

#Load data
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


#Build spatial mesh
mesh <- inla.mesh.create.helper( loc,
                                 cutoff=0,
                                 max.edge=c(1,1),
                                 offset=c(-0.1,-0.1),
                                 min.angle=20)


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
s  
  
    

