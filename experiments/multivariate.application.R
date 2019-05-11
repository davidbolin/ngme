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

df = data.frame(x = loc[,1],y=loc[,2],z=pres)
p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
df = data.frame(x = loc[,1],y=loc[,2],z=temp)
p2 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
grid.arrange(p1,p2,ncol=2)

########################################################
#estimate independent normal to each dimension
#######################################################
mesh <- inla.mesh.create.helper( loc,
                               cutoff=0,
                               max.edge=c(1,1),
                               offset=c(-0.1,-0.1),
                               min.angle=20)

noise="Normal"
#create observation locations, both fields are observed at the same locations
obs.loc <- list(loc)
B_fixed <- list(matrix(rep(1, n.obs)))
mixedEffect_list  <- list(B_fixed  = B_fixed,
                          beta_fixed  = as.matrix(0),
                          noise = "Normal")

mError_list <- list(noise = "Normal", 
                    theta = matrix(10))

operator_list <- create_operator_matern2D(mesh)

operator_list$kappa <- 1
operator_list$tau   <- 0.1

processes_list = list(noise = noise, V <- list())
processes_list$V[[1]] <- c(operator_list$h[[1]])  
processes_list$X[[1]] <- 0*processes_list$V[[1]]


res.est.pres <- estimateLong(Y                = list(pres),
                          nIter            = 10000,
                          nSim             = 5,
                          locs             = obs.loc,
                          mixedEffect_list = mixedEffect_list,
                          measurment_list  = mError_list,
                          processes_list   = processes_list,
                          operator_list    = operator_list,
                          learning_rate = 0.9,
                          step0 = 1,
                          alpha = 0.3,
                          nBurnin_learningrate = 50,
                          silent = FALSE)
  
  cat("beta = ", res.est.pres$mixedEffect_list$beta_fixed, "kappa = ", res.est.pres$operator_list$kappa,
      "sigma.e = ", res.est.pres$measurementError_list$sigma)
  
  #kappa = res.est$operator_list$kappa
  #tau = res.est$operator_list$tau
  #G = operator_list$G[[1]]
  #C = operator_list$C[[1]]
  ###
  
  #K = 0.5*(tau/kappa^1.5)*(G + kappa^2*C)
  #Q = K%*%solve(C,K)
  #Sigma = solve(Q)
  #A = inla.spde.make.A(mesh,loc=obs.loc[[1]])
  #df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=as.double(sqrt(A%*%matrix(diag(Sigma)))))
  #ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
  
  ##3
  par(mfrow = c(2,2))
  matplot(res.est.pres$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1,xlab="",ylab="")
  matplot(res.est.pres$measurementError_list$sigma_vec,type="l",main="noise",col=1,xlab="",ylab="")
  plot(res.est.pres$operator_list$tauVec,type="l",main="process tau")
  plot(res.est.pres$operator_list$kappaVec,type="l",main="process kappa")


  proj <- inla.mesh.projector(mesh,dims=c(80,80))
  
  locs.pred <- list(proj$lattice$loc)
  Bfixed.pred <- list(matrix(rep(1, dim(locs.pred[[1]])[1])))
  
  
  res.pred.pres <- predictLong( Y                = list(pres),
                      locs.pred        = locs.pred,
                      Bfixed.pred      = Bfixed.pred,
                      type             = "Smoothing",
                      nSim             = 100,
                      locs             = obs.loc,
                      mixedEffect_list = res.est.pres$mixedEffect_list,
                      measurment_list  = res.est.pres$measurementError_list,
                      processes_list   = res.est.pres$processes_list,
                      operator_list    = res.est.pres$operator_list)
  
  
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=pres)
  p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
  df <- expand.grid(x= proj$x, y = proj$y)
  df$z <- res.pred.pres$X.summary[[1]]$Mean
  p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
  grid.arrange(p1,p2,ncol=2)

  #compare with direct CV using predictLong
  res.cv.pres <- predictLong(Y                = list(pres),
                      type             = "LOOCV",
                      nSim             = 100,
                      locs             = obs.loc,
                      mixedEffect_list = res.est.pres$mixedEffect_list,
                      measurment_list  = res.est.pres$measurementError_list,
                      processes_list   = res.est.pres$processes_list,
                      operator_list    = res.est.pres$operator_list,
                      crps             = TRUE,
                      return.samples = TRUE)
  
  cat(c(median(abs(res.cv.pres$Y.summary[[1]]$Mean - pres)),median(res.cv.pres$Y.summary[[1]]$crps)))
  
  ########################################
  #same for temp
  ########################################
  
  res.est.temp <- estimateLong(Y          = list(temp),
                          nIter            = 10000,
                          nSim             = 2,
                          locs             = obs.loc,
                          mixedEffect_list = mixedEffect_list,
                          measurment_list  = mError_list,
                          processes_list   = processes_list,
                          operator_list    = operator_list,
                          learning_rate = 0.9,
                          step0 = 1,
                          alpha = 0.3,
                          nBurnin_learningrate = 50,
                          silent = FALSE)
  
  par(mfrow = c(2,2))
  matplot(res.est.temp$mixedEffect_list$betaf_vec,type="l",main="fixed effects",col=1,xlab="",ylab="")
  matplot(res.est.temp$measurementError_list$sigma_vec,type="l",main="noise",col=1,xlab="",ylab="")
  plot(res.est.temp$operator_list$tauVec,type="l",main="process tau")
  plot(res.est.temp$operator_list$kappaVec,type="l",main="process kappa")
  
  res.pred.temp <- predictLong( Y                = list(temp),
                      locs.pred        = locs.pred,
                      Bfixed.pred      = Bfixed.pred,
                      type             = "Smoothing",
                      nSim             = 100,
                      locs             = obs.loc,
                      mixedEffect_list = res.est.temp$mixedEffect_list,
                      measurment_list  = res.est.temp$measurementError_list,
                      processes_list   = res.est.temp$processes_list,
                      operator_list    = res.est.temp$operator_list)
  
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=pres)
  p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
  df <- expand.grid(x= proj$x, y = proj$y)
  df$z <- res.pred.temp$X.summary[[1]]$Mean
  p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
  grid.arrange(p1,p2,ncol=2)
  
  #compare with direct CV using predictLong
  res.cv.temp <- predictLong(Y                = list(temp),
                          type             = "LOOCV",
                          nSim             = 1000,
                          locs             = obs.loc,
                          mixedEffect_list = res.est.temp$mixedEffect_list,
                          measurment_list  = res.est.temp$measurementError_list,
                          processes_list   = res.est.temp$processes_list,
                          operator_list    = res.est.temp$operator_list,
                          crps             = TRUE,
                          return.samples = TRUE)
  
  cat(c(median(abs(res.cv.temp$Y.summary[[1]]$Mean - temp)),median(res.cv.temp$Y.summary[[1]]$crps)))
  
  
  
  
  ####################################################
  # Multivariate gaussian
  #################################################
  
  
  #create fixed effects-list with one intercept per field
  B_fixed <- list(kronecker(diag(2),matrix(rep(1, n.obs))))
  mixedEffect_list  <- list(B_fixed  = B_fixed,
                            beta_fixed  = as.matrix(c(res.est.pres$mixedEffect_list$beta_fixed,
                                                      res.est.temp$mixedEffect_list$beta_fixed)),
                            noise = "Normal")
  
  #create measurement error list, with one sigma per field
  mError_list <- list(noise = "nsNormal", 
                      B = list(kronecker(diag(2),matrix(rep(1, n.obs)))),
                      theta = matrix(c(log(res.est.pres$measurementError_list$sigma),
                                       log(res.est.temp$measurementError_list$sigma))))
  
  operator_list <- create_operator_matern2Dbivariate(mesh)
  
  operator_list$kappa1 <- res.est.pres$operator_list$kappa
  operator_list$kappa2 <- res.est.temp$operator_list$kappa
  operator_list$tau1   <- res.est.pres$operator_list$tau
  operator_list$tau2   <- res.est.temp$operator_list$tau
  operator_list$rho    <- 0
  operator_list$theta  <- 0
  
  n.grid <- length(operator_list$h[[1]])/2
  processes_list = list(noise = noise, 
                        V = list(c(operator_list$h[[1]])), 
                        X = list(0*c(operator_list$h[[1]])))
  
  res.est.gauss <- estimateLong(Y                = list(cbind(pres,temp)),
                            nIter            = 1000,
                            nSim             = 5,
                            locs             = obs.loc,
                            mixedEffect_list = mixedEffect_list,
                            measurment_list  = mError_list,
                            processes_list   = processes_list,
                            operator_list    = operator_list,
                            learning_rate = 0.9,
                            step0 = 1,
                            alpha = 0.3,
                            nBurnin_learningrate = 50,
                            silent = FALSE)
    
  
  par(mfrow = c(2,3))
  matplot(res.est.gauss$mixedEffect_list$betaf_vec,type="l",main="fixed effects",xlab="",ylab="")
  matplot(res.est.gauss$measurementError_list$theta_vec,type="l",main="noise",xlab="",ylab="")
  matplot(cbind(res.est.gauss$operator_list$tau1Vec,res.est.gauss$operator_list$tau2Vec),
          type="l",main="process tau",xlab="",ylab="")
  matplot(cbind(res.est.gauss$operator_list$kappa1Vec,res.est.gauss$operator_list$kappa2Vec),
          type="l",main="process kappa",xlab="",ylab="")
  plot(res.est.gauss$operator_list$rhoVec,type="l",main="process rho",xlab="",ylab="")
  
  
  Bfixed.pred <- list(kronecker(diag(2),matrix(rep(1, dim(locs.pred[[1]])[1]))))
  Be.pred <- list(kronecker(diag(2),matrix(rep(1, dim(locs.pred[[1]])[1]))))
  res.est.gauss$measurementError_list$Bpred <- Be.pred
  
  res.pred.gauss <- predictLong( Y                = list(cbind(pres,temp)),
                      locs.pred        = locs.pred,
                      Bfixed.pred      = Bfixed.pred,
                      type             = "Smoothing",
                      nSim             = 100,
                      locs             = obs.loc,
                      mixedEffect_list = res.est.gauss$mixedEffect_list,
                      measurment_list  = res.est.gauss$measurementError_list,
                      processes_list   = res.est.gauss$processes_list,
                      operator_list    = res.est.gauss$operator_list)
  
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=pres)
  df2 = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=temp)
  df3 = data.frame(x = locs.pred[[1]][,1],y=locs.pred[[1]][,2],z=res.pred.gauss$X.summary[[1]]$Mean[1:length(locs.pred[[1]][,1])])
  df4 = data.frame(x = locs.pred[[1]][,1],y=locs.pred[[1]][,2],z=res.pred.gauss$X.summary[[1]]$Mean[(length(locs.pred[[1]][,1])+1):(2*length(locs.pred[[1]][,1]))])
  p1 <- ggplot(df, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  p2 <- ggplot(df2, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  p3 <- ggplot(df3,aes(x,y))+geom_raster(aes(fill=z))+ scale_fill_gradientn(colours=tim.colors(100)) 
  p4 <- ggplot(df4,aes(x,y))+geom_raster(aes(fill=z))+ scale_fill_gradientn(colours=tim.colors(100)) 
  grid.arrange(p3,p4,p1,p2,ncol=2)
  
  res.cv.gauss <- predictLong(Y                = list(cbind(pres,temp)),
                              type             = "LOOCV",
                              nSim             = 100,
                              locs             = obs.loc,
                              mixedEffect_list = res.est.gauss$mixedEffect_list,
                              measurment_list  = res.est.gauss$measurementError_list,
                              processes_list   = res.est.gauss$processes_list,
                              operator_list    = res.est.gauss$operator_list,
                              crps = TRUE)
  
  cat(c(median(abs(res.cv.gauss$Y.summary[[1]]$Mean[1:n.obs] - pres)),
        median(res.cv.gauss$Y.summary[[1]]$crps[1:n.obs]),
        median(abs(res.cv.gauss$Y.summary[[1]]$Mean[(n.obs+1):(2*n.obs)] - temp)),
        median(res.cv.gauss$Y.summary[[1]]$crps[(n.obs+1):(2*n.obs)])))
  
  
  
  ########################################################
  #univariate NIG for preassure
  #######################################################
  
  
  obs.loc <- list(loc)
  B_fixed <- list(matrix(rep(1, n.obs)))
  mixedEffect_list  <- list(B_fixed  = B_fixed,
                            beta_fixed  = as.matrix(0),
                            noise = "Normal")
  
  mError_list <- list(noise = "Normal", 
                      theta = matrix(10))
  
  operator_list <- create_operator_matern2D(mesh)
  
  operator_list$kappa <- 1
  operator_list$tau   <- 0.1
  
  processes_list = list(noise = "NIG", V <- list())
  processes_list$V[[1]] <- c(operator_list$h[[1]])  
  processes_list$X[[1]] <- 0*processes_list$V[[1]]
  
  
  res.est.nig.pres <- estimateLong(Y                = list(pres),
                               nIter            = 10000,
                               nSim             = 5,
                               locs             = obs.loc,
                               mixedEffect_list = mixedEffect_list,
                               measurment_list  = mError_list,
                               processes_list   = processes_list,
                               operator_list    = operator_list,
                               learning_rate = 0.9,
                               step0 = 1,
                               alpha = 0.3,
                               nBurnin_learningrate = 50,
                               silent = FALSE)
  
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
  
  Bfixed.pred <- list(matrix(rep(1, dim(locs.pred[[1]])[1])))
  
  
  res.pred.nig.pres <- predictLong( Y                = list(pres),
                                locs.pred        = locs.pred,
                                Bfixed.pred      = Bfixed.pred,
                                type             = "Smoothing",
                                nSim             = 100,
                                locs             = obs.loc,
                                mixedEffect_list = res.est.nig.pres$mixedEffect_list,
                                measurment_list  = res.est.nig.pres$measurementError_list,
                                processes_list   = res.est.nig.pres$processes_list,
                                operator_list    = res.est.nig.pres$operator_list)
  
  
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=pres)
  p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
  df <- expand.grid(x= proj$x, y = proj$y)
  df$z <- res.pred.nig.pres$X.summary[[1]]$Mean
  p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
  grid.arrange(p1,p2,ncol=2)
  
  #compare with direct CV using predictLong
  res.cv.nig.pres <- predictLong(Y                = list(pres),
                             type             = "LOOCV",
                             nSim             = 100,
                             locs             = obs.loc,
                             mixedEffect_list = res.est.nig.pres$mixedEffect_list,
                             measurment_list  = res.est.nig.pres$measurementError_list,
                             processes_list   = res.est.nig.pres$processes_list,
                             operator_list    = res.est.nig.pres$operator_list,
                             crps             = TRUE)
  
  cat(c(median(abs(res.cv.nig.pres$Y.summary[[1]]$Mean - pres)),median(res.cv.nig.pres$Y.summary[[1]]$crps)))
  
  ########################################################
  #univariate NIG for temperature
  #######################################################
  
  B_fixed <- list(matrix(rep(1, n.obs)))
  mixedEffect_list  <- list(B_fixed  = B_fixed,
                            beta_fixed  = as.matrix(0),
                            noise = "Normal")
  
  mError_list <- list(noise = "Normal", 
                      theta = matrix(10))
  
  operator_list <- create_operator_matern2D(mesh)
  
  operator_list$kappa <- 1
  operator_list$tau   <- 0.1
  
  processes_list = list(noise = "NIG", V <- list())
  processes_list$V[[1]] <- c(operator_list$h[[1]])  
  processes_list$X[[1]] <- 0*processes_list$V[[1]]
  
  
  res.est.nig.temp <- estimateLong(Y                = list(temp),
                                   nIter            = 10000,
                                   nSim             = 5,
                                   locs             = obs.loc,
                                   mixedEffect_list = mixedEffect_list,
                                   measurment_list  = mError_list,
                                   processes_list   = processes_list,
                                   operator_list    = operator_list,
                                   learning_rate = 0.9,
                                   step0 = 1,
                                   alpha = 0.3,
                                   nBurnin_learningrate = 50,
                                   silent = FALSE)
  
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
  
  Bfixed.pred <- list(matrix(rep(1, dim(locs.pred[[1]])[1])))
  
  
  res.pred.nig.temp <- predictLong( Y                = list(temp),
                                    locs.pred        = locs.pred,
                                    Bfixed.pred      = Bfixed.pred,
                                    type             = "Smoothing",
                                    nSim             = 100,
                                    locs             = obs.loc,
                                    mixedEffect_list = res.est.nig.temp$mixedEffect_list,
                                    measurment_list  = res.est.nig.temp$measurementError_list,
                                    processes_list   = res.est.nig.temp$processes_list,
                                    operator_list    = res.est.nig.temp$operator_list)
  
  
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=pres)
  p1 <- ggplot(df) + geom_point(aes(x,y,colour=z), size=1, alpha=1) + scale_colour_gradientn(colours=tim.colors(100)) 
  df <- expand.grid(x= proj$x, y = proj$y)
  df$z <- res.pred.nig.temp$X.summary[[1]]$Mean
  p2 <- ggplot(df, aes(x, y, fill = z)) + geom_raster() + scale_fill_gradientn(colours=tim.colors(100)) 
  grid.arrange(p1,p2,ncol=2)
  
  #compare with direct CV using predictLong
  res.cv.nig.temp <- predictLong(Y                = list(temp),
                                 type             = "LOOCV",
                                 nSim             = 100,
                                 locs             = obs.loc,
                                 mixedEffect_list = res.est.nig.temp$mixedEffect_list,
                                 measurment_list  = res.est.nig.temp$measurementError_list,
                                 processes_list   = res.est.nig.temp$processes_list,
                                 operator_list    = res.est.nig.temp$operator_list,
                                 crps             = TRUE)
  
  cat(c(median(abs(res.cv.nig.temp$Y.summary[[1]]$Mean - temp)),median(res.cv.nig.temp$Y.summary[[1]]$crps)))
  
  
  ####################################################
  # Multivariate NIG general
  #################################################
  
  
  #create fixed effects-list with one intercept per field
  B_fixed <- list(kronecker(diag(2),matrix(rep(1, n.obs))))
  mixedEffect_list  <- list(B_fixed  = B_fixed,
                            beta_fixed  = as.matrix(c(res.est.pres$mixedEffect_list$beta_fixed,
                                                      res.est.temp$mixedEffect_list$beta_fixed)),
                            noise = "Normal")
  
  #create measurement error list, with one sigma per field
  mError_list <- list(noise = "nsNormal", 
                      B = list(kronecker(diag(2),matrix(rep(1, n.obs)))),
                      theta = matrix(c(log(res.est.pres$measurementError_list$sigma),
                                       log(res.est.temp$measurementError_list$sigma))))
  
  operator_list <- create_operator_matern2Dbivariate(mesh)
  
  operator_list$kappa1 <- res.est.pres$operator_list$kappa
  operator_list$kappa2 <- res.est.temp$operator_list$kappa
  operator_list$tau1   <- res.est.pres$operator_list$tau
  operator_list$tau2   <- res.est.temp$operator_list$tau
  operator_list$rho    <- 0
  operator_list$theta  <- 0
  
  n.grid <- length(operator_list$h[[1]])/2
  
  processes_list = list(noise = "MultiGH", 
                        V = list(c(operator_list$h[[1]])), 
                        X = list(0*c(operator_list$h[[1]])),
                        Bmu = list(kronecker(diag(2),matrix(rep(1, n.grid)))), 
                        Bnu = list(kronecker(diag(2),matrix(rep(1, n.grid)))),
                        mu = as.matrix(c(0,0)), 
                        nu = as.matrix(c(1,1)))
  
  
  res.est.nig.full <- estimateLong(Y                = list(cbind(pres,temp)),
                                nIter            = 1000,
                                nSim             = 5,
                                locs             = obs.loc,
                                mixedEffect_list = mixedEffect_list,
                                measurment_list  = mError_list,
                                processes_list   = processes_list,
                                operator_list    = operator_list,
                                learning_rate = 0.9,
                                step0 = 1,
                                alpha = 0.3,
                                nBurnin_learningrate = 50,
                                silent = FALSE)
  
  
  par(mfrow = c(3,3))
  matplot(res.est.nig.full$mixedEffect_list$betaf_vec,type="l",main="fixed effects",xlab="",ylab="")
  matplot(res.est.nig.full$measurementError_list$theta_vec,type="l",main="noise",xlab="",ylab="")
  matplot(cbind(res.est.nig.full$operator_list$tau1Vec,res.est.nig.full$operator_list$tau2Vec),
          type="l",main="process tau",xlab="",ylab="")
  matplot(cbind(res.est.nig.full$operator_list$kappa1Vec,res.est.nig.full$operator_list$kappa2Vec),
          type="l",main="process kappa",xlab="",ylab="")
  plot(res.est.nig.full$operator_list$rhoVec,type="l",main="process rho",xlab="",ylab="")
  plot(res.est.nig.full$operator_list$thetaVec,type="l",main="process theta",xlab="",ylab="")  
  matplot(res.est.nig.full$processes_list$nu_vec, type="l",main="process nu",xlab="",ylab="")
  matplot(res.est.nig.full$processes_list$mu_vec, type="l",main="process mu",xlab="",ylab="")
  
  Bfixed.pred <- list(kronecker(diag(2),matrix(rep(1, dim(locs.pred[[1]])[1]))))
  Be.pred <- list(kronecker(diag(2),matrix(rep(1, dim(locs.pred[[1]])[1]))))
  res.est.nig.full$measurementError_list$Bpred <- Be.pred
  
  res.pred.nig.full <- predictLong( Y                = list(cbind(pres,temp)),
                                 locs.pred        = locs.pred,
                                 Bfixed.pred      = Bfixed.pred,
                                 type             = "Smoothing",
                                 nSim             = 100,
                                 locs             = obs.loc,
                                 mixedEffect_list = res.est.nig.full$mixedEffect_list,
                                 measurment_list  = res.est.nig.full$measurementError_list,
                                 processes_list   = res.est.nig.full$processes_list,
                                 operator_list    = res.est.nig.full$operator_list)
  
  df = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=pres)
  df2 = data.frame(x = obs.loc[[1]][,1],y=obs.loc[[1]][,2],z=temp)
  df3 = data.frame(x = locs.pred[[1]][,1],y=locs.pred[[1]][,2],z=res.pred.nig.full$X.summary[[1]]$Mean[1:length(locs.pred[[1]][,1])])
  df4 = data.frame(x = locs.pred[[1]][,1],y=locs.pred[[1]][,2],z=res.pred.nig.full$X.summary[[1]]$Mean[(length(locs.pred[[1]][,1])+1):(2*length(locs.pred[[1]][,1]))])
  p1 <- ggplot(df, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  p2 <- ggplot(df2, aes(x, y,color=z)) + geom_point() + scale_color_gradientn(colours=tim.colors(100)) 
  p3 <- ggplot(df3,aes(x,y))+geom_raster(aes(fill=z))+ scale_fill_gradientn(colours=tim.colors(100)) 
  p4 <- ggplot(df4,aes(x,y))+geom_raster(aes(fill=z))+ scale_fill_gradientn(colours=tim.colors(100)) 
  grid.arrange(p3,p4,p1,p2,ncol=2)
  
  res.cv.nig.full <- predictLong(Y                = list(cbind(pres,temp)),
                              type             = "LOOCV",
                              nSim             = 100,
                              locs             = obs.loc,
                              mixedEffect_list = res.est.nig.full$mixedEffect_list,
                              measurment_list  = res.est.nig.full$measurementError_list,
                              processes_list   = res.est.nig.full$processes_list,
                              operator_list    = res.est.nig.full$operator_list,
                              crps = TRUE)
  
  cat(c(median(abs(res.cv.nig.full$Y.summary[[1]]$Mean[1:n.obs] - pres)),
        median(res.cv.nig.full$Y.summary[[1]]$crps[1:n.obs]),
        median(abs(res.cv.nig.full$Y.summary[[1]]$Mean[(n.obs+1):(2*n.obs)] - temp)),
        median(res.cv.nig.full$Y.summary[[1]]$crps[(n.obs+1):(2*n.obs)])))
  
  
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
  