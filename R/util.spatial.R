#'
#' @title Collects emperical covariances, and cross covarianes
#'
#' @description Using bin data to emperical estimate covariances
#' @param loc  A list of locations      (d x 2)
#' @param Y    A list of observations   (m x 1-2)
#' @param Bins right locations of the bins {0},(0,Bins[1]], ... (Bins[n-1], Bins[n]] (Bins[n], \inf)
#' @param nBins create bins by equally spaced location
#' @details Standard empirical covariance estimation.
#' @return A list with bins and values
emp.cov <- function(loc, Y, Bins = NULL, nBins = 10){
  N <- length(loc) 
  d <- dim(Y[[1]])[2]
  data <- c()
  for(i in 1:N){
    D <- rdist(loc[[i]],loc[[i]])
    D[upper.tri(D)] <- -1
    if(d == 1){
      y <- Y[[i]] - mean(Y[[i]],  na.rm=T)
      yyt <- y%*%t(y)
      data <- rbind(data, cbind(D[D>=0], yyt[D>=0]))
    }else{
      y1    <- Y[[i]][ ,1] - mean(Y[[i]][ ,1], na.rm=T)
      y2    <- Y[[i]][ ,2] - mean(Y[[i]][ ,2], na.rm=T)
      yyt1  <- y1%*%t(y1)
      yyt2  <- y2%*%t(y2)
      yyt12 <- y1%*%t(y2)
      yyt21 <- t(yyt12)
      data <- rbind(data, cbind(D[D>=0], yyt1[D>=0], yyt2[D>=0], yyt12[D>=0], yyt21[D>=0]))
    }
  }
  if(is.null(Bins)){
    Bins <- quantile(data[data[,1]>0,1], probs = seq(0, 0.1, length.out = nBins+2))
    Bins <- Bins[2:(nBins+1)]
  }
  ns  <- matrix(0, nrow=1, ncol=length(Bins)+1)
  if(d==1){
    res <- matrix(0, nrow=1, ncol=length(Bins)+1)
    res[1] <- mean(data[data[,1]==0, 2], na.rm=T)
  }else{
    res <- matrix(0, nrow=4, ncol=nBins+1)
    res[,1] <- colMeans(data[data[,1]==0, 2:5], na.rm=T)
  }
  Bins <- c(0, Bins, max(D)+1)
  for(i in 2:(nBins+1))
  {
    index <- (Bins[i] < data[,1])* (data[,1] <= Bins[i+1])==1
    if(d==1){
      res[i] <- mean(data[index,2])
    }else{
      res[, i] <- colMeans(data[index,2:5], na.rm=T)
    }
  }
  return(list(Bins = Bins[1:(length(Bins)-1)]), res = res)
}