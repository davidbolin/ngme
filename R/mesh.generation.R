#' @title Create meshes.
#'
#' @description A function to create mesh for FEM discretization.
#' @inheritParams create_matrices_FD2
#'
#' @return Returns a list of meshes.
#'
#' @details This is a supplementary function to be used internally by other functions.
#'
#' @seealso \code{\link{create_matrices_Matern}}, \code{\link{create_matrices_FD2}}
#'
#' @examples
#'   \dontrun{
#'   create_meshes.1d(...)
#'   }
#'

create.meshes.1d <- function(locs, n, common.grid, extend = NULL)
{
  loc <- h <- n.list <- list()
  if(missing(extend) | is.null(extend)){
    extend  = c(0,0)
  } else if(length(extend) == 1) {
    extend = rep(extend,2)
  }
  if(common.grid || length(locs)==1){
    min_l <- min(locs[[1]])
    max_l <- max(locs[[1]])
    if(length(locs) > 1){
      for(i in 2:length(locs))
      {
        min_l <- min(min_l, min(locs[[i]]))
        max_l <- max(max_l, max(locs[[i]]))
      }
    }
    loc_len = max_l - min_l
    if(loc_len == 0){
      stop("min(locs) = max(locs)")
    }
    loc[[1]] <- seq(min_l - extend[1]*loc_len, max_l + extend[2]*loc_len, length.out = n)

    hi = (c(diff(loc[[1]]),Inf)+c(Inf,diff(loc[[1]])))/2
    hi[1] <- (loc[[1]][2]-loc[[1]][1])/2
    hi[n] <- (loc[[1]][n]-loc[[1]][n-1])/2

    h[[1]] <- hi
    n.list[[1]] <- n
  } else {
    if(length(n) == 1){
      n = rep(n,length(locs))
    }
    loc <- list()
    for(i in 1:length(locs)){
      min_l <- min(locs[[i]])
      max_l <- max(locs[[i]])
      loc_len = max_l - min_l
      loc[[i]] <- seq(min_l - extend[1]*loc_len, max_l + extend[2]*loc_len, length.out = n[i])

      hi = (c(diff(loc[[i]]),Inf)+c(Inf,diff(loc[[i]])))/2
      hi[1] <- (loc[[i]][2]-loc[[i]][1])/2
      hi[n] <- (loc[[i]][n]-loc[[i]][n-1])/2
      h[[i]] <- hi
      n.list[[i]] <- n
    }
  }
  return(list(loc = loc,h = h,n = n.list))
}


#' @title Create 1d mesh from observation locations.
#'
#' @description A function to create mesh for 1d FEM discretization.
#' @param x A list of measurement locations.
#' @param max.dist The largest distance between nodes in the mesh
#' @param cutoff Merge nodes in x that are closer than cutoff
#'
#' @return Returns a list with mesh nodes, distances, and number of nodes.
#'
#' @details This is a supplementary function to be used internally by other functions.
#'
#' @seealso \code{\link{create_matrices_Matern}}, \code{\link{create_matrices_FD2}}
#'
#' @examples
#'   \dontrun{
#'   x = c(0,0.05,0.2,0.25,0.5,0.7,0.8,0.95,0.99,1)
#'   x <- 2*sort(runif(n))
#'   m <- generate.1d.mesh(x,max.dist = 0.2,cutoff = 0.1)
#'  plot(x,0*x)
#'  points(m$s,0*m$s,pch=4,col=2)
#'   }
#'
generate.1d.mesh <- function(x, 
                             max.dist, 
                             cutoff = 1e-10, 
                             extend=NULL,
                             Y = NULL,
                             loc.Y = NULL,
                             max.dY = -1,
                             nJump  = 3){
  if(missing(x))
    stop('Must supply x')
  refine = TRUE
  if(missing(max.dist) | is.null(max.dist)){
    refine= FALSE
  }
  extend_grid = TRUE
  if(missing(extend) | is.null(extend)){
    extend_grid = FALSE
  } else if(length(extend) == 1) {
    extend = rep(extend,2)
  }
  s = sort(x)

  if(length(s)<3){
    s = c(s[1]-2*cutoff,s,s[length(s)]+2*cutoff)
  }
  max_l = max(s)
  min_l = min(s)

  loc_len = max_l - min_l
  if(extend_grid){
    s = c(min_l-extend[1]*loc_len,s,max_l+extend[2]*loc_len)
  }

  #merge nodes closer than cutoff
  d  = diff(s)
  n.mesh = length(s)
  while(min(d)<cutoff & n.mesh > 3){
    i = which.min(d)
    if(!extend_grid & i==1){ #make sure that mesh covers observation locations
      s = c(s[i],s[(i+1)+seq(length=n.mesh-(i+1))])
    } else if(!extend_grid & i == n.mesh-1){
      s = c(s[seq(length=(i-1))],s[i+1])
    } else {
      s = c(s[seq(length=(i-1))],mean(s[i:(i+1)]),s[(i+1)+seq(length=n.mesh-(i+1))])
    }
    d = diff(s)
    n.mesh= n.mesh-1
  }

  #add nodes to satify max.dist
  if(refine){
    while(max(d)>max.dist){
      i = which.max(d)
      si = seq(from=s[i],to=s[i+1],length.out = ceiling(d[i]/max.dist)+1)
      s = c(s[seq(length=(i-1))],si,s[(i+1)+seq(length=n.mesh-(i+1))])
      d = diff(s)
      n.mesh = length(s)
    }
  }

  if(max.dY> 0){
    nY = length(Y)
    dY= diff(Y)
    index = which( abs(dY)>max.dY)
    for(i in index){
      l1 <- loc.Y[i]
      l2 <- loc.Y[i+1]
      index2 <- l1 <= s &  s <= l2
      if(sum(index2) <= nJump){
        s <- s[!index2]
        s <- c(s, seq(l1, l2, length.out = nJump))
      }
      
      if(i +2 <= nY   ){
        l3 <- loc.Y[i+2]
        index2 <- l2 <= s &  s <= l3
        if(sum(index2) <= nJump){
          s <- s[!index2]
          s <- c(s, seq(l2, l3, length.out = nJump))
        }
        
      }
    }
    
    s <- sort(s)
    n.mesh = length(s) 
  }
  h = (c(diff(s),Inf) + c(Inf,diff(s)))/2
  h[1] <- (s[2]-s[1])/2
  h[n.mesh] <- (s[n.mesh]-s[n.mesh-1])/2

  return(list(loc = s,
              h = h,
              hs = diff(s),
              n = n.mesh))
}


#' @title Create 1d mesh from observation locations adaptively.
#'
#' @description A function to create mesh for 1d FEM discretization adaptively.
#' @param locs A list of measurement locations.
#' @param max.dist The largest distance between nodes in the mesh
#' @param cutoff Merge nodes in x that are closer than cutoff
#' @inheritParams create.meshes.1d
#' @param n.cores A numeric value for the number of cores to be used.
#'
#' @return Returns a list of output.
#'
#' @details This is a supplementary function to be used internally by other functions.
#'
#' @seealso \code{\link{create.meshes.1d}}
#'
#' @examples
#'   \dontrun{
#'   generate.adaptive.meshes.1d(...)
#'   }
#'
generate.adaptive.meshes.1d <- function(locs, 
                                        max.dist = NULL,
                                        cutoff = 1e-10,
                                        common.grid=FALSE,
                                        extend = NULL,
                                        n.cores = 1,
                                        Y = NULL,
                                        loc.Y = NULL,
                                        max.dY = -1,
                                        nJump  = 3)
{
  loc <- h <- hs <- n <- list()
  loc.length <- length(locs)
  if(common.grid || loc.length==1){
    locs <- unlist(locs)
    m <- generate.1d.mesh(x = locs,
                          max.dist = max.dist,
                          cutoff = cutoff,
                          extend = extend,
                          Y = Y[[1]],
                          loc.Y = loc.Y[[1]],
                          max.dY = max.dY,
                          nJump  = nJump)
    for(i in 1:loc.length){
      loc[[i]] <- m$loc
      h[[i]] <- m$h
      hs[[i]] <- m$hs
      n[[i]] <- m$n
    }
  } else {
    if(n.cores > 1){
      cl <- makeCluster(n.cores)
      registerDoSNOW(cl)
      clusterExport(cl, list = c('locs'),envir=environment())
      mesh.list <- foreach(i = 1:length(locs)) %dopar%
      {
        m <- generate.1d.mesh(x = locs[[i]],
                              max.dist = max.dist,
                              cutoff = cutoff,
                              extend = extend,
                              Y = Y[[i]],
                              loc.Y = loc.Y[[i]],
                              max.dY = max.dY,
                              nJump  = nJump)
        m$i = i
        return(m)
      }
      stopCluster(cl)

      h <- lapply(mesh.list,function(x) x$h)
      hs <- lapply(mesh.list,function(x) x$hs)
      n <- lapply(mesh.list,function(x) x$n)
      loc <- lapply(mesh.list,function(x) x$loc)

      ind <- unlist(lapply(mesh.list,function(x) x$i))

      h <- h[ind]
      hs <- hs[ind]
      n <- n[ind]
      loc <- loc[ind]

    } else {
      for(i in 1:length(locs)){
        m <- generate.1d.mesh(x = locs[[i]],
                              max.dist = max.dist,
                              cutoff = cutoff,
                              extend = extend,
                              Y = Y[[i]],
                              loc.Y = loc.Y[[i]],
                              max.dY = max.dY,
                              nJump  = nJump)
        loc[[i]] <- m$loc
        h[[i]] <- m$h
        hs[[i]] <- m$hs
        n[[i]] <- m$n
      }
    }
  }
  return(list(loc=loc,h=h,hs=hs,n=n))
}
