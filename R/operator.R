
#' @title Observation matrix computation.
#'
#' @description Computes observation matrix linking basis functions to
#'    measurement locations.
#'
#' @param operator_list A list for the operator created using \code{"create_operator"}.
#' @param locs A list of measurement locations.
#' @param i A numerical value for the index in \code{"locs"} for
#'    which the observation matrix should be computed.
#' @return An observation matrix.
#'
#' @details This function is supplementary and internally used.
#'
#' @seealso \code{\link{estimateLong}}
#'
#' @examples
#'   \dontrun{
#'   build.A.matrix(...)
#'   }
#'

build.A.matrix <- function(operator_list, locs, i)
{
  if(operator_list$manifold == "R2"){
    if(operator_list$common.grid){
      return(A = INLA::inla.spde.make.A(mesh=operator_list$mesh[[1]],loc=locs[[i]]))
    } else {
      return(A = INLA::inla.spde.make.A(mesh=operator_list$mesh[[i]],loc=locs[[i]]))
    }
  } else {
      return(spde.A(locs[[i]],
                    operator_list$loc[[i]],
                    right.boundary = operator_list$right.boundary,
                    left.boundary = operator_list$left.boundary))
  }
}

#' @title Computes observation matrix.
#'
#' @description A function to compute observation matrix for 1D problems.
#'
#' @param loc A numeric vector of measurement locations.
#' @param x A numeric vector of node locations for the basis functions.
#' @param right.boundary A character string denoting the boundary condition
#'    for the right boundary.
#' @param left.boundary A character string denoting the boundary condition
#'    for the left boundary.
#'
#' @return Returns an observation matrix.
#'
#' @details This is a supplementary function and internally used.
#'
#' @seealso \code{\link{build.A.matrix}}
#'
#' @examples
#'   \dontrun{
#'   spde.A(...)
#'   }
#'

spde.A <- function(loc, x, right.boundary = 'neumann', left.boundary = 'neumann')
{
  if(min(loc)< min(x) || max(loc) > max(x))
    stop("locations outside support of basis")

  if(is.null(right.boundary))
    right.boundary='neumann'

  if(is.null(left.boundary))
    left.boundary='neumann'

  n.x  <- length(x)
  n.loc <- length(loc)
  i <- as.vector(cBind(1:n.loc,1:n.loc))
  j <- matrix(0,n.loc,2)
  vals <- matrix(1,n.loc,2)
  for(ii in seq_len(n.loc)){
    j[ii,1] <- sum(sum((loc[ii] - x)>=0))
    vals[ii,1] <- loc[ii] - x[j[ii,1]]
    j[ii,2] <- j[ii,1] + 1
    if(j[ii,2]<=n.x){
      vals[ii,2] <- x[j[ii,2]] - loc[ii]
    } else {
      j[ii,2] = j[ii,2] -2
    }
  }
  j <- as.vector(j)
  vals <- as.vector(matrix(1-vals/rowSums(vals)))

  A <- sparseMatrix(i=i,j=j,x=vals, dims=c(n.loc,n.x))
  if(left.boundary=='dirichlet'){
    A = A[,-1,drop=FALSE]
  }
  if(right.boundary=='dirichlet'){
    A = A[,-dim(A)[2],drop=FALSE]
  }
  return(A)
}

#' @title Compute FEM matrices.
#'
#' @description A function to compute FEM matrices for 1D (longitudinal) problems.
#'
#' @param x A numeric vector containing node locations.
#' @param right.boundary A character string denoting the
#'    boundary condition for the right boundary.
#' @param left.boundary A character string denoting the
#'    boundary condition for the left boundary.
#'
#' @return Returns a list that contains mass matrix, stiffness matrix,
#'    and some other related quantities.
#'
#' @details This is a supplementary function and internally used.
#'
#' @seealso \code{\link{create_matrices_Matern}}
#'
#' @examples
#'   \dontrun{
#'   spde.basis(...)
#'   }
#'

spde.basis <- function(x, right.boundary = 'neumann', left.boundary = 'neumann',compute.Ce = FALSE)
{
  n = length(x)
  dx = diff(x)
  d <- c(Inf,dx)
  dm1 = c(d[-1],Inf)
  #d1  = c(Inf,d[-n])

  #G = bandSparse(n=n,m=n,k=c(-1,0,1),diagonals=cBind(-1/dm1, (1/dm1 + 1/d), -1/dm1))
  G = as(sparseMatrix(i = c(1:n,2:n),j=c(1:n,1:(n-1)),x=c((1/dm1 + 1/d),-1/dm1[-n]),dims=c(n,n),symmetric=TRUE),"dgCMatrix")

  if(compute.Ce){
    Ce = bandSparse(n=n,m=n,k=c(-1,0,1),diagonals=cBind(dm1/6, (dm1+d)/3, d/6))
    Ce[1,1] <- d[2]/3
    Ce[1,2] <- d[2]/6
    Ce[n,n] <- d[n]/3
    Ce[n,n-1] <- d[n]/6
  }


  h = (c(dx,Inf)+d)/2
  h[1] <- (x[2]-x[1])/2
  h[n] <- (x[n]-x[n-1])/2
  C = sparseMatrix(i=1:n,j=1:n,x=h,dims = c(n,n))
  Ci = sparseMatrix(i=1:n,j=1:n,x=1/h,dims = c(n,n))

  if(left.boundary=='dirichlet'){
    C = C[-1,-1]
    Ci = Ci[-1,-1]
    G = G[-1,-1]
    Ce = Ce[-1,-1]
    n = n-1
    h = h[-1]
  }
  if(right.boundary=='dirichlet'){
    C = C[-n,-n]
    Ci = Ci[-n,-n]
    G = G[-n,-n]
    Ce = Ce[-n,-n]
    n = n-1
    h = h[-n]
  }
out.list <- list(G=G,
                 C=C,
                 Ci=Ci,
                 loc=x,
                 h = h)
  if(compute.Ce)
    out.list$Ce = Ce
  return(out.list)
}

#' @title Compute FEM matrices - 2D.
#'
#' @description A function to compute FEM matrices for 2D (spatial) problems.
#'
#' @param mesh An \code{inla.mesh} object.
#'
#' @details This is a supplementary function to be used internally by other functions.
#'
#' @examples
#'   \dontrun{
#'   create_operator_matern2D(...)
#'   }
#'

create_operator_matern2D <- function(mesh)
{
  INLA:::inla.require.inherits(mesh, c("inla.mesh", "inla.mesh.1d"), "'mesh'")
  fem = INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2,
                           output = list("c0", "c1", "g1", "g2","dx","dy","dz"),
                           gradients=TRUE)
  n = mesh$n
  h = (fem$c1%*%matrix(rep(1,n)))@x

  out <- list(type = "matern",
              mesh = list(mesh),
              C = list(as(fem$c1,"CsparseMatrix")),
              G = list(as(fem$g1,"CsparseMatrix")),
              Ci = list(Matrix::Diagonal(n,1/h)),
              h = list(h),
              loc   = list(mesh$loc),
              common.grid = TRUE,
              manifold ="R2")
}

#' @title Create operator components.
#'
#' @description A function to compute a list of objects for the operator.
#'
#' @param locs A numeric list of measurement locations.
#' @param n A numeric value for the number of FEM basis functions that should be used.
#' @param name A character string for the operator type,
#'    possible options are \code{"matern"} and \code{"fd2"}.
#' @param right.boundary A character string denoting the boundary condition
#'    for the right boundary.
#' @param left.boundary A character string denoting the boundary condition
#'    for the left boundary.
#' @param common.grid A logical variable for using a common grid for all subjects,
#'    \code{"TRUE"} indicates using a common grid,
#'    \code{"FALSE"} uncommon grids.
#' @param extend A numeric vector with two elements specifying the amount of extension
#'    of the grid to the left and right beyondthe measurement locations.
#'
#' @details This is a supplementary function to be used internally by other functions.
#'
#' @seealso \code{\link{estimateLong}}
#'
#' @examples
#'   \dontrun{
#'   create_operator(...)
#'   }
#'

create_operator <- function(locs,
                            n,
                            name = "matern",
                            right.boundary = 'neumann',
                            left.boundary='neumann',
                            common.grid = FALSE,
                            extend  = NULL,
                            max.dist,
                            cutoff = 1e-10,
                            n.cores = 1)
{
  if(!missing(n) & !missing(max.dist)){
    stop("Supply either n or max.dist")
  }
  if(tolower(name) == "matern"){
    return(create_matrices_Matern(locs = locs,
                                  right.boundary = right.boundary,
                                  left.boundary = left.boundary,
                                  common.grid = common.grid,
                                  extend = extend,
                                  max.dist = max.dist,
                                  cutoff = cutoff,
                                  n.cores = n.cores))
  }else{
    return(create_matrices_FD2(locs = locs,
                               common.grid = common.grid,
                               extend = extend,
                               max.dist = max.dist,
                               cutoff = cutoff))
  }

}

#' @title Create matrices for Matern.
#'
#' @description A function to create matrices for Matern 1D (longitudinal) operator.
#'
#' @param locs A numeric list for the meansurement locations.
#' @param n A numeric value for the number of FEM basis functions that should be used.
#' @param right.boundary A character string denoting the boundary condition
#'    for the right boundary.
#' @param left.boundary A character string denoting the boundary condition
#'    for the left boundary.
#' @inheritParams create_operator
#'
#' @return Returns matrices.
#'
#' @details This is a supplementary function to be used internally by other functions.
#'
#' @seealso \code{\link{create_operator}}
#'
#' @examples
#'   \dontrun{
#'   create_matrices_Matern(...)
#'   }
#'

create_matrices_Matern <- function(locs,
                                   n,
                                   right.boundary = 'neumann',
                                   left.boundary ='neumann',
                                   common.grid = FALSE,
                                   extend = NULL,
                                   max.dist,
                                   cutoff = 1e-10,
                                   n.cores = 1)
{

  meshes <- generate.adaptive.meshes.1d(locs,
                                        max.dist = max.dist,
                                        cutoff = cutoff,
                                        common.grid=common.grid,
                                        extend = extend,
                                        n.cores = n.cores)

  operator_List <- list()
  C <- Ci <- G <- Ce <- h <- list()
  if(common.grid || length(locs) == 1){
    MatrixBlock <- spde.basis(meshes$loc[[1]],right.boundary=right.boundary,left.boundary=left.boundary)
    for(i in 1:length(locs))
    {
      C[[i]] = MatrixBlock$C
      Ci[[i]] = MatrixBlock$Ci
      G[[i]] = MatrixBlock$G
      Ce[[i]] = MatrixBlock$Ce
      h[[i]] = MatrixBlock$h
    }
  } else {
    if(n.cores>1){
      cl <- makeCluster(n.cores)
      registerDoSNOW(cl)
      clusterExport(cl, list = c('locs'),envir=environment())
      b.list <- foreach(i = 1:length(locs)) %dopar%
        {
          m <- spde.basis(meshes$loc[[i]],right.boundary=right.boundary,left.boundary=left.boundary)
          m$i = i
          return(m)
        }
      stopCluster(cl)

      C <- lapply(b.list,function(x) x$C)
      Ci <- lapply(b.list,function(x) x$Ci)
      G <- lapply(b.list,function(x) x$G)
      h <- lapply(b.list,function(x) x$h)

      ind <- unlist(lapply(b.list,function(x) x$i))

      C <- C[ind]
      Ci <- Ci[ind]
      G <- G[ind]
      h <- h[ind]

    } else {
      for(i in 1:length(locs))
      {
        MatrixBlock <- spde.basis(meshes$loc[[i]],right.boundary=right.boundary,left.boundary=left.boundary)
        C[[i]] = MatrixBlock$C
        Ci[[i]] = MatrixBlock$Ci
        G[[i]] = MatrixBlock$G
        h[[i]] = MatrixBlock$h
      }
    }
  }

  operator_List <- list(type = 'Matern',
                        C = C,
                        Ci = Ci,
                        G = G,
                        h = h,
                        kappa = 0,
                        loc   = meshes$loc,
                        right.boundary=right.boundary,
                        left.boundary=left.boundary,
                        manifold = "R",
                        common.grid = common.grid)
  return(operator_List)
}

#' @title Create matrices for Finite difference.
#' @description A function to create matrices for Finite difference operator, one sided.
#' @inheritParams create_matrices_Matern
#' @return Returns matrices.
#'
#' @details This is a supplementary function to be used internally by other functions.
#'
#' @seealso \code{\link{create_operator}}
#'
#' @examples
#'   \dontrun{
#'   create_matrices_FD2(...)
#'   }
#'

create_matrices_FD2 <- function(locs,
                                common.grid = FALSE,
                                extend = NULL,
                                max.dist,
                                cutoff = 1e-10)
{

  meshes <- generate.adaptive.meshes.1d(locs,
                                        max.dist = max.dist,
                                        cutoff = cutoff,
                                        common.grid=common.grid,
                                        extend = extend)

  Q <- list()
  if(common.grid || length(locs) == 1) {
    h <- meshes$hs[[1]]
    Operator_1D <- bandSparse(n=meshes$n[[1]]-1,m=meshes$n[[1]]-1,k=c(-1,0),diagonals=cbind(-1/h,1/h))
    hOperator_1D <- bandSparse(n=meshes$n[[1]]-1,m=meshes$n[[1]]-1,k=c(-1,0),
                               diagonals=cbind(-rep(1,length(h)),rep(1,length(h))))
    Operator_2D <- (hOperator_1D) %*% Operator_1D #mult with h for the first operator is scaled
    # due to W(t+h) - W(t) = N(0,h), not (W(t+h) - W(t))/h = N(0,h)
    for(i in 1:length(locs)){
      Q[[i]] <-as(Operator_2D, "dgCMatrix")
    }
  } else {
    for(i in 1:length(locs)){
      h <- meshes$hs[[i]]
      Operator_1D <- bandSparse(n=meshes$n[[i]]-1,m=meshes$n[[i]]-1,k=c(-1,0),diagonals=cbind(-1/h,1/h))
      hOperator_1D <- bandSparse(n=meshes$n[[i]]-1,m=meshes$n[[i]]-1,k=c(-1,0),
                                 diagonals=cbind(-rep(1,length(h)),rep(1,length(h))))
      Operator_2D <- (hOperator_1D) %*% Operator_1D

      Q[[i]] <- as(Operator_2D, "dgCMatrix")
    }
  }
    operator_List <- list(type   = 'fd2',
                          Q      = Q,
                          h      = meshes$hs,
                          loc   = meshes$loc,
                          right.boundary = 'neumann',
                          left.boundary='dirichlet',
                          manifold = "R",
                          common.grid = common.grid)
  return(operator_List)
}


create_matrices_FD2_fem <- function(locs,
                                    common.grid = TRUE,
                                    extend = NULL,
                                    max.dist,
                                    cutoff = 1e-10)
{

  meshes <- generate.adaptive.meshes.1d(locs=locs,
                                        max.dist = max.dist,
                                        cutoff = cutoff,
                                        common.grid=common.grid,
                                        extend = extend)

  Q <- h <- list()
  if(common.grid || length(locs) == 1) {
    MatrixBlock <- spde.basis(meshes$loc[[1]],left.boundary = "dirichlet")
    for(i in 1:length(locs)){
      #Q[[1]] = as(MatrixBlock$B - MatrixBlock$G,"dgCMatrix")
      Q[[i]] = as(-MatrixBlock$G,"dgCMatrix")
      h[[i]] <- MatrixBlock$h
      #Q[[1]] = Matrix::t(L)%*%MatrixBlock$Ci%*%L,"dgCMatrix")
    }
  } else {
    n = unlist(meshes$n)
    if(length(n) == 1){
      n = rep(n,length(locs))
    }
    for(i in 1:length(locs)){
      MatrixBlock <- spde.basis(meshes$loc[[i]],left.boundary = "dirichlet")
      #L = MatrixBlock$B - MatrixBlock$G
      #Q[[i]] = as(MatrixBlock$B - MatrixBlock$G,"dgCMatrix")
      Q[[i]] = as(-MatrixBlock$G,"dgCMatrix")
      h[[i]] <- MatrixBlock$h
      #Q[[i]] = as(t(L)%*%MatrixBlock$Ci%*%L,"dgCMatrix")
    }
  }
  operator_List <- list(type   = 'fd2',
                        Q      = Q,
                        h      = h,
                        loc   = meshes$loc,
                        right.boundary = 'dirichlet',
                        left.boundary='neumann',
                        manifold = "R",
                        common.grid = common.grid)
  return(operator_List)
}
