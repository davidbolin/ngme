
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
    if(operator_list$common.grid){
      return(spde.A(locs[[i]],
                    operator_list$loc[[1]],
                    right.boundary = operator_list$right.boundary,
                    left.boundary = operator_list$left.boundary))
    }  else {
      return(spde.A(locs[[i]],
                    operator_list$loc[[i]],
                    right.boundary = operator_list$right.boundary,
                    left.boundary = operator_list$left.boundary))
    }
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
    A = A[,-1]
  }
  if(right.boundary=='dirichlet'){
    A = A[,-dim(A)[2]]
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

spde.basis <- function(x, right.boundary = 'neumann', left.boundary = 'neumann')
{
  n = length(x)
  d <- c(Inf,diff(x))
  dm1 = c(d[2:n],Inf)
  d1  = c(Inf,d[1:(n-1)])

  D = cBind(1/dm1, -(1/dm1 + 1/d), 1/dm1)
  G = -bandSparse(n=n,m=n,k=c(-1,0,1),diagonals=D)

  D = cBind(dm1/6, (dm1+d)/3, d/6)
  Ce = bandSparse(n=n,m=n,k=c(-1,0,1),diagonals=D)
  Ce[1,1] <- d[2]/3
  Ce[1,2] <- d[2]/6
  Ce[n,n] <- d[n]/3
  Ce[n,n-1] <- d[n]/6

  h = (c(diff(x),Inf)+d)/2
  h[1] <- (x[2]-x[1])/2
  h[n] <- (x[n]-x[n-1])/2
  C = Diagonal(n,h)
  Ci = Diagonal(n,1/h)

  if(left.boundary=='dirichlet'){
    C = C[-1,-1]
    Ci = Ci[-1,-1]
    G = G[-1,-1]
    Ce = Ce[-1,-1]
    n = n-1
  }
  if(right.boundary=='dirichlet'){
    C = C[-n,-n]
    Ci = Ci[-n,-n]
    G = G[-n,-n]
    Ce = Ce[-n,-n]
    n = n-1
  }

  return(list(G=G,
              C=C,
              Ce=Ce,
              Ci=Ci,
              loc=x,
              h = h))
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
#' @seealso \link{\code{estimateLong}}
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
                            common.grid = TRUE,
                            extend  = NULL)
{
  if(tolower(name) == "matern"){
    return(create_matrices_Matern(locs = locs,
                                  n = n,
                                  right.boundary = right.boundary,
                                  left.boundary = left.boundary,
                                  common.grid = common.grid,
                                  extend = extend))
  }else{
    return(create_matrices_FD2(locs = locs,
                               n = n,
                               common.grid = common.grid,
                               extend = extend))
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
#' @param common.grid should a common grid be used for all subjects?
#' @param extend should the grid be extended beyond the measurement locations? 
#' vector with two values specifying the amount to extend in each direction.
#' 
#' @return operator_List
#' #' @details This is a supplementary function to be used internally by other functions.
#' 
#' @seealso \link{\code{create_operator}}
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
                                   common.grid = TRUE,
                                   extend = NULL)
{
  meshes <- create.meshes.1d(locs,n,common.grid,extend)
  operator_List <- list()
  if(common.grid || length(locs) == 1){
    MatrixBlock <- spde.basis(meshes$loc[[1]],right.boundary=right.boundary,left.boundary=left.boundary)
    C = list(as(as(MatrixBlock$C,"CsparseMatrix"), "dgCMatrix"))
    Ci = list(as(as(MatrixBlock$Ci,"CsparseMatrix"), "dgCMatrix"))
    G = list(MatrixBlock$G)
    Ce = list(MatrixBlock$Ce)
  } else {
    C <- Ci <- G <- Ce <- list()
    if(length(n) == 1){
      n <- rep(n,length(locs))
    }
    for(i in 1:length(locs))
    {
      MatrixBlock <- spde.basis(meshes$loc[[i]],right.boundary=right.boundary,left.boundary=left.boundary)
      C[[i]] = as(as(MatrixBlock$C,"CsparseMatrix"), "dgCMatrix")
      Ci[[i]] = as(as(MatrixBlock$Ci,"CsparseMatrix"), "dgCMatrix")
      G[[i]] = MatrixBlock$G
      Ce[[i]] = MatrixBlock$Ce
    }
  }

  operator_List <- list(type = 'Matern',
                        C = C,
                        Ci = Ci,
                        G = G,
                        Ce = Ce,
                        h = meshes$h,
                        kappa = 0,
                        loc   = meshes$loc,
                        right.boundary=right.boundary,
                        left.boundary=left.boundary,
                        manifold = "R",
                        common.grid = common.grid)
  return(operator_List)
}

#' creates matrices for Finite difference operator, one sided
#'@param locs meansurement locations
#'@param n number of FEM basis functions that should be used
#' @param right.boundary boundary condition at right boundary
#' @param left.boundary boundary condition at left boundary
#' @param common.grid should a common grid be used for all subjects?
#' @param extend should the grid be extended beyond the measurement locations? vector with two values specifying the amount to extend in each direction.
#' @return operator_List

create_matrices_FD2 <- function(locs,
                                n,
                                right.boundary = 'neumann',
                                left.boundary='neumann',
                                common.grid = TRUE,
                                extend = NULL)
{
  meshes <- create.meshes.1d(locs,n,common.grid,extend)
  Q <- list()
  if(common.grid || length(locs) == 1) {
    vec_toeplitz <- rep(0, length=n)
    h <- meshes$h[[1]][1]
    vec_toeplitz[1] <- -1 / h # -1/h
    vec_toeplitz[2] <- 1  / h  # 1/h
    Operator_1D <- Matrix(toeplitz(vec_toeplitz), sparse=T)
    Operator_1D[upper.tri(Operator_1D)] = 0
    Operator_2D <- (h*Operator_1D) %*% Operator_1D #mult with h for the first operator is scaled
    # due to W(t+h) - W(t) = N(0,h), not (W(t+h) - W(t))/h = N(0,h)
    Q[[1]] <-as(Operator_2D, "dgCMatrix")
  } else {
    if(length(n) == 1){
      n = rep(n,length(locs))
    }
    for(i in 1:length(locs)){
      vec_toeplitz <- rep(0, length=n[i])
      h <- meshes$h[[i]][1]
      vec_toeplitz[1] <- -1 / h
      vec_toeplitz[2] <- 1  / h
      Operator_1D <- Matrix(toeplitz(vec_toeplitz), sparse=T)
      Operator_1D[upper.tri(Operator_1D)] = 0
      Operator_2D <- (h*Operator_1D) %*% Operator_1D
      Q[[i]] <- as(Operator_2D, "dgCMatrix")
    }
  }
    operator_List <- list(type   = 'fd2',
                          Q      = Q,
                          h      = meshes$h,
                          loc   = meshes$loc,
                          right.boundary = 'neumann',
                          left.boundary='neumann',
                          manifold = "R",
                          common.grid = common.grid)
  return(operator_List)
}

#' creates mesh for FEM discretization
#'@param locs meansurement locations
#'@param n number of FEM basis functions that should be used
#' @param common.grid should a common grid be used for all subjects?
#' @param extend should the grid be extended beyond the measurement locations? vector with two values specifying the amount to extend in each direction.
#' @return list of meshes

create.meshes.1d <- function(locs,n,common.grid,extend = NULL)
{
  loc <- h <- list()
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
    h[[1]] <- rep(loc[[1]][2] - loc[[1]][1],n)
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
      h[[i]] <- rep(loc[[i]][2] - loc[[i]][1],n[i])
    }
  }
  return(list(loc = loc,h = h))
}
