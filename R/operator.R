spde.A <- function(loc,x,right.boundary='neumann',left.boundary='neumann')
{
  if(min(loc)< min(x) || max(loc) > max(x))
    stop("locations outside support of basis")

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

spde.basis <- function(x,right.boundary = 'neumann',left.boundary = 'neumann')
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



create_operator <- function(locs,
                            n,
                            name = "matern",
                            right.boundary = 'neumann',
                            left.boundary='neumann',
                            common.grid = TRUE,
                            extend)
{
  if(tolower(name) == "matern"){
    return(create_matrices_Matern(locs,
                                  n,
                                  right.boundary = right.boundary,
                                  left.boundary = left.boundary,
                                  common.grid,
                                  extend))
  }else{
    return(create_matrices_FD2(locs,
                               n,
                               common.grid = common.grid,
                               extend))
  }

}
#' creates matrices for Matern 1D operator
#'
#' @return operator_List list to to use in simulation and estimation object
create_matrices_Matern <- function(locs,
                                   n,
                                   right.boundary = 'neumann',
                                   left.boundary='neumann',
                                   common.grid,
                                   extend)
{
  meshes <- create.meshes.1d(locs,n,common.grid,extend)
  operator_List <- list()
  if(common.grid || length(locs) == 1){

    P <- seq(min_l, max_l, length.out = n)
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
                        common.grid)
  return(operator_List)
}
#' creates matrices for Finite difference operator, one sided
#'
create_matrices_FD2 <- function(locs,
                                n,
                                right.boundary = 'neumann',
                                left.boundary='neumann',
                                common.grid,
                                extend)
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
                          common.grid)
  return(operator_List)
}

create.meshes.1d <- function(locs,n,common.grid,extend)
{
  loc <- h <- list()
  if(missing(extend)){
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
