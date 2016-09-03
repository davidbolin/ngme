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
                            common.grid = TRUE)
{
  if(tolower(name) == "matern"){
    return(create_matrices_Matern(locs, n, right.boundary = right.boundary,left.boundary = left.boundary,common.grid))
  }else{
    return(create_matrices_FD2(locs, n, common.grid = common.grid))
  }

}
#' creates matrices for Matern 1D operator
#'
#' @return operator_List list to to use in simulation and estimation object
create_matrices_Matern <- function(locs, n, right.boundary = 'neumann',left.boundary='neumann',common.grid)
{
  operator_List <- list()
  if(common.grid || length(locs) == 1){
    min_l <- min(locs[[1]])
    max_l <- max(locs[[1]])
    if(length(locs) > 1){
      for(i in 2:length(locs))
      {
        min_l <- min(min_l, min(locs[[i]]))
        max_l <- max(max_l, max(locs[[i]]))
      }
    }
    P <- seq(min_l, max_l, length.out = n)
    MatrixBlock <- spde.basis(P,right.boundary=right.boundary,left.boundary=left.boundary)
    C = list(as(as(MatrixBlock$C,"CsparseMatrix"), "dgCMatrix"))
    Ci = list(as(as(MatrixBlock$Ci,"CsparseMatrix"), "dgCMatrix"))
    G = list(MatrixBlock$G)
    Ce = list(MatrixBlock$Ce)
    h = list(MatrixBlock$h)
    loc = list(P)
  } else {
    C <- Ci <- G <- Ce <- h <- list()
    if(length(n) == 1){
      n <- rep(n,length(locs))
    }
    for(i in 1:length(locs))
    {
      P <- seq(min(locs[[i]],max(locs[[i]]),length.out = n[i]))
      MatrixBlock <- spde.basis(P,right.boundary=right.boundary,left.boundary=left.boundary)
      C[[i]] = as(as(MatrixBlock$C,"CsparseMatrix"), "dgCMatrix")
      Ci[[i]] = as(as(MatrixBlock$Ci,"CsparseMatrix"), "dgCMatrix")
      G[[i]] = MatrixBlock$G
      Ce[[i]] = MatrixBlock$Ce
      h[[i]] = MatrixBlock$h
    }
  }

  operator_List <- list(type = 'Matern',
                        C = C,
                        Ci = Ci,
                        G = G,
                        Ce = Ce,
                        h = h,
                        kappa = 0,
                        loc   = loc,
                        right.boundary=right.boundary,
                        left.boundary=left.boundary,
                        common.grid)
  return(operator_List)
}
#' creates matrices for Finite difference operator, one sided
#'
create_matrices_FD2 <- function(locs, n,right.boundary = 'neumann',left.boundary='neumann',common.grid)
{
  if(common.grid || length(locs) == 1) {
    min_l <- min(locs[[1]])
    max_l <- max(locs[[1]])
    if(length(locs) > 1){
      for(i in 2:length(locs))
      {
        min_l <- min(min_l, min(locs[[i]]))
        max_l <- max(max_l, max(locs[[i]]))
      }
    }
    P <- seq(min_l, max_l, length.out = n)
    vec_toeplitz <- rep(0, length=n)
    h <- (P[2] - P[1])
    vec_toeplitz[1] <- -1 / h # -1/h
    vec_toeplitz[2] <- 1  / h  # 1/h
    Operator_1D <- Matrix(toeplitz(vec_toeplitz), sparse=T)
    Operator_1D[upper.tri(Operator_1D)] = 0
    Operator_2D <- (h*Operator_1D) %*% Operator_1D #mult with h for the first operator is scaled
    # due to W(t+h) - W(t) = N(0,h), not (W(t+h) - W(t))/h = N(0,h)
    Q <- list(as(Operator_2D, "dgCMatrix"))
    h <- list(rep(h,n))
    loc <- list(P)
  } else {
    Q <- h <- list()
    for(i in 1:length(locs)){
      vec_toeplitz <- rep(0, length=n)
      hi <- (locs[[i]][2] - locs[[i]][1])
      vec_toeplitz[1] <- -1 / hi # -1/h
      vec_toeplitz[2] <- 1  / hi  # 1/h
      Operator_1D <- Matrix(toeplitz(vec_toeplitz), sparse=T)
      Operator_1D[upper.tri(Operator_1D)] = 0
      Operator_2D <- (hi*Operator_1D) %*% Operator_1D #mult with h for the first operator is scaled
      # due to W(t+h) - W(t) = N(0,h), not (W(t+h) - W(t))/h = N(0,h)
      Q[[i]] <- as(Operator_2D, "dgCMatrix")
      h[[i]] <- rep(hi,n)
      loc <- locs
    }
  }
    operator_List <- list(type   = 'fd2',
                          Q      = Q,
                          h      = h,
                          loc   = loc,
                          right.boundary = 'neumann',
                          left.boundary='neumann',
                          common.grid)
  return(operator_List)
}
