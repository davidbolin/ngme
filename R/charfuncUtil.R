#'
#'  computes the centralized moment of X+Y
#'  @param moment1 - (4 x 1) four centralized moments of X
#'  @param moment2 - (4 x 1) four centralized moments of Y
#'
#'
#'
sum_moment <- function(moment1, moment2){
  M1 <- moment1[1] + moment2[1]
  M2 <- moment1[2] + moment2[2]
  M3 <- moment1[3] + moment2[3]
  M4 <- moment1[4] + moment2[4] + 6 * moment1[2]*moment2[2]
  return(c(M1,M2,M3,M4))
}
#'
#' computes mean, variance, skewness, kurtosis from the four moments
#'
#'
moment_transform <- function(moment){
  EX <- moment[1]
  VX <- moment[2] 
  SX <- moment[3]/moment[2]^(3/2)
  KX <- moment[4]/moment[2]^2
  return(c(EX,VX,SX,KX))
}


##
#' takes output from characteristic_function_to_density and computes
#' mean, variance, skewness, kurtosis
#' 
#' @param dens - (list) output from characteristic_function_to_density
#' 
#' @return  vec - (4x1) expectiation, variance, skewness, kurtosis
##
calc_moment_density <- function(dens, central=F){
  
  h <- dens$x[2]-dens$x[1]
  M1 <- sum(dens$x*dens$density)*h
  M2 <- sum((dens$x-M1)^2*dens$density)*h 
  M3 <- sum((dens$x-M1)^3*dens$density)*h 
  M4 <- sum((dens$x-M1)^4*dens$density)*h 
  if(central)
    return(c(M1,M2,M3,M4))
  return(moment_transform(c(M1,M2,M3,M4)))
}

##
#' takes output from characteristic_function_to_density and computes
#' mean, variance, skewness, kurtosis
#' 
#' @param dens - (list) output from characteristic_function_to_density
#' 
#' @return  vec - list [[1]] dim one expectiation, variance, skewness, kurtosis
#'                     [[2]] dim two expectiation, variance, skewness, kurtosis
#'                     [[3]] covariance
##
calc_moment_density_2d <- function(dens, central=F){
  res <- list()
  dens1 <- list(x = dens$x, density = colSums(dens$density)*(dens$y[2]-dens$y[1]))
  dens2 <- list(x = dens$y, density = rowSums(dens$density)*(dens$x[2]-dens$x[1]))
  res[[1]] <- calc_moment_density(dens1,central)
  res[[2]] <- calc_moment_density(dens2,central)
  
  mesh_xy <- meshgrid(dens$x - res[[1]][1],
                      dens$y - res[[2]][1])
  h1 <- dens$x[2]-dens$x[1]
  h2 <- dens$y[2]-dens$y[1]
  res[[3]] <- sum(mesh_xy$X*mesh_xy$Y*dens$density) * h1 * h2
  return(res)
}

##
#' takes sample vector and computes
#' mean, variance, skewness, kurtosis
#' 
#' @param  X  - (n x 1) output from iid sample of densit
#' 
#' @return  vec - (4x1) emperical expectiation, variance, skewness, kurtosis
##
calc_moment_emperical <- function(X, central=F){
  
  M1 <- mean(X)
  M2 <- mean((X-M1)^2)
  M3 <- mean((X-M1)^3)
  M4 <- mean((X-M1)^4)
  if(central)
    return(c(M1,M2,M3,M4))
  return(moment_transform(c(M1,M2,M3,M4)))
}

##
#' log density from char func using fft
#'
#' @param logphi   - (function) should return expontial term of the char func
#' @param n2       - (1x1) power of 2 for number of point evaluations
#' @param interv   - (2x1) evalution interval 
###
characteristic_function_to_density <- function( logphi, n2, interv){
  n = 2**n2
  a = interv[1]
  b = interv[2]
  i <- 0:(n-1)            # Indices
  dx <- (b-a)/n           # Step size, for the density
  x <- a + i * dx         # Grid, for the density
  dt <- 2*pi / ( n * dx ) # Step size, frequency space
  c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
  d <-  n/2 * dt          # (center the interval on zero)
  t <- c + i * dt         # Grid, frequency space
  
  #no idea what I am doing below but it works!
  X <- exp( -(0+1i) * i * dt * a + logphi(t)) 
  Y <- fft(X)
  density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
  return(data.frame(x = x, density = Re(density)))
}

meshgrid <- function(x, y = x) {
  if (!is.numeric(x) || !is.numeric(y))
    stop("Arguments 'x' and 'y' must be numeric vectors.")
  
  x <- c(x); y <- c(y)
  n <- length(x)
  m <- length(y)
  
  X <- matrix(rep(x, each = m),  nrow = m, ncol = n)
  Y <- matrix(rep(y, times = n), nrow = m, ncol = n)
  
  return(list(X = X, Y = Y))
}

##
#' log density from char func using fft \phi(x,y)
#' on the grid [a_x,b_x] x [a_y, b_y]
#'
#' @param logphi   - (function) logphi(x,y) should return log( \phi(x,y)) 
#'                              the grid of kron(x,1) and kron(y,1)
#' @param n2_x     - (1) power of 2 for number of point evaluations of x dir
#' @param n2_y     - (1) power of 2 for number of point evaluations of y dir
#' @param interv_x - (2x1) [a_x,b_x]
#' @param interv_y - (2x1) [a_y,b_y]
###
characteristic_function_to_density2d <- function( logphi,
                                                  n2_x,
                                                  n2_y,
                                                  interv_x,
                                                  interv_y){
  n_x = 2**n2_x
  n_y = 2**n2_y
  a_x = interv_x[1]
  b_x = interv_x[2]
  a_y = interv_y[1]
  b_y = interv_y[2]
  # x - grid
  i_x <- 0:(n_x-1)            # Indices
  dx <- (b_x-a_x)/n_x         # Step size, for the density
  x <- a_x + i_x * dx           # Grid, for the density
  dt_x <- 2*pi / ( n_x * dx ) # Step size, frequency space
  c_x <- -n_x/2 * dt_x          # Evaluate the characteristic function on [c,d]
  d_x <-  n_x/2 * dt_x          # (center the interval on zero)
  t_x <- c_x + i_x * dt_x     # Grid, frequency space
  # y grid
  i_y   <- 0:(n_y-1)            # Indices
  dy    <- (b_y-a_y)/n_y         # Step size, for the density
  y     <- a_y + i_y * dy           # Grid, for the density
  dt_y  <- 2*pi / ( n_y * dy ) # Step size, frequency space
  c_y   <- -n_y/2 * dt_y          # Evaluate the characteristic function on [c,d]
  d_y   <-  n_y/2 * dt_y          # (center the interval on zero)
  t_y   <- c_y + i_y * dt_y     # Grid, frequency space
  
  mesh_i <- meshgrid(i_x, i_y)
  #no idea what I am doing below but it works!
  X <- exp(-(0+1i) * mesh_i$Y * dt_y * a_y
           -(0+1i) * mesh_i$X * dt_x * a_x 
           + logphi(t_x, t_y))
  mesh_xy <- meshgrid(x,y)
  density <- fft(X)
  density <- dt_y / (2*pi) * exp( - (0+1i) * c_x * mesh_xy$X ) * density
  density <- dt_x / (2*pi) * exp( - (0+1i) * c_y * mesh_xy$Y ) * density
  return(list(y = y,x = x, density = Re(density)))
}

