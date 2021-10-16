

#' @title Marginal moments and density for ngme.spatial object
#' 
#' @description Returns density and moments of various components of the model
#' in the form of mean, variance, skewness and kuriotos
#' @param  obj - ngme.spatial
#' @param  n2x - number 2^n2x grid point in char func calc
#' @param  n2y - number 2^n2y grid point in char func calc (multivariate only)
#' @param  sdx - number of +- standard devations to calc the density around
#' @param  sdy - number of +- standard devations to calc the density around
#'
ngme.spatial.moment <-function(obj, n2x = 7, n2y = 7, sdx = 6, sdy= 6){
  
  stopifnot(inherits(obj, "ngme.spatial"))
  operator_list <- obj$operator_list
  process_list  <- obj$processes_list
  res <- list()
  dens <- NULL
  if(grepl('bivariate',operator_list$type)){
    
    maternParam <- list()
    maternParam[[1]] <- c(2, operator_list$tau1, operator_list$kappa1)
    maternParam[[2]] <- c(2, operator_list$tau2, operator_list$kappa2)
    maternParam[[3]] <- c(operator_list$rho,operator_list$theta)
    
    # creates function to compute range
    f_1 = function(x){maternkernelMulti(x,2,operator_list$tau1,operator_list$kappa1,2)}
    range1 = uniroot(function(x) {abs(f_1(x)/f_1(0)) - 10^(-6)}, c(0, 10^20))$root
    f_2 = function(x){maternkernelMulti(x,2,operator_list$tau2,operator_list$kappa2,2)}
    range2 = uniroot(function(x) {abs(f_2(x)/f_2(0)) - 10^(-6)}, c(0, 10^20))$root
    # radial integration 
    f_1 = function(x){2*pi*x*maternkernelMulti(x,2,operator_list$tau1,operator_list$kappa1,2)^2}
    f_2 = function(x){2*pi*x*maternkernelMulti(x,2,operator_list$tau2,operator_list$kappa2,2)^2}
    rho = maternParam[[3]][1]
    theta = maternParam[[3]][2]
    B = c(cos(theta)+rho*sin(theta) , -sin(theta)*sqrt(1+rho^2)  ,
          sin(theta)-rho*cos(theta),cos(theta)*sqrt(1+rho^2))
    B = matrix(B, nrow=2,ncol=2)
    Binv = solve(B)
    VX1 <- integrate(f_1,0,range1)$value*(solve(B)[1,1]^2+solve(B)[2,1]^2)
    VX2 <- integrate(f_2,0,range2)$value*(solve(B)[2,2]^2+solve(B)[1,2]^2)
    if(process_list$noise=="MultiGH"){
      
      dens.p <- density_2d_nig_multivariate(c(-sdx*sqrt(VX1),sdx*sqrt(VX1)),
                                            c(-sdy*sqrt(VX2),sdy*sqrt(VX2)),
                                            n2x,
                                            n2y,
                                            c(0,process_list$mu[1], 1, exp(process_list$nu[1])),
                                            c(0,process_list$mu[2], 1, exp(process_list$nu[2])), 
                                            maternParam)
    }else{
      dens.p <- density_2d_normal_multivariate(c(-sdx*sqrt(VX1),sdx*sqrt(VX1)),
                                               c(-sdy*sqrt(VX2),sdy*sqrt(VX2)),
                                               n2x,
                                               n2y,
                                               c(0, 1),
                                               c(0, 1), 
                                               maternParam)
    }
    moment.p <- calc_moment_density_2d(dens.p)
    moment_raw.p <- calc_moment_density_2d(dens.p, central=T)
    if(grepl('Normal',obj$measurementError_list$noise)){
      sigma <- exp(obj$measurementError_list$theta)
      moment.e <- list()
      moment.e[[1]] <- c(0, sigma[1]^2,0,3)
      moment.e[[2]] <- c(0, sigma[2]^2,0,3)
      moment_raw.e <- list()
      moment_raw.e[[1]] <- c(0, sigma[1]^2,0,3*sigma[1]^4)
      moment_raw.e[[2]] <- c(0, sigma[2]^2,0,3*sigma[2]^4)
    }
    moment <- list()
    moment[[1]] <- moment_transform(sum_moment(moment_raw.p[[1]],moment_raw.e[[1]]))
    moment[[2]] <- moment_transform(sum_moment(moment_raw.p[[2]],moment_raw.e[[2]]))
  }else{
    # moment and desnity for normal
    VX1 <- 0
    if(grepl('Normal',obj$measurementError_list$noise)){
      sigma <- obj$measurementError_list$sigma
      moment.e     <-  c(0, sigma^2,0,3)
      moment_raw.e <- c(0, sigma[1]^2,0,3*sigma[1]^4)
      VX1 <- sigma^2
    }
    
    maternParam <- list()
    maternParam <- c(2, operator_list$tau, operator_list$kappa)
    f_1 = function(x){maternkernel(x,2,operator_list$tau,operator_list$kappa,2)}
    range = uniroot(function(x) {abs(f_1(x)/f_1(0)) - 10^(-6)}, c(0, 10^20))$root
    f_1 = function(x){2*pi*x*maternkernel(x,2,operator_list$tau,operator_list$kappa,2)^2}
    VX1 <- VX1 + integrate(f_1,0,range)$value
    if(process_list$noise=="NIG"){
      
      dens.p <- density_2d_nig(c(-sdx*sqrt(VX1),sdx*sqrt(VX1)),
                               n2x,
                               c(0,process_list$mu, 1, process_list$nu), 
                               maternParam)
    }else{
      dens.p <- density_2d_normal(c(-sdx*sqrt(VX1),sdx*sqrt(VX1)),
                                  n2x,
                                  c(0, 1), 
                                  maternParam)
    } 
    moment.p <- calc_moment_density(dens.p)
    moment_raw.p <- calc_moment_density(dens.p, central=T)
    
    moment <- moment_transform(sum_moment(moment_raw.p,moment_raw.e))
    if(grepl('Normal',obj$measurementError_list$noise)){
      dens <- characteristic_function_to_density(
        function(t,mu=0,sigma=obj$measurementError_list$sigma) 
        {1i*t*mu - sigma^2/2*t^2 },
        n2x,
        c(-sdx*sqrt(VX1),sdx*sqrt(VX1)),
        logphi_prev = dens.p$logphi
      )
    }
    
  }
  res$dens                    <- dens
  res$moment                  <- moment
  res$measurementError$moment <- moment.e
  res$process$dens            <- dens.p
  res$process$moment          <- moment.p
  return(res)
}

#' @title Density for bivariate NIG random fields on R2
#' 
#' @param interv_x        - (2x1) [a_x,b_x]
#' @param interv_y       - (2x1) [a_y,b_y]
#' @param n2_x           - (1) power of 2 for number of point evaluations of x dir
#' @param n2_y           - (1) power of 2 for number of point evaluations of y dir
#' @param  param_x       - 
#'                      - (1) delta - location   parameter
#'                      - (2) mu    - assymteric parameter
#'                      - (3) sigma - scale      parameter
#'                      - (4) nu    - shape      parameter
#' @param  param_y       - 
#'                      - (1) delta - location   parameter
#'                      - (2) mu    - assymteric parameter
#'                      - (3) sigma - scale      parameter
#'                      - (4) nu    - shape      parameter
#' @param  maternParam - list
#'                          [[1]] (3 x 1) -   nu_x, tau_x, kappa_x
#'                          [[2]] (3 x 1) -   nu_y, tau_y, kappa_y
#'                          [[3]] (2)    -   rho, theta
##
density_2d_nig_multivariate<-function(interv_x,
                                      interv_y,
                                      n2_x,
                                      n2_y,
                                      param_x,
                                      param_y, 
                                      maternParam){

  rho = maternParam[[3]][1]
  theta = maternParam[[3]][2]
  B = c(cos(theta)+rho*sin(theta) , -sin(theta)*sqrt(1+rho^2)  ,
        sin(theta)-rho*cos(theta),cos(theta)*sqrt(1+rho^2))
  B = matrix(B, nrow=2,ncol=2)
  Binv = solve(B)
  logphi_Wx <- function(x, u, f){2*pi*x*logchar_nig_f_eval(x,u, param_x, f)}
  logphi_Wy <- function(x, u, f){2*pi*x*logchar_nig_f_eval(x,u, param_y, f)}
  
  logphi_1 <- function(t){ matern1  = maternParam[[1]];
  matern1[2] = matern1[2]/Binv[1,1];
  logchar_dens_matern(t, 
                      logphi_W  = logphi_Wx, 
                      maternParam = matern1,
                      d= 2,
                      typeMatern = 2)}
  logphi_W12 <- function(x, u_1, u_2, f_x,f_y){2*pi*x*logchar_nig_multi_f_eval(x,
                                                                                  u_1,
                                                                                  u_2,
                                                                                  param_x,
                                                                                  f_x,
                                                                                  f_y)}
  logphi_W21 <- function(x, u_1, u_2, f_x,f_y){2*pi*x*logchar_nig_multi_f_eval(x,
                                                                                  u_1,
                                                                                  u_2,
                                                                                  param_y,
                                                                                  f_x,
                                                                                  f_y)}
  logphi_12 <- function(t_x,t_y){
    matern1  = maternParam[[1]];
    matern1[2] = matern1[2]/Binv[1,1];
    matern2  = maternParam[[2]];
    matern2[2] = matern2[2]/Binv[1,2];
    logchar_dens_multi_matern(t_x,
                              t_y,
                              logphi_W12,
                              maternParam_x = matern1,
                              maternParam_y = matern2, 
                              2)}
  logphi_2 <- function(t){matern2  = maternParam[[2]];
  matern2[2] = matern2[2]/Binv[2,2];
  logchar_dens_matern(t, 
                      logphi_W  = logphi_Wy, 
                      maternParam = matern2,
                      d= 2,
                      typeMatern = 2)}
  logphi_21 <- function(t_x,t_y){
    matern1  = maternParam[[2]];
    matern1[2] = matern1[2]/Binv[2,2];
    matern2  = maternParam[[1]];
    matern2[2] = matern2[2]/Binv[2,1];
    logchar_dens_multi_matern(t_x,
                              t_y,
                              logphi_W21,
                              maternParam_x = matern1,
                              maternParam_y = matern2, 
                              2)}
  
  logphi <- function(t_x, t_y){
    n_x = length(t_x)
    n_y = length(t_y)
    if(abs(Binv[1,2]) < 10^-15){
      phi_1 <- logphi_1(t_x)
      phi_1 <- kronecker(rep(1,n_y), t(phi_1))
    }else{
      phi_1 <- t(matrix(logphi_12(t_x, t_y), nrow=n_x, ncol=n_y))
    }
    if(abs(Binv[2,1]) < 10^-15){
      phi_2 <- logphi_2(t_y)
      phi_2 <- kronecker(t(rep(1,n_x)), phi_2)
    }else{
      phi_2 <- matrix(logphi_21(t_y, t_x), nrow=n_y, ncol=n_x)
    }
    return(phi_1 + phi_2)
  }
  return(characteristic_function_to_density2d(logphi,
                                              n2_x,
                                              n2_y,
                                              interv_x,
                                              interv_y))
}


#' @title Density for bivariate Gaussian random fields on a square in R2
#' 
#' @param interv_x        - (2x1) [a_x,b_x]
#' @param interv_y       - (2x1) [a_y,b_y]
#' @param n2_x           - (1) power of 2 for number of point evaluations of x dir
#' @param n2_y           - (1) power of 2 for number of point evaluations of y dir
#' @param  param_x       - (2 x 1)  mu, sigma
#' @param  param_y       - (2 x 1)  mu, sigma
#' @param  maternParam - list
#'                          [[1]] (3 x 1) -   nu_x, tau_x, kappa_x
#'                          [[2]] (3 x 1) -   nu_y, tau_y, kappa_y
#'                          [[3]] (2)    -   rho, theta
##
density_2d_normal_multivariate<-function(interv_x,
                                         interv_y,
                                        n2_x,
                                        n2_y,
                                        param_x,
                                        param_y, 
                                        maternParam){
  #logchar_normal_f_eval
  rho = maternParam[[3]][1]
  theta = maternParam[[3]][2]
  B = c(cos(theta)+rho*sin(theta) , -sin(theta)*sqrt(1+rho^2)  ,
        sin(theta)-rho*cos(theta),cos(theta)*sqrt(1+rho^2))
  B = matrix(B, nrow=2,ncol=2)
  Binv = solve(B)
  logphi_Wx <- function(x, u, f){2*pi*x*logchar_normal_f_eval(x,u, param_x, f)}
  logphi_Wy <- function(x, u, f){2*pi*x*logchar_normal_f_eval(x,u, param_y, f)}
  
  logphi_1 <- function(t){ matern1  = maternParam[[1]];
                          matern1[2] = matern1[2]/Binv[1,1];
                          logchar_dens_matern(t, 
                                             logphi_W  = logphi_Wx, 
                                             maternParam = matern1,
                                             d= 2,
                                             typeMatern = 2)}
  logphi_W12 <- function(x, u_1, u_2, f_x,f_y){2*pi*x*logchar_normal_multi_f_eval(x,
                                                                                  u_1,
                                                                                  u_2,
                                                                                  param_x,
                                                                                  f_x,
                                                                                  f_y)}
  logphi_W21 <- function(x, u_1, u_2, f_x,f_y){2*pi*x*logchar_normal_multi_f_eval(x,
                                                                                  u_1,
                                                                                  u_2,
                                                                                  param_y,
                                                                                  f_x,
                                                                                  f_y)}
  logphi_12 <- function(t_x,t_y){
                                matern1  = maternParam[[1]];
                                matern1[2] = matern1[2]/Binv[1,1];
                                matern2  = maternParam[[2]];
                                matern2[2] = matern2[2]/Binv[1,2];
                                 logchar_dens_multi_matern(t_x,
                                                           t_y,
                                                           logphi_W12,
                                                           maternParam_x = matern1,
                                                           maternParam_y = matern2, 
                                                           2)}
  logphi_2 <- function(t){matern2  = maternParam[[2]];
                          matern2[2] = matern2[2]/Binv[2,2];
                            logchar_dens_matern(t, 
                                               logphi_W  = logphi_Wy, 
                                               maternParam = matern2,
                                               d= 2,
                                               typeMatern = 2)}
  logphi_21 <- function(t_x,t_y){
                                      matern1  = maternParam[[2]];
                                      matern1[2] = matern1[2]/Binv[2,2];
                                      matern2  = maternParam[[1]];
                                      matern2[2] = matern2[2]/Binv[2,1];
                                      logchar_dens_multi_matern(t_x,
                                                           t_y,
                                                           logphi_W21,
                                                           maternParam_x = matern1,
                                                           maternParam_y = matern2, 
                                                           2)}

  logphi <- function(t_x, t_y){
    n_x = length(t_x)
    n_y = length(t_y)
    if(abs(Binv[1,2]) < 10^-15){
      phi_1 <- logphi_1(t_x)
      phi_1 <- kronecker(rep(1,n_y), t(phi_1))
    }else{
      phi_1 <- t(matrix(logphi_12(t_x, t_y), nrow=n_x, ncol=n_y))
    }
    if(abs(Binv[2,1]) < 10^-15){
      phi_2 <- logphi_2(t_y)
      phi_2 <- kronecker(t(rep(1,n_x)), phi_2)
    }else{
      phi_2 <- matrix(logphi_21(t_y, t_x), nrow=n_y, ncol=n_x)
    }
    return(phi_1 + phi_2)
  }
  return(characteristic_function_to_density2d(logphi,
                                              n2_x,
                                              n2_y,
                                              interv_x,
                                              interv_y))
}



#' @title Density for NIG stochastic process on an interval
#' 
#' @description Computes the density for a random field on \eqn{[interv[0], interv[1]]}
#' 
#' @param interv       - start and endpoint to evalute density
#' @param n2           - number of grid points in the interval
#' @param param        - (4 x 1)
#'                      - (1) delta - location   parameter
#'                      - (2) mu    - assymteric parameter
#'                      - (3) sigma - scale      parameter
#'                      - (4) nu    - shape      parameter
#' @param  maternParam - (3 x 1)  nu, tau, kappa
density_1d_nig <-function(interv,
                          n2,
                          param, 
                          maternParam){
  
  logphi_W <- function(x, u, f){2*logchar_nig_f_eval(x,u, param, f)}
  logphi <- function(t){ logchar_dens_matern(t, 
                                             logphi_W  = logphi_W, 
                                            maternParam = MaternParameter,
                                               1)}
  return(characteristic_function_to_density( logphi, n2, interv))
}

#' @title Density for NIG random field on R2
#' 
#' @description Computes the density for a random field on
#' 
#' 
#' @param interv       - start and endpoint to evalute density
#' @param n2           - number of grid points in the interval
#' @param param        - (4 x 1)
#'                      - (1) delta - location   parameter
#'                      - (2) mu    - assymteric parameter
#'                      - (3) sigma - scale      parameter
#'                      - (4) nu    - shape      parameter
#' @param  maternParam - (3 x 1)  nu, tau, kappa
##
density_2d_nig<-function(interv,
                         n2,
                         param, 
                         maternParam){
  #logchar_normal_f_eval
  logphi_W <- function(x, u, f){2*pi*x*logchar_nig_f_eval(x,u, param, f)}
  logphi <- function(t){ logchar_dens_matern(t, 
                                             logphi_W  = logphi_W, 
                                             maternParam = maternParam,
                                             2)}
  return(characteristic_function_to_density( logphi, n2, interv))
}


#' @title Density for bivariate Gaussian stochastic process on an interval
#' computes density for normal univarate random fields in 2d
#' 
#' @param interv       - start and endpoint to evalute density
#' @param n2           - number of grid points in the interval
#' @param  param       - (2 x 1)  mu, sigma
#' @param  maternParam - (3 x 1)  nu, tau, kappa
##
density_2d_normal<-function(interv,
                            n2,
                            param, 
                            maternParam){
  #logchar_normal_f_eval
  logphi_W <- function(x, u, f){2*pi*x*logchar_normal_f_eval(x,u, param, f)}
  logphi <- function(t){ logchar_dens_matern(t, 
                                             logphi_W  = logphi_W, 
                                             maternParam = maternParam,
                                             2)}
  return(characteristic_function_to_density( logphi, n2, interv))
}
##
#' computes density for normal univarate random fields in 1d
#' between [interv[0], interv[1]]
#' @param interv       - start and endpoint to evalute density
#' @param n2           - number of grid points in the interval
#' @param  param       - (2 x 1)  mu, sigma
#' @param  maternParam - (3 x 1)  nu, tau, kappa
##
density_1d_normal<-function(interv,
                            n2,
                            param, 
                            maternParam){
  #logchar_normal_f_eval
  logphi_W <- function(x, u, f){2*logchar_normal_f_eval(x,u, param, f)}
  logphi <- function(t){ logchar_dens_matern(t, 
                                             logphi_W  = logphi_W, 
                                             maternParam = maternParam,
                                             1)}
  return(characteristic_function_to_density( logphi, n2, interv))
}

#' @title log-characterisct function of integral
#' 
#' @description log charateristic function of \eqn{\int f(x)dW_x} assuming that f is istorpic,
#' where W_x brownian sheet.
#'
#' t           - (vector) location of chara eval
#' logphi_W    - (function)   evaluates log of char function of \eqn{f(x)*W(dx)}
#' maternParam - (3 x 1)  nu, tau, kappa
#' d           - (1 x 1)  dimension
logchar_dens_matern <- function(t, 
                                logphi_W, 
                                maternParam,
                                d,
                                typeMatern = 1){
  
  if(typeMatern==1){
  f     =  function(x) {maternkernel(x,
                        maternParam[1],
                        maternParam[2],
                        maternParam[3],
                       d)}
  }else{
    f     =  function(x) {maternkernelMulti(x,
                                       maternParam[1],
                                       maternParam[2],
                                       maternParam[3],
                                       d)}
  }
  range = uniroot(function(x) {abs(f(x)/f(0)) - 10^(-6)}, c(0, 10^20))$root

  logPhi = logchar_f_W(t,
                      logphi_W = logphi_W, 
                      f= f,
                      range = c(0,range))
  return(logPhi)
}


#' @title log-characterisct function of bivariate integral
#' 
#' @description log charateristic function of \eqn{[X,Y] = \int [f_1(x), f_2(x)]dW(x)} assuming that f is istorpic,
#' where W_x brownian sheet.
#'
#' t_y           - (vector) location of chara eval
#' t_x           - (vector) location of chara eval
#' logphi_W    - (function)   evaluates log of char function of \eqn{\int [f_1(x), f_2(x)]dW(x)}
#' maternParam_x - (3 x 1)  nu, tau, kappa
#' maternParam_y - (3 x 1)  nu, tau, kappa
#' d           - (1 x 1)  dimension
logchar_dens_multi_matern <- function(t_y,
                                t_x,
                                logphi_W, 
                                maternParam_x,
                                maternParam_y,
                                d){
  
  
  f_x     =  function(x) {maternkernelMulti(x,
                                     maternParam_x[1],
                                     maternParam_x[2],
                                     maternParam_x[3],
                                     d)}
  f_y     =  function(x) {maternkernelMulti(x,
                                       maternParam_y[1],
                                       maternParam_y[2],
                                       maternParam_y[3],
                                       d)}
  range_x = uniroot(function(x) {abs(f_x(x)/f_x(0)) - 10^(-6)}, c(0, 10^15))$root
  range_y = uniroot(function(x) {abs(f_y(x)/f_y(0)) - 10^(-6)}, c(0, 10^25))$root

  logPhi = logchar_f_multi_W(t_x,
                             t_y,
                             logphi_W = logphi_W, 
                             f_x= f_x,
                             f_y= f_y,
                           range = c(0,max(range_x,range_y)))
  return(logPhi)
}

#' @title exponent in characteristic function
#' 
#' @description Calculate exponent in characteristic function of \eqn{\int_{a}^b f(x) W(dx)}.
#'
#' @param u - point where to evalute \eqn{log(\phi(x))}
#' @param logphi_W - (function) evaluates log of Char function of f(x)*W(dx)
#' @param f        - (function) the kernel
#' @param range    - (2 x 1)    a,b 
logchar_f_W <- function(u, logphi_W, f, range){
  
  
  int_fun<- function(u){
    Real <- integrate(function(x) {Re(logphi_W( x = x,
                                                u = u,
                                                f = f))}, 
                      lower = range[1], 
                      upper = range[2])
    Imag <- integrate(function(x) {Im(logphi_W( x = x,
                                                u = u,
                                                f = f))}, 
                      lower = range[1], 
                      upper = range[2])
    return(Real$value + 1i*Imag$value)
  }
  res <- apply(as.matrix(u), 1, int_fun)
  return(res)
}


#' @title exponent in bivariate characteristic function
#' 
#' calculate exponent in characteristic function of \eqn{\int_{a}^b [f_1(x), f_2(x)] W(dx)}.
#'
#' @param u - point where to evalute \eqn{log(\phi(x))}
#' @param logphi_W - (function) evaluates log of Char function of \eqn{[f_1(x), f_2(x)]*W(dx)}
#' @param f        - (function) the kernel
#' @param range    - (2 x 1)    a,b 
logchar_f_multi_W <- function(u_x, u_y, logphi_W, f_x, f_y, range){
 
  
  int_fun<- function(u){
    u_1 = u[1]
    u_2 = u[2]
    Real <- integrate(function(x) {Re(logphi_W( x = x,
                                                u_1 = u_1,
                                                u_2 = u_2,
                                                f_x = f_x,
                                                f_y = f_y))}, 
                      lower = range[1], 
                      upper = range[2]/4,
                      subdivisions = 200L,
                      rel.tol = .Machine$double.eps^0.5)
    Real2 <- integrate(function(x) {Re(logphi_W( x = x,
                                                u_1 = u_1,
                                                u_2 = u_2,
                                                f_x = f_x,
                                                f_y = f_y))}, 
                      lower = range[2]/4, 
                      upper = range[2],
                      subdivisions = 200L,
                      rel.tol = .Machine$double.eps^0.5)
    Imag <- integrate(function(x) {Im(logphi_W( x = x,
                                                u_1 = u_1,
                                                u_2 = u_2,
                                                f_x = f_x,
                                                f_y = f_y))}, 
                      lower = range[1], 
                      upper = range[2]/4,
                      subdivisions = 200L,
                      rel.tol = .Machine$double.eps^0.5)
    Imag2 <- integrate(function(x) {Im(logphi_W( x = x,
                                                u_1 = u_1,
                                                u_2 = u_2,
                                                f_x = f_x,
                                                f_y = f_y))}, 
                      lower = range[2]/4, 
                      upper = range[2],
                      subdivisions = 200L,
                      rel.tol = .Machine$double.eps^0.5)
    return((Real$value + Real2$value ) + 1i*(Imag$value + Imag2$value ))
  }
  mesh_xy <- meshgrid(u_x, u_y)
  res <- apply(cbind(c(mesh_xy$X),c(mesh_xy$Y)), 1, int_fun)
  return(res)
}



#' @title calculate exponent in characteristic function.
#'
#' @param x -
#' @param u -
#' @param param (1) delta - location   parameter
#'            - (2) mu    - assymteric parameter
#'            - (3) sigma - scale      parameter
#'            - (4) nu    - shaper     parameter
#' @param f       - function kernel
logchar_nig_f_eval<-function(x, u, param, f){
  delta = param[1]
  mu    = param[2]
  sigma = param[3]
  nu    = param[4]
  fx = f(x)
  h = 1i * (delta-mu) *  u * fx +
      sigma * sqrt(nu)* (sqrt(nu)/sigma - sqrt(nu/sigma^2 + 
                   -2 * 1i * mu/sigma^2 *u *fx  + fx^2 * u^2))
 
  return(h)
}


#' @title calculate exponent in characteristic function for normal
#' f(x) dW(x) where W is a brownian motion
#'
#' @param x  - (real) location
#' @param u  - (real) char func value
#' @param param (2 x 1) parameter values
#'                   - (1) delta - location   parameter
#'                  -  (2) sigma - scale      parameter
#' @param f       - function kernel
###
logchar_normal_f_eval <-function(x, u, param, f){
  
  delta = param[1]
  sigma = param[2]
  fx     = f(x)
  h = 1i * delta *  u * fx - 0.5 * sigma^2 * fx^2 *u^2
  
  return(h)
}


#' @title calculate exponent in characteristic function for bivariate normal
#' 
#' [f_1(x) f_2(x)] dW(x) where W is a brownian motion
#'
#' @param x  - (real) location
#' @param u_x  - (real) char func value
#' @param u_y  - (real) char func value
#' @param param (2 x 1) parameter values
#'                   - (1) delta - location   parameter
#'                  -  (2) sigma - scale      parameter
#' @param f_x       - function kernel
#' @param f_y       - function kernel
logchar_normal_multi_f_eval <-function(x, u_x, u_y, param, f_x, f_y){
  
  delta = param[1]
  sigma = param[2]
  fy     = f_x(x)
  fx     = f_y(x)
  uf = u_x * fx + u_y * fy
  h = 1i * delta *  uf- 0.5 * sigma^2 * uf^2
  
  #h = h * 2 * pi * x
  return(h)
}


#' @title calculate exponent in characteristic function.
#' 
#' @description \eqn{[f_1(x) f_2(x)] dW(x)} where W is a brownian motion
#'
#' @param x  - (real) location
#' @param u_x  - (real) char func value
#' @param u_y  - (real) char func value
#' @param param (4 x 1)
#'            - (1) delta - location   parameter
#'            - (2) mu    - assymteric parameter
#'            - (3) sigma - scale      parameter
#'            - (4) nu    - shaper     parameter
#' @param f_x       - function kernel
#' @param f_y       - function kernel
#' @name logchar_nig_multi_f_eval
logchar_nig_multi_f_eval <- function(x, u_x, u_y, param, f_x, f_y){
  delta = param[1]
  mu    = param[2]
  sigma = param[3]
  nu    = param[4]
  fy     = f_x(x)
  fx     = f_y(x)
  uf = u_x * fx + u_y * fy
  h = 1i * (delta-mu ) * uf + sigma * sqrt(nu)* (sqrt(nu)/sigma - sqrt(nu/sigma^2 + 
                                             -2 * 1i * mu/sigma^2 *uf + uf^2))
  return(h)
}



#' @title Matern kernel
#'
#' @param x     - (n x 1) distance between points
#' @param nu    - (>0)    smoothnes
#' @param tau   - (>0)    precision parameter
#' @param kappa - (>0)    range-scale parameter
#' @param d     - (int>0) dimension of x
#'
#' @return f(x) - (n x 1) kernel evaluted at x
#' @name maternkernel
maternkernel <- function(x, alpha, tau, kappa, d)
  {
  c_k =  0.5*(tau/kappa^(1.5)) #tau/kappa
  C = 1/gamma(alpha/2)
  C = C /(c_k * (4*pi)^(d/2) * kappa^(alpha - d) )
  nu = (alpha - d)/2
  if(nu ==0)
    x[abs(x)<0.001/kappa]= 0.001/kappa
  M = 2^(1-nu) * (kappa*abs(x))^(nu) * besselK(kappa * abs(x) , nu)
  if(nu>0)
    M[x==0] = 2^(1-nu) * gamma(nu) * 2^(nu-1)
  f = C * M
  return(f)
}


#' @title kernel function producing the Matern covariance
#'
#' @param x     - (n x 1) distance between points
#' @param nu    - (>0)    smoothnes
#' @param tau   - (>0)    precision parameter
#' @param kappa - (>0)    range-scale parameter
#' @param d     - (int>0) dimension of x
#' 
#' @return f(x) - (n x 1) kernel evaluted at x
#' @name maternkernelMulti
maternkernelMulti <- function(x, alpha, tau, kappa, d)
  {
  c_k =  tau/kappa
  C = 1/gamma(alpha/2)
  C = C /(c_k * (4*pi)^(d/2) * kappa^(alpha - d) )
  nu = (alpha - d)/2
  if(nu ==0)
    x[abs(x)<0.001/kappa]= 0.001/kappa
  M = 2^(1-nu) * (kappa*abs(x))^(nu) * besselK(kappa * abs(x) , nu)
  if(nu>0)
    M[x==0] = 2^(1-nu) * gamma(nu) * 2^(nu-1)
  f = C * M
  return(f)
}


#' @title Matern correlation function
#'
#' @param x     - (n x 1) distance 
#' @param alpha - (>0)    smoothness
#' @param kappa - (>0)    range-scale parameter
#' @param d     - (int>0) dimension of x
#' 
#' @return M(x) - (n x 1) correlation for Matern
#' @name materncorr
materncorr <- function(x, alpha, kappa, d)
  {
  x = abs(x)
  nu = alpha - d/2
  f  = 2^(1-nu)/gamma(nu) * (kappa * x)^nu *  besselK(kappa * x, nu)
  f[is.na(f)==T] = 1
  return(f)
}





