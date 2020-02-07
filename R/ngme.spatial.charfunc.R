

#' Create margninal moment and density for ngme.spatial object
#' returns density and moments of various components of the model
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
    VX1 <- integrate(f_1,0,range1)$value*(solve(B)[1,1]^2+solve(B)[1,2]^2)
    VX2 <- integrate(f_2,0,range2)$value*(solve(B)[2,2]^2+solve(B)[2,1]^2)
    if(process_list$noise=="MultiGH"){
      
      dens.p <- density_2d_nig_multivariate(c(-sdx*sqrt(VX1),sdx*sqrt(VX1)),
                                              c(-sdy*sqrt(VX2),sdy*sqrt(VX2)),
                                              n2x,
                                              n2y,
                                              c(0,process_list$mu[1], 1, exp(process_list$nu[1])),
                                              c(0,process_list$mu[2], 1, exp(process_list$nu[2])), 
                                        maternParam)
    }else{
      dens.p <- density_2d_normal_multivariate(c(-4*sqrt(VX1),4*sqrt(VX1)),
                                          c(-4*sqrt(VX2),4*sqrt(VX2)),
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