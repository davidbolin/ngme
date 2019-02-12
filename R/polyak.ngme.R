#' @title Post-processing of ngme parameter estimates
#'
#' @description Updates the results using Polyak averaging and estimates the MC variance of the parameter estimates.
#'
#' @param object A fitted object returned by the \code{"ngme"} function.
#' @param param A character string for which parameter set the averaging should be performed;
#'   \code{"fixed"} for fixed effects,
#'   \code{"random"} for time-invariant random effects,
#'   \code{"process"} for stochastic process,
#'   \code{"error"} for error term.
#' @param plot Plot the results or not?
#' @param polyak.rate - ([0,1]) how much much of moving average should be used
#'                      (using weighted average where weight is polyak.rate)       
#' @param burnin.end  - (int) where from start the variance calc and polyak averaging
#' @seealso \code{\link{ngme}}
#'
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   fitted(fit)
#'   }
#'

polyak.ngme <- function(object,
                        param = "fixed",
                        polyak.rate = 1,
                        plot = TRUE,
                        burnin.end = NULL){
  
  if(param == "fixed" & !is.null(object$polyaked_fixed)){
    stop("Fixed effects have already been Polyak averaged")
  }
  if(param == "random" & !is.null(object$polyaked_random)){
    stop("Random effects have already been Polyak averaged")
  }
  if(param == "process" & !is.null(object$polyaked_process)){
    stop("Process parameters have already been Polyak averaged")
  }
  if(param == "error" & !is.null(object$polyaked_error)){
    stop("Error parameters have already been Polyak averaged")
  }
  
  # if process was not specified in the model
  if(param == "process" && object$use_process == FALSE){
    stop("No process was included in the model fit")
  }
  
  if(param == "fixed"){
    
    index_fixed  <- object$index_fixed
    index_random <- object$index_random
    
    #n.random <- dim(object$mixedEffect_list$beta_random)[2]
    #n.fixed <- dim(object$mixedEffect_list$beta_fixed)[2]
    #pars <- object$fixed_est_vec[,(n.random+1):(n.random+n.fixed)]
    pars <- object$fixed_est_vec
    pars.smooth <- smooth.trajectory(pars,
                                     polyak.rate = polyak.rate, 
                                     burnin.end = burnin.end)
    colnames(pars.smooth$x) <- colnames(pars.smooth$xe) <- colnames(object$fixed_est_vec)
    
    object$mixedEffect_list$betaf_vec  <- pars.smooth$x[, index_fixed]
    object$mixedEffect_list$beta_fixed <- pars.smooth$xe[, index_fixed]
    
    object$mixedEffect_list$betar_vec   <- pars.smooth$x[, index_random]
    object$mixedEffect_list$beta_random <- pars.smooth$xe[, index_random]
    
    #update summaries (which also contain the random effects)
    #object$fixed_est_vec[,(n.random+1):(n.random+n.fixed)] <- pars.smooth$x
    #object$fixed_est[(n.random+1):(n.random+n.fixed)] <- pars.smooth$xe
    
    object$fixed_est_vec <- pars.smooth$x
    object$fixed_est     <- as.numeric(pars.smooth$xe)
    names(object$fixed_est) <- colnames(pars.smooth$xe)
    
    #add estimate of variance
    if(is.null(object$fixed_est_var)){
      object$fixed_est_var <- rep(NA, ncol(object$fixed_est_vec))
      names(object$fixed_est_var) <- colnames(object$fixed_est_vec)
    }
    
    object$fixed_est_var <- compute.polyak.variance(pars,
                                                    polyak.rate,
                                                    i.start = burnin.end)
    
    if(plot){
      plot.trajectories2(x = pars, y = pars.smooth$x)
    }
    
  } else if(param == "random"){
    # n.fixed <- dim(object$mixedEffect_list$beta_fixed)[2]
    # n.random <- dim(object$mixedEffect_list$beta_random)[2]
    # names <- colnames(object$ranef_Sigma)
    # 
    # pars <- object$mixedEffect_list$betar_vec
    # colnames(pars) <- names
    # pars.smooth <- smooth.trajectory(pars,
    #                                  polyak.rate = polyak.rate,
    #                                  burnin.end = burnin.end)
    # object$mixedEffect_list$betar_vec <- pars.smooth$x
    # object$mixedEffect_list$beta_random <- pars.smooth$xe
    # 
    # #add to summaries (with fixed effects)
    # object$fixed_est_vec[,1:n.random] <- pars.smooth$x
    # object$fixed_est[1:n.random] <- pars.smooth$xe
    # 
    # #add estimate of variance
    # if(is.null(object$fixed_est_var)){
    #   n.fixed <- dim(object$mixedEffect_list$beta_fixed)[2]
    #   object$fixed_est_var <- rep(NA,n.random+n.fixed)
    #   names(object$fixed_est_var) <- names(object$fixed_est)
    # }
    # 
    # object$fixed_est_var[1:n.random] <- compute.polyak.variance(pars,polyak.rate,i.start = burnin.end)
    # pars.full <- pars
    # pars.full.smooth <- pars.smooth$x
    
    # process the covariance
    
    n.random <- length(object$mixedEffect_list$beta_random)
    
    ncol_Sigma <- dim(object$ranef_Sigma)[1]
    
    pars <- object$ranef_Sigma_vec
    #colnames(pars) <- rep("Sigma", dim(object$ranef_Sigma_vec)[2])
    
    pars.smooth <- smooth.trajectory(pars,
                                     polyak.rate = polyak.rate,
                                     burnin.end = burnin.end)
    
    object$ranef_Sigma_vec <- pars.smooth$x
    object$ranef_Sigma <- matrix(pars.smooth$xe, n.random, n.random, 
                                 dimnames = list(rownames(object$ranef_Sigma), 
                                                 colnames(object$ranef_Sigma)))
    
    object$mixedEffect_list$Sigma_vec <- pars.smooth$x
    object$mixedEffect_list$Sigma <- matrix(pars.smooth$xe, n.random, n.random,
                                            dimnames = list(rownames(object$ranef_Sigma), 
                                                            colnames(object$ranef_Sigma)))
    
    object$ranef_Sigma_var <- compute.polyak.variance(pars,
                                                      polyak.rate, 
                                                      i.start = burnin.end)
    
    if(plot){
      pars.full <- pars #cbind(pars.full, pars)
      pars.full.smooth <- pars.smooth$x #cbind(pars.full.smooth, pars.smooth$x)
      
      if(ncol_Sigma > 1){
        cols_omit         <- unlist(lapply(1:(ncol_Sigma - 1), 
                                           function(i) (ncol_Sigma * i + 1):(ncol_Sigma * i + i)))
        pars.full <- pars[, -cols_omit]
        pars.full.smooth <- pars.smooth$x[, -cols_omit]
        
        colnames <- colnames(summary(object)$random_results$Sigma_est)
        colnames_ext <- paste(rep(colnames, each = ncol_Sigma), colnames, sep = "_")
        colnames(pars.full) <- colnames(pars.full.smooth) <- colnames_ext[-cols_omit]
      }else{
        pars.full <- pars
        pars.full.smooth <- pars.smooth$x
        colnames(pars.full) <- colnames(pars.full.smooth) <- colnames(summary(object)$random_results$Sigma_est)
      }
      
    }
    
    if(object$random_distr %in% ("NIG")){
      
      mu <- object$ranef_mu_vec
      mu.smooth <- smooth.trajectory(mu, 
                                     polyak.rate = polyak.rate,
                                     burnin.end = burnin.end)
      
      object$ranef_mu_vec <- mu.smooth$x
      object$ranef_mu <- mu.smooth$xe
      
      object$mixedEffect_list$mu_vec <- mu.smooth$x
      object$mixedEffect_list$mu <- mu.smooth$xe
      
      object$ranef_mu_var <- compute.polyak.variance(mu, 
                                                     polyak.rate, 
                                                     i.start = burnin.end)
      
      nu <- object$ranef_nu_vec
      nu.smooth <- smooth.trajectory(nu,
                                     polyak.rate = polyak.rate,
                                     burnin.end = burnin.end)
      
      object$ranef_nu_vec <- nu.smooth$x
      object$ranef_nu <- nu.smooth$xe
      
      object$mixedEffect_list$nu_vec <- nu.smooth$x
      object$mixedEffect_list$nu <- nu.smooth$xe
      
      object$ranef_nu_var <- compute.polyak.variance(nu, 
                                                     polyak.rate, 
                                                     i.start = burnin.end)
      
      if(plot){
        pars.full <- cbind(pars.full, mu, nu)
        pars.full.smooth <- cbind(pars.full.smooth, mu.smooth$x, nu.smooth$x)
        
        name_index <- (ncol(pars.full) - ncol_Sigma ) : (ncol(pars.full))
        
        colnames(pars.full)[name_index] <- colnames(pars.full.smooth)[name_index] <- 
          c(paste("mu", colnames(object$ranef_Sigma), sep = "_"), "nu")
      }
    }
    if(plot){
      plot.trajectories2(pars.full, pars.full.smooth)
    }
    
  } else if(param == "process"){
    
    tau <- object$operator_tau_vec
    
    tau.smooth <- smooth.trajectory(tau,
                                    polyak.rate = polyak.rate,
                                    burnin.end = burnin.end)
    
    tau.smooth$xe <- as.numeric(tau.smooth$xe)
    
    object$operator_tau_vec <- tau.smooth$x
    object$operator_tau     <- tau.smooth$xe
    
    object$operator_list$tauVec <- tau.smooth$x
    object$operator_list$tau    <- tau.smooth$xe
    
    object$operator_tau_var <- compute.polyak.variance(tau,
                                                       polyak.rate,
                                                       i.start = burnin.end)
    
    pars <- matrix(tau, ncol = 1)
    pars.smooth <- matrix(tau.smooth$x, ncol = 1)
    colnames(pars) <- colnames(pars.smooth) <- "tau"
    
    if(object$operator_type == "matern" || object$operator_type == "exponential"){
      
      kappa <- object$operator_kappa_vec
      
      kappa.smooth <- smooth.trajectory(kappa,
                                        polyak.rate = polyak.rate,
                                        burnin.end = burnin.end)
      
      kappa.smooth$xe <- as.numeric(kappa.smooth$xe)
      
      object$operator_kappa_vec <- kappa.smooth$x
      object$operator_kappa     <- kappa.smooth$xe
      
      object$operator_list$kappaVec <- kappa.smooth$x
      object$operator_list$kappa    <- kappa.smooth$xe
      
      object$operator_kappa_var <- compute.polyak.variance(kappa, 
                                                           polyak.rate, 
                                                           i.start = burnin.end)
      pars <- cbind(pars, kappa)
      pars.smooth <- cbind(pars.smooth, kappa.smooth$x)
      colnames(pars)[2] <- colnames(pars.smooth)[2] <- "kappa"
      
    }
    
    if(object$process_distr %in% c("NIG", "GAL")){
      
      mu <- object$process_mu_vec
      
      mu.smooth <- smooth.trajectory(mu,
                                     polyak.rate = polyak.rate,
                                     burnin.end = burnin.end)
      
      mu.smooth$xe <- as.numeric(mu.smooth$xe)
      
      object$process_mu_vec <- mu.smooth$x
      object$process_mu     <- mu.smooth$xe
      
      object$processes_list$mu_vec <- mu.smooth$x
      object$processes_list$mu     <- mu.smooth$xe
      
      object$process_mu_var <- compute.polyak.variance(mu, 
                                                       polyak.rate,
                                                       i.start = burnin.end)
      
      nu <- object$process_nu_vec
      
      nu.smooth <- smooth.trajectory(nu,
                                     polyak.rate = polyak.rate,
                                     burnin.end = burnin.end)
      
      nu.smooth$xe <- as.numeric(nu.smooth$xe)
      
      object$process_nu_vec <- nu.smooth$x
      object$process_nu     <- nu.smooth$xe
      
      object$processes_list$nu_vec <- nu.smooth$x
      object$processes_list$nu     <- nu.smooth$xe
      
      object$process_nu_var <- compute.polyak.variance(nu, 
                                                       polyak.rate, 
                                                       i.start = burnin.end)
      
      pars <- cbind(pars, mu, nu)
      pars.smooth <- cbind(pars.smooth, mu.smooth$x, nu.smooth$x)
      
      colnames(pars)[(dim(pars)[2]-1):dim(pars)[2]] <- colnames(pars.smooth)[(dim(pars)[2]-1):dim(pars)[2]] <- c("mu", "nu")
      
    }
    if(plot){
      plot.trajectories2(pars, pars.smooth)
    }
    
  }else if(param == "error"){
    
    sigma <- object$meas_error_sigma_vec
    
    sigma.smooth <- smooth.trajectory(sigma, 
                                      polyak.rate = polyak.rate,
                                      burnin.end = burnin.end)
    
    sigma.smooth$xe <- as.numeric(sigma.smooth$xe)
    
    object$meas_error_sigma_vec <- sigma.smooth$x
    object$meas_error_sigma     <- sigma.smooth$xe
    
    object$measurementError_list$sigma_vec <- sigma.smooth$x
    object$measurementError_list$sigma     <- sigma.smooth$xe
    
    object$meas_error_sigma_var <- compute.polyak.variance(sigma,
                                                           polyak.rate,
                                                           i.start = burnin.end)
    
    pars <- matrix(sigma, ncol = 1)
    pars.smooth <- matrix(sigma.smooth$x, ncol = 1)
    
    colnames(pars) <- colnames(pars.smooth) <- "sigma"
    
    if(object$error_distr %in% c("NIG", "tdist")){
      
      nu <- object$meas_error_nu_vec
      
      nu.smooth <- smooth.trajectory(nu,
                                     polyak.rate = polyak.rate,
                                     burnin.end = burnin.end)
      
      nu.smooth$xe <- as.numeric(nu.smooth$xe)
      
      object$meas_error_nu_vec <- nu.smooth$x
      object$meas_error_nu     <- nu.smooth$xe
      
      object$measurementError_list$nu_vec <- nu.smooth$x
      object$measurementError_list$nu     <- nu.smooth$xe
      
      object$meas_error_nu_var <- compute.polyak.variance(nu,
                                                          polyak.rate,
                                                          i.start = burnin.end)
      
      pars <- cbind(pars, nu)
      pars.smooth <- cbind(pars.smooth, nu.smooth$x)
      
      colnames(pars)[2] <- colnames(pars.smooth)[2] <- "nu"
    }
    if(plot){
      plot.trajectories2(pars, pars.smooth)
    }
  }
  
  if(param == "fixed"){
    object$polyaked_fixed <- "yes"
  }else if(param == "random"){
    object$polyaked_random <- "yes"
  }else if(param == "process"){
    object$polyaked_process <- "yes"
  }else if(param == "error"){
    object$polyaked_error <- "yes"
  }
  
  return(object)
  
}



WLS.loss <- function(p,r){
  c <- exp(2*p[1])*exp(-exp(p[2])*r$lag[,1,1])
  w <- r$n.used - (1:(dim(r$acf)[1]))
  return(sum(w*(c-r$acf)^2))
}

compute.polyak.variance <- function(x,polyak.rate,i.start){
  
  if(length(dim(x))==2){
    n = dim(x)[2]
    n.iter = dim(x)[1]
    if(is.null(i.start)){
      i.start <- round(n.iter)/2 #base covariance esitmation on second half of the data
    }
    v <- x[1,]
    for(i in 1:n){
      r <- acf(x[i.start:n.iter,i],type="covariance",
               lag.max <- min(round(length(x[,i])/2),n.iter-i.start-1),
               plot=FALSE)
      p0 <- c(log(sqrt(var(x[i.start:n.iter,i]))), -log(max(r$lag[,1,1])))
      p <- optim(p0,WLS.loss, r = r)
      c <- exp(2*p$par[1])*exp(-exp(p$par[2])*r$lag[,1,1])
      v[i] <- (polyak.rate/(2-polyak.rate))*(c[1] + 2*sum(c[-1]))
    }
  } else {
    n.iter = length(x)
    if(is.null(i.start)){
      i.start <- round(n.iter)/2 #base covariance esitmation on second half of the data
    }
    r <- acf(x[i.start:n.iter],type="covariance",
             lag.max <- min(round(length(x)/2),n.iter-i.start-1),
             plot = FALSE)
    p0 <- c(log(sqrt(var(x[i.start:n.iter]))), -log(max(r$lag[,1,1])))
    p <- optim(p0,WLS.loss, r = r)
    c <- exp(2*p$par[1])*exp(-exp(p$par[2])*r$lag[,1,1])
    v <- (polyak.rate/(2-polyak.rate))*(c[1] + 2*sum(c[-1]))
    names(v) <- colnames(x)
  }
  
  return(v)
}



smooth.trajectory <- function(x,polyak.rate, burnin.end = NULL){
  if(is.null(burnin.end))
    burnin.end <- 2
  if(length(dim(x))==2){
    n = dim(x)[1]
    for(i in burnin.end:n){
      x[i,] <- polyak.rate * x[i,] + (1 - polyak.rate) * x[i-1,];
    }
    xe = matrix(x[n,],1,dim(x)[2])
  } else {
    n = length(x)
    for(i in burnin.end:n){
      x[i] <- polyak.rate * x[i] + (1 - polyak.rate) * x[i-1];
    }
    xe = matrix(x[n],1,1)
  }
  return(list(x=x,xe = xe))
}

plot.trajectories <- function(x, y){
  n <- 1
  if(length(dim(x)) == 2){
    n <- dim(x)[2]
  }
  m <- ceiling(sqrt(n))
  par(mfcol = c(m, m), mai = c(0, 0.2, 0.2, 0.1))
  for (i in 1:n){
    plot(x[, i], col = 1, type="l", main = dimnames(x)[[2]][i], xlab = "", ylab = "", xaxt = "n")
    lines(y[, i], col = 2, xlab = "", ylab = "")
  }
}

plot.trajectories2 <- function(x, y){
  
  if(ncol(x) == 1){
    ncol_plot <- 1
  }else if(ncol(x) < 5){
    ncol_plot <- 2
  }else{
    ncol_plot <- floor(ncol(x)/5) + 1
  }
  
  df <- data.frame(value = c(as.numeric(x), as.numeric(y)), 
                   var_name = rep(rep(colnames(x), each = nrow(x)), 2), 
                   series_name = rep(c("original", "polyak"), each = nrow(x) * ncol(x)),
                   iteration = rep(1:(nrow(x) * ncol(x)), 2))
  
  g <- ggplot(df, aes(x = iteration, y = value, group = series_name)) + 
    geom_line(aes(color = series_name)) + 
    facet_wrap(var_name ~ ., scales = "free", ncol = ncol_plot) + 
    xlab("Iteration") + ylab("Value") + 
    theme(legend.title=element_blank())
  
  print(g)
  
}
