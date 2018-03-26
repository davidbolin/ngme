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
#' @param param plot Plot the results or not?
#'
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
                        polyak.rate = 0,
                        plot = TRUE,
                        burnin.end = NULL){

  # if process was not specified in the model
  if(param == "process" && object$use_process == FALSE){
    stop("No process was included in the model fit")
  }

  if(param == "fixed"){
    n.random <- dim(object$mixedEffect_list$beta_random)[2]
    n.fixed <- dim(object$mixedEffect_list$beta_fixed)[2]
    pars <- object$fixed_est_vec[,(n.random+1):(n.random+n.fixed)]
    pars.smooth <- smooth.trajectory(pars,polyak.rate = polyak.rate)
    object$mixedEffect_list$betaf_vec <- pars.smooth$x
    object$mixedEffect_list$beta_fixed <- pars.smooth$xe

    #update summaries (which also contain the random effects)
    object$fixed_est_vec[,(n.random+1):(n.random+n.fixed)] <- pars.smooth$x
    object$fixed_est[(n.random+1):(n.random+n.fixed)] <- pars.smooth$xe

    #add estimate of variance
    if(is.null(object$fixed_est_var)){
      object$fixed_est_var <- rep(NA,n.random+n.fixed)
      names(object$fixed_est_var) <- names(object$fixed_est)
    }

    object$fixed_est_var[(n.random+1):(n.random+n.fixed)] <- compute.polyak.variance(pars,
                                                                                     polyak.rate,
                                                                                     i.start = burnin.end)

    if(plot){
      plot.trajectories(pars,pars.smooth$x)
    }
  } else if(param == "random"){
    n.fixed <- dim(object$mixedEffect_list$beta_fixed)[2]
    n.random <- dim(object$mixedEffect_list$beta_random)[2]
    names <- names(object$fixed_est)[1:n.random]

    pars <- object$mixedEffect_list$betar_vec
    colnames(pars) <- names
    pars.smooth <- smooth.trajectory(pars,polyak.rate = polyak.rate)
    object$mixedEffect_list$betar_vec <- pars.smooth$x
    object$mixedEffect_list$beta_random <- pars.smooth$xe

    #add to summaries (with fixed effects)
    object$fixed_est_vec[,1:n.random] <- pars.smooth$x
    object$fixed_est[1:n.random] <- pars.smooth$xe

    #add estimate of variance
    if(is.null(object$fixed_est_var)){
      n.fixed <- dim(object$mixedEffect_list$beta_fixed)[2]
      object$fixed_est_var <- rep(NA,n.random+n.fixed)
      names(object$fixed_est_var) <- names(object$fixed_est)
    }

    object$fixed_est_var[1:n.random] <- compute.polyak.variance(pars,polyak.rate,i.start = burnin.end)
    pars.full <- pars
    pars.full.smooth <- pars.smooth$x

    # process the covariance
    ncol_Sigma        <- dim(object$ranef_Sigma)[1]
    pars <- object$ranef_Sigma_vec
    colnames(pars) <- rep("Sigma",dim(object$ranef_Sigma_vec)[2])

    pars.smooth <- smooth.trajectory(pars,polyak.rate = polyak.rate)
    object$ranef_Sigma_vec <- pars.smooth$x
    object$ranef_Sigma <- pars.smooth$xe
    object$mixedEffect_list$Sigma <- pars.smooth$xe
    object$ranef_Sigma_var <- compute.polyak.variance(pars,polyak.rate,i.start = burnin.end)

    if(plot){
      pars.full <- cbind(pars.full,pars)
      pars.full.smooth <- cbind(pars.full.smooth,pars.smooth$x)
    }

    if(object$random_distr %in% ("NIG")){
      mu <- object$ranef_mu_vec
      mu.smooth <- smooth.trajectory(mu,polyak.rate = polyak.rate)
      object$ranef_mu_vec <- mu.smooth$x
      object$ranef_mu <- mu.smooth$xe
      object$processes_list$mu <- mu.smooth$xe
      object$ranef_mu_var <- compute.polyak.variance(mu,polyak.rate,i.start = burnin.end)

      nu <- object$ranef_nu_vec
      nu.smooth <- smooth.trajectory(nu,polyak.rate = polyak.rate)
      object$ranef_nu_vec <- nu.smooth$x
      object$ranef_nu <- nu.smooth$xe
      object$processes_list$nu <- nu.smooth$xe
      object$ranef_nu_var <- compute.polyak.variance(nu,polyak.rate,i.start = burnin.end)

      if(plot){
        pars.full <- cbind(pars.full, mu, nu)
        pars.full.smooth <- cbind(pars.full.smooth,mu.smooth$x, nu.smooth$x)
        colnames(pars.full)[(ncol(pars.full)-ncol_Sigma ) : (ncol(pars.full))] <-
          c(paste("mu", colnames(object$ranef_Sigma), sep = "_"), "nu")
      }
    }
    if(plot){
      plot.trajectories(pars.full,pars.full.smooth)
    }
  } else if(param == "process"){
    tau <- object$operator_tau_vec
    tau.smooth <- smooth.trajectory(tau,polyak.rate = polyak.rate)
    object$operator_tau_vec <-tau.smooth$x
    object$operator_tau <- tau.smooth$xe
    object$operator_list$tau <- tau.smooth$xe
    object$operator_tau_var <- compute.polyak.variance(tau,polyak.rate,i.start = burnin.end)
    pars <- matrix(tau, ncol = 1)
    pars.smooth <- matrix(tau.smooth$x, ncol = 1)
    colnames(pars) <- "tau"
    if(object$operator_type == "matern"){
      kappa <- object$operator_kappa_vec
      kappa.smooth <- smooth.trajectory(kappa,polyak.rate = polyak.rate)
      object$operator_kappa_vec <-kappa.smooth$x
      object$operator_kappa <- kappa.smooth$xe
      object$operator_list$kappa <- kappa.smooth$xe
      object$operator_kappa_var <- compute.polyak.variance(kappa,polyak.rate,i.start = burnin.end)
      pars <- cbind(pars,kappa)
      pars.smooth <- cbind(pars.smooth,kappa.smooth$x)
      colnames(pars) <- c("tau", "kappa")
    }

    if(object$process_distr %in% c("NIG", "GAL")){
      mu <- object$process_mu_vec
      mu.smooth <- smooth.trajectory(mu,polyak.rate = polyak.rate)
      object$process_mu_vec <-mu.smooth$x
      object$process_mu <- mu.smooth$xe
      object$processes_list$mu <- mu.smooth$xe
      object$process_mu_var <- compute.polyak.variance(mu,polyak.rate,i.start = burnin.end)

      nu <- object$process_nu_vec
      nu.smooth <- smooth.trajectory(nu,polyak.rate = polyak.rate)
      object$process_nu_vec <-nu.smooth$x
      object$process_nu <- nu.smooth$xe
      object$processes_list$nu <- nu.smooth$xe
      object$process_nu_var <- compute.polyak.variance(nu,polyak.rate,i.start = burnin.end)
      pars <- cbind(pars, mu,nu)
      pars.smooth <- cbind(pars.smooth,mu.smooth$x,nu.smooth$x)
    }
    if(plot){
      plot.trajectories(pars,pars.smooth)
    }
    } else if(param == "error"){

      sigma <- object$meas_error_sigma_vec
      sigma.smooth <- smooth.trajectory(sigma,polyak.rate = polyak.rate)
      object$meas_error_sigma_vec <-sigma.smooth$x
      object$meas_error_sigma <- sigma.smooth$xe
      object$measurementError_list$sigma <- sigma.smooth$xe
      object$meas_error_sigma_var <- compute.polyak.variance(sigma,polyak.rate,i.start = burnin.end)
      pars <- matrix(sigma, ncol = 1)
      pars.smooth <- matrix(sigma.smooth$x, ncol = 1)
      colnames(pars) <- "sigma"

      if(object$error_distr %in% c("NIG", "tdist")){
        nu <- object$meas_error_nu_vec
        nu.smooth <- smooth.trajectory(nu,polyak.rate = polyak.rate)
        object$meas_error_nu_vec <-nu.smooth$x
        object$meas_error_nu <- nu.smooth$xe
        object$measurementError_list$nu <- nu.smooth$xe
        object$meas_error_nu_var <- compute.polyak.variance(nu,polyak.rate,i.start = burnin.end)
        pars <- cbind(pars, nu)
        pars.smooth <- cbind(pars.smooth,nu.smooth$x)
        colnames(pars)[2] <- c("nu")
      }
      if(plot){
        plot.trajectories(pars,pars.smooth)
      }
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



smooth.trajectory <- function(x,polyak.rate){
  if(length(dim(x))==2){
    n = dim(x)[1]
    for(i in 2:n){
      x[i,] <- polyak.rate * x[i,] + (1 - polyak.rate) * x[i-1,];
    }
    xe = matrix(x[n,],1,dim(x)[2])
  } else {
    n = length(x)
    for(i in 2:n){
      x[i] <- polyak.rate * x[i] + (1 - polyak.rate) * x[i-1];
    }
    xe = matrix(x[n],1,1)
  }
  return(list(x=x,xe = xe))
}

plot.trajectories <- function(x,y){
  n = 1
  if(length(dim(x))==2){
    n = dim(x)[2]
  }
  m <- ceiling(sqrt(n))
  par(mfcol=c(m,m),mai = c(0,0.2,0.2,0.1))
  for (i in 1:n) {
    plot(x[,i],col=1,type="l",main = dimnames(x)[[2]][i],xlab="",ylab="",xaxt="n")
    lines(y[,i],col=2,xlab="",ylab="")
  }
}

