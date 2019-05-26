
ngme.par <- function(n.cores, std.lim,max.rep,controls, controls.init = NULL,nIter,silent,...)
{
  output <- NULL
  step0 <- controls$step0
  if(is.null(step0))
    step0 <- 1
  
  alpha <- controls$alpha
  if(is.null(alpha))
    alpha = 0.6
  
  
  for(ii in 1:max.rep){
    cl <- makeCluster(n.cores)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = n.cores, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    parallel::clusterExport(cl, varlist = c(), envir = environment())
    if(ii>1){
      controls$iter.start <- (ii-1)*nIter
      controls$nBurnin = 5
    }
      
    est.list <- foreach(i = 1:n.cores, .options.snow = opts,.packages = c("ngme","Matrix")) %dopar%
    {
      if(ii==1){
        est <- ngme::ngme(controls=controls,
                          controls.init = controls.init,
                          nIter = nIter,silent = TRUE,...)
      } else {
        est <- ngme::ngme(init.fit = est.list.old[[i]],
                          controls=controls,
                          nIter = nIter,silent=TRUE,...)
      }
      return(est)  
    }
    close(pb)
    stopCluster(cl)
    est.list.old <- est.list
    est.merge <- merge.ngme.outputs(est.list)
    output <- attach.ngme.output(output,est.merge)
    
    n.fixed = dim(output$mixedEffect_list$betaf_vec)[2]
    n.random = dim(output$mixedEffect_list$betar_vec)[2]
    n.p = ceiling(sqrt(n.fixed+n.random))
    par(mfcol = c(n.p, n.p), mai = c(0.05, 0.3, 0.05, 0.05))
    for(k in 1:n.fixed){
      #fint plot limits
      min.y <- min(output$mixedEffect_list$betaf_vec[,k])
      max.y <- max(output$mixedEffect_list$betaf_vec[,k])
      for(i in 1:n.cores){
        min.i <- min(est.list[[i]]$mixedEffect_list$betaf_vec[,k])
        max.i <- max(est.list[[i]]$mixedEffect_list$betaf_vec[,k])
        min.y <- min(min.y,min.i)
        max.y <- max(max.y,max.i)
      }
        
      plot(output$mixedEffect_list$betaf_vec[,k],type="l", xlab="",ylab="",ylim=c(min.y,max.y),xaxt="n")  
      for(i in 1:n.cores)
        lines((ii-1)*nIter + (1:nIter),est.list[[i]]$mixedEffect_list$betaf_vec[,k],col="gray")  
      
      points((1:ii)*nIter,output$fixed_est_vec[,output$index_fixed[k]],col=2)
      points((1:ii)*nIter,output$fixed_est_vec[,output$index_fixed[k]]+2*sqrt(output$fixed_est_var[,output$index_fixed[k]]),col=3)
      points((1:ii)*nIter,output$fixed_est_vec[,output$index_fixed[k]]-2*sqrt(output$fixed_est_var[,output$index_fixed[k]]),col=3)  
    }
    for(k in 1:n.random){
      min.y <- min(output$mixedEffect_list$betar_vec[,k])
      max.y <- max(output$mixedEffect_list$betar_vec[,k])
      for(i in 1:n.cores){
        min.i <- min(est.list[[i]]$mixedEffect_list$betar_vec[,k])
        max.i <- max(est.list[[i]]$mixedEffect_list$betar_vec[,k])
        min.y <- min(min.y,min.i)
        max.y <- max(max.y,max.i)
      }
      plot(output$mixedEffect_list$betar_vec[,k],type="l",  xlab="",ylab="",ylim=c(min.y,max.y),xaxt="n")   
      for(i in 1:n.cores)
        lines((ii-1)*nIter + (1:nIter),est.list[[i]]$mixedEffect_list$betar_vec[,k],col="gray")  
      
      points((1:ii)*nIter,output$fixed_est_vec[,output$index_random[k]],col=2)
      points((1:ii)*nIter,output$fixed_est_vec[,output$index_random[k]]+2*sqrt(output$fixed_est_var[,output$index_random[k]]),col=3)
      points((1:ii)*nIter,output$fixed_est_vec[,output$index_random[k]]-2*sqrt(output$fixed_est_var[,output$index_random[k]]),col=3)  
    }
    
    if(0){
      plot(output$mixedEffect_list$Sigma_vec[,1],type="l", main = "Sigma")  
      for(i in 1:n.cores)
        lines((ii-1)*nIter + (1:nIter), est.list[[i]]$mixedEffect_list$Sigma_vec[,1],col="gray")  
      
      points((1:ii)*nIter,output$ranef_Sigma_vec[,1],col=2)
      points((1:ii)*nIter,output$ranef_Sigma_vec[,1]+2*sqrt(output$ranef_Sigma_var[,1]),col=3)
      points((1:ii)*nIter,output$ranef_Sigma_vec[,1]-2*sqrt(output$ranef_Sigma_var[,1]),col=3)
      if(!is.null(est.list[[1]]$operator_list)){
        plot(output$operator_list$tauVec,type="l", main = "tau")  
        for(i in 1:n.cores)
          lines((ii-1)*nIter + (1:nIter), est.list[[i]]$operator_list$tauVec,col="gray")  
        
        points((1:ii)*nIter,output$operator_tau_vec,col=2)
        points((1:ii)*nIter,output$operator_tau_vec+2*sqrt(output$operator_tau_var),col=3)
        points((1:ii)*nIter,output$operator_tau_vec-2*sqrt(output$operator_tau_var),col=3)
        
        if(!is.null(est.list[[1]]$processes_list$nu)){
          plot(output$processes_list$nu_vec,type="l", main = "nu process")    
          for(i in 1:n.cores)
            lines((ii-1)*nIter + (1:nIter), est.list[[i]]$processes_list$nu_vec,col="gray")
          
          points((1:ii)*nIter,output$process_nu_vec,col=2)
          points((1:ii)*nIter,output$process_nu_vec+2*sqrt(output$process_nu_var),col=3)
          points((1:ii)*nIter,output$process_nu_vec-2*sqrt(output$process_nu_var),col=3)
        }
      }  
    }
    
    
    
    converged <- check.convergence(output,std.lim)
    
    cat("fixed = ", converged$fixed==0,
        ", random =",converged$random==0,
        ", operator = ",converged$operator==0,
        ", process =",converged$process==0,"\n")
    if(converged$fixed>0){
      cat("fixed not converged:", converged$fixed.not.converged, "\n")
      cat("standard deviations:", converged$fixed.not.converged.std, "\n")
    }
      
    if(converged$fixed==0 && converged$random == 0 & converged$operator == 0 && converged$process == 0)
          break
  }
  return(output)
}

merge.ngme.outputs <- function(est.list){
  n.cores <- length(est.list)
  est.merge <- est.list[[1]]
  #merge mixedEffects
  if(n.cores>1){
    betaf_vec <- est.list[[1]]$mixedEffect_list$betaf_vec/n.cores
    beta_fixed_samples <- est.list[[1]]$mixedEffect_list$beta_fixed
    n.fixed <- dim( est.list[[1]]$mixedEffect_list$beta_fixed)[2]
    beta_f <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$mixedEffect_list$beta_fixed)),n.fixed,n.cores)
    beta_fixed <- apply(beta_f,1,mean)
    beta_fixed_var <- apply(beta_f,1,var)/n.cores
    
    use.random = FALSE
    if(!is.null(est.list[[1]]$mixedEffect_list$betar_vec)){
      betar_vec <- est.list[[1]]$mixedEffect_list$betar_vec/n.cores
      n.random <- dim( est.list[[1]]$mixedEffect_list$beta_random)[2]
      beta_r <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$mixedEffect_list$beta_random)),n.random,n.cores)
      beta_random <- apply(beta_r,1,mean)
      beta_random_var <- apply(beta_r,1,var)/n.cores
      beta_random_samples <- est.list[[1]]$mixedEffect_list$beta_random
      Sigma_vec <- est.list[[1]]$mixedEffect_list$Sigma_vec/n.cores
      sigma_v <- matrix(unlist(lapply(1:n.cores,function(x) c(est.list[[x]]$mixedEffect_list$Sigma))),n.random^2,n.cores)
      Sigma <- matrix(apply(sigma_v,1,mean),n.random,n.random)
      Sigma_var <- matrix(apply(sigma_v,1,var)/n.cores,n.random,n.random)
      use.random = TRUE
    }
    use.nu = FALSE
    if(!is.null(est.list[[1]]$mixedEffect_list$nu_vec)){
      nu_vec <- est.list[[1]]$mixedEffect_list$nu_vec/n.cores
      use.nu = TRUE
      nu_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$mixedEffect_list$nu)),1,n.cores)
      nu <- apply(nu_v,1,mean)
      nu_var <- apply(nu_v,1,var)/n.cores
    }
    use.mu = FALSE
    if(!is.null(est.list[[1]]$mixedEffect_list$mu_vec)){
      mu_vec <- est.list[[1]]$mixedEffect_list$mu_vec/n.cores
      mu <- est.list[[1]]$mixedEffect_list$mu/n.cores
      use.mu = TRUE
      mu_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$mixedEffect_list$mu)),1,n.cores)
      mu <- apply(mu_v,1,mean)
      mu_var <- apply(mu_v,1,var)/n.cores
    }
    
    for(i in 2:n.cores){
      betaf_vec <- betaf_vec + est.list[[i]]$mixedEffect_list$betaf_vec/n.cores
      beta_fixed_samples <- rbind(beta_fixed_samples,est.list[[i]]$mixedEffect_list$beta_fixed)
      if(use.random){
        beta_random_samples <- rbind(beta_random_samples,est.list[[i]]$mixedEffect_list$beta_random)
        betar_vec <- betar_vec + est.list[[i]]$mixedEffect_list$betar_vec/n.cores
        Sigma_vec <- Sigma_vec + est.list[[i]]$mixedEffect_list$Sigma_vec/n.cores
      }
      if(use.nu)
        nu_vec <- nu_vec + est.list[[i]]$mixedEffect_list$nu_vec/n.cores
      if(use.mu)
        mu_vec <- mu_vec + est.list[[i]]$mixedEffect_list$mu_vec/n.cores
    }  
    est.merge$mixedEffect_list$betaf_vec <- betaf_vec
    est.merge$mixedEffect_list$beta_fixed_samples <- list(beta_fixed_samples)
    est.merge$fixed_est_vec <- t(matrix(rep(0,length(est.merge$index_fixed)+length(est.merge$index_random))))
    est.merge$fixed_est_vec[est.merge$index_fixed] <- beta_fixed#betaf_vec
    est.merge$fixed_est[est.merge$index_fixed] <- beta_fixed
    est.merge$fixed_est_var <- t(matrix(rep(0,length(est.merge$index_fixed)+length(est.merge$index_random))))
    est.merge$fixed_est_var[est.merge$index_fixed] <- beta_fixed_var
    est.merge$mixedEffect_list$beta_fixed <- beta_fixed
    est.merge$mixedEffect_list$beta_fixed_var <- beta_fixed_var
    if(use.random){
      est.merge$mixedEffect_list$betar_vec <- betar_vec
      est.merge$mixedEffect_list$beta_random <- beta_random
      est.merge$mixedEffect_list$beta_random_samples <- list(beta_random_samples)
      est.merge$mixedEffect_list$beta_random_var <- beta_random_var
      est.merge$fixed_est_vec[est.merge$index_random] <- beta_random#betar_vec
      est.merge$fixed_est[est.merge$index_random] <- beta_random
      est.merge$fixed_est_var[est.merge$index_random] <- beta_random_var
      est.merge$mixedEffect_list$Sigma_vec <- Sigma_vec
      est.merge$mixedEffect_list$Sigma <- Sigma
      est.merge$mixedEffect_list$Sigma_var <- Sigma_var
      est.merge$ranef_Sigma <- Sigma
      est.merge$ranef_Sigma_var <- Sigma_var
      est.merge$ranef_Sigma_vec <- t(matrix(c(Sigma)))#Sigma_vec
    }
    if(use.nu){
      est.merge$mixedEffect_list$nu_vec <- nu_vec 
      est.merge$mixedEffect_list$nu <- nu
      est.merge$ranef_nu_vec <- nu#nu_vec 
      est.merge$ranef_nu <- nu
      est.merge$ranef_nu_var <- nu_var
      est.merge$mixedEffect_list$nu_var <- nu_var
    }
    if(use.mu){
      est.merge$mixedEffect_list$mu_vec <- mu_vec 
      est.merge$mixedEffect_list$mu <- mu 
      est.merge$mixedEffect_list$mu_var <- mu_var 
      est.merge$ranef_mu_vec <- mu#mu_vec 
      est.merge$ranef_mu <- mu
      est.merge$ranef_mu_var <- mu_var
    }
    #merge error list
    sigma_vec <- est.list[[1]]$mixedEffect_list$betaf_vec/n.cores
    sigma_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$measurementError_list$sigma)),1,n.cores)
    sigma <- apply(sigma_v,1,mean)
    sigma_var <- apply(sigma_v,1,var)/n.cores
    
    use.nu = FALSE
    if(!is.null(est.list[[1]]$measurementError_list$nu_vec)){
      nu_vec <- est.list[[1]]$measurementError_list$nu_vec/n.cores
      use.nu = TRUE
      nu_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$measurementError_list$nu)),1,n.cores)
      nu <- apply(nu_v,1,mean)
      nu_var <- apply(nu_v,1,var)/n.cores
    }
    use.mu = FALSE
    if(!is.null(est.list[[1]]$measurementError_list$mu_vec)){
      mu_vec <- est.list[[1]]$measurementError_list$mu_vec/n.cores
      mu <- est.list[[1]]$measurementError_list$mu/n.cores
      use.mu = TRUE
      mu_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$measurementError_list$mu)),1,n.cores)
      mu <- apply(mu_v,1,mean)
      mu_var <- apply(mu_v,1,var)/n.cores
    }
    
    for(i in 2:n.cores){
      sigma_vec <- sigma_vec + est.list[[i]]$measurementError_list$sigma_vec/n.cores
      if(use.nu)
        nu_vec <- nu_vec + est.list[[i]]$measurementError_list$nu_vec/n.cores
      if(use.mu)
        mu_vec <- mu_vec + est.list[[i]]$measurementError_list$mu_vec/n.cores
    }  
    est.merge$measurementError_list$sigma <- sigma
    est.merge$measurementError_list$sigma_var <- sigma_var
    est.merge$meas_error_sigma <- sigma
    est.merge$meas_error_sigma_var <- sigma_var
    est.merge$meas_error_sigma_vec <- sigma#sigma_vec
    if(use.nu){
      est.merge$measurementError_list$nu_vec <- nu_vec 
      est.merge$measurementError_list$nu <- nu 
      est.merge$measurementError_list$nu_var <- nu_var 
      est.merge$meas_error_nu <- nu
      est.merge$meas_error_nu_var <- nu_var
      est.merge$meas_error_nu_vec <- nu#nu_vec
    }
    if(use.mu){
      est.merge$measurementError_list$mu_vec <- mu_vec 
      est.merge$measurementError_list$mu <- mu 
      est.merge$measurementError_list$mu_var <- mu_var 
      est.merge$meas_error_mu <- mu
      est.merge$meas_error_mu_var <- mu_var
      est.merge$meas_error_mu_vec <- mu#mu_vec
    }
    #merge operator list
    if(!is.null(est.list[[1]]$operator_list)){
      tau_vec <- est.list[[1]]$operator_list$tauVec/n.cores
      tau_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$operator_list$tau)),1,n.cores)
      tau <- apply(tau_v,1,mean)
      tau_var <- apply(tau_v,1,var)/n.cores  
      use.kappa = FALSE
      if(!is.null(est.list[[1]]$operator_list$kappa)){
        use.kappa = TRUE
        kappa_vec <- est.list[[1]]$operator_list$kappaVec/n.cores
        kappa_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$operator_list$kappa)),1,n.cores)
        kappa <- apply(kappa_v,1,mean)
        kappa_var <- apply(kappa_v,1,var)/n.cores  
      }
      for(i in 2:n.cores){
        tau_vec <- tau_vec + est.list[[i]]$operator_list$tauVec/n.cores
        if(use.kappa)
          kappa_vec <- kappa_vec + est.list[[i]]$operator_list$kappaVec/n.cores
      }  
      est.merge$operator_list$tau <- tau
      est.merge$operator_list$tauVec <- tau_vec
      est.merge$operator_list$tau_var <- tau_var
      est.merge$operator_tau <- tau
      est.merge$operator_tau_var <- tau_var
      est.merge$operator_tau_vec <- tau#tau_vec
      if(use.kappa){
        est.merge$operator_list$kappa <- kappa
        est.merge$operator_list$kappaVec <- tau_vec
        est.merge$operator_list$kappa_var <- kappa_var
        est.merge$operator_kappa <- kappa
        est.merge$operator_kappa_var <- kappa_var
        est.merge$operator_kappa_vec <- kappa#kappa_vec
      }
      #merge process list
      use.nu = FALSE
      if(!is.null(est.list[[1]]$processes_list$nu_vec)){
        nu_vec <- est.list[[1]]$processes_list$nu_vec/n.cores
        use.nu = TRUE
        nu_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$processes_list$nu)),1,n.cores)
        nu <- apply(nu_v,1,mean)
        nu_var <- apply(nu_v,1,var)/n.cores
      }
      use.mu = FALSE
      if(!is.null(est.list[[1]]$processes_list$mu_vec)){
        mu_vec <- est.list[[1]]$processes_list$mu_vec/n.cores
        mu <- est.list[[1]]$processes_list$mu/n.cores
        use.mu = TRUE
        mu_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$processes_list$mu)),1,n.cores)
        mu <- apply(mu_v,1,mean)
        mu_var <- apply(mu_v,1,var)/n.cores
      }
      
      for(i in 2:n.cores){
        if(use.nu)
          nu_vec <- nu_vec + est.list[[i]]$processes_list$nu_vec/n.cores
        if(use.mu)
          mu_vec <- mu_vec + est.list[[i]]$processes_list$mu_vec/n.cores
      }  
      if(use.nu){
        est.merge$processes_list$nu_vec <- nu_vec 
        est.merge$processes_list$nu <- nu 
        est.merge$processes_list$nu_var <- nu_var 
        est.merge$process_nu <- nu
        est.merge$process_nu_var <- nu_var
        est.merge$process_nu_vec <- nu#nu_vec
      }
      if(use.mu){
        est.merge$processes_list$mu_vec <- mu_vec 
        est.merge$processes_list$mu <- mu 
        est.merge$processes_list$mu_var <- mu_var 
        est.merge$process_mu <- mu
        est.merge$process_mu_var <- mu_var
        est.merge$process_mu_vec <- mu#mu_vec
      }
    }
  }
  return(est.merge)
}

attach.ngme.output <- function(obj1,obj2){
  if(is.null(obj1)){
    return(obj2)
  } else {
    output <- obj2
    #merge mixedEffects
    output$mixedEffect_list$betaf_vec <- rbind(obj1$mixedEffect_list$betaf_vec,
                                               obj2$mixedEffect_list$betaf_vec)
    l <- length(obj1$mixedEffect_list$beta_fixed_samples)
    output$mixedEffect_list$beta_fixed_samples <- obj1$mixedEffect_list$beta_fixed_samples
    output$mixedEffect_list$beta_fixed_samples[[l+1]] <- obj2$mixedEffect_list$beta_fixed_samples[[1]]
    n <- dim(output$mixedEffect_list$betaf_vec)[1]
    #output$fixed_est_vec <- matrix(0,n,length(obj2$index_fixed)+length(obj2$index_random))
    #output$fixed_est_vec[,obj2$index_fixed] <- output$mixedEffect_list$betaf_vec
    
    output$fixed_est_vec <- rbind(obj1$fixed_est_vec,obj2$fixed_est_vec)
    output$fixed_est_var <- rbind(obj1$fixed_est_var,obj2$fixed_est_var)
    
    if(!is.null(obj2$mixedEffect_list$betar_vec)){
      
      output$mixedEffect_list$beta_random_samples <- obj1$mixedEffect_list$beta_random_samples
      output$mixedEffect_list$beta_random_samples[[l+1]] <- obj2$mixedEffect_list$beta_random_samples[[1]]
      output$mixedEffect_list$betar_vec <- rbind(obj1$mixedEffect_list$betar_vec,
                                                 obj2$mixedEffect_list$betar_vec)
      #output$fixed_est_vec[,obj2$index_random] <- output$mixedEffect_list$betar_vec
      output$mixedEffect_list$Sigma_vec <- rbind(obj1$mixedEffect_list$Sigma_vec,
                                                 obj2$mixedEffect_list$Sigma_vec)
      #output$ranef_Sigma_vec <- output$mixedEffect_list$Sigma_vec
      output$ranef_Sigma_vec <- rbind(obj1$ranef_Sigma_vec,obj2$ranef_Sigma_vec)
      output$ranef_Sigma_var <- rbind(obj1$ranef_Sigma_var,obj2$ranef_Sigma_var)
      
    }
    
    if(!is.null(obj2$mixedEffect_list$nu_vec)){
      output$mixedEffect_list$nu_vec <- rbind(obj1$mixedEffect_list$nu_vec,
                                              obj2$mixedEffect_list$nu_vec)
      #output$ranef_nu_vec <- output$mixedEffect_list$nu_vec
      output$ranef_nu_vec <- rbind(obj1$ranef_nu_vec,obj2$ranef_nu_vec)
      output$ranef_nu_var <- rbind(obj1$ranef_nu_var,obj2$ranef_nu_var)
    }
    
    if(!is.null(obj2$mixedEffect_list$mu_vec)){
      output$mixedEffect_list$mu_vec <- rbind(obj1$mixedEffect_list$mu_vec,
                                              obj2$mixedEffect_list$mu_vec)
      #output$ranef_mu_vec <- output$mixedEffect_list$mu_vec
      output$ranef_mu_vec <- rbind(obj1$ranef_mu_vec,obj2$ranef_mu_vec)
      output$ranef_mu_var <- rbind(obj1$ranef_mu_var,obj2$ranef_mu_var)
    }
    #merge error list
    output$measurementError_list$sigma_vec <- rbind(obj1$measurementError_list$sigma_vec,
                                                    obj2$measurementError_list$sigma_vec)
    #output$meas_error_sigma_vec <- output$measurementError_list$sigma_vec
    output$meas_error_sigma_vec <- rbind(obj1$meas_error_sigma_vec,obj2$meas_error_sigma_vec)
    output$meas_error_sigma_var <- rbind(obj1$meas_error_sigma_var,obj2$meas_error_sigma_var)
    
    if(!is.null(obj2$measurementError_list$nu_vec)){
      output$measurementError_list$nu_vec <- rbind(obj1$measurementError_list$nu_vec,
                                                   obj2$measurementError_list$nu_vec)
      #output$meas_error_nu_vec <- output$measurementError_list$nu_vec
      output$meas_error_nu_vec <- rbind(obj1$meas_error_nu_vec,obj2$meas_error_nu_vec)
      output$meas_error_nu_var <- rbind(obj1$meas_error_nu_var,obj2$meas_error_nu_var)
      
    }
    if(!is.null(obj2$measurementError_list$mu_vec)){
      output$measurementError_list$mu_vec <- rbind(obj1$measurementError_list$mu_vec,
                                                   obj2$measurementError_list$mu_vec)
      #output$meas_error_mu_vec <- output$measurementError_list$mu_vec
      output$meas_error_mu_vec <- rbind(obj1$meas_error_mu_vec,obj2$meas_error_mu_vec)
      output$meas_error_mu_var <- rbind(obj1$meas_error_mu_var,obj2$meas_error_mu_var)
    }
    
    #merge operator list
    if(!is.null(obj2$operator_list$tauVec)){
      output$operator_list$tauVec <- rbind(matrix(obj1$operator_list$tauVec),
                                           matrix(obj2$operator_list$tauVec))
      #output$operator_tau_vec <- output$operator_list$tau_vec
      output$operator_tau_vec <- rbind(obj1$operator_tau_vec,obj2$operator_tau_vec)
      output$operator_tau_var <- rbind(obj1$operator_tau_var,obj2$operator_tau_var)  
    }
    
    
    if(!is.null(obj2$operator_list$kappaVec)){
      output$operator_list$kappa_vec <- rbind(matrix(obj1$operator_list$kappaVec),
                                                   matrix(obj2$operator_list$kappaVec))
      #output$operator_kappa_vec <- output$operator_list$kappa_vec
      output$operator_kappa_vec <- rbind(obj1$operator_kappa_vec,obj2$operator_kappa_vec)
      output$operator_kappa_var <- rbind(obj1$operator_kappa_var,obj2$operator_kappa_var)
    }
    
    #merge process list 
    if(!is.null(obj2$processes_list$mu_vec)){
      output$processes_list$mu_vec <- rbind(obj1$processes_list$mu_vec,
                                            obj2$processes_list$mu_vec)
      #output$process_mu_vec <- output$processes_list$mu_vec
      output$process_mu_vec <- rbind(obj1$process_mu_vec,obj2$process_mu_vec)
      output$process_mu_var <- rbind(obj1$process_mu_var,obj2$process_mu_var)
    }
    if(!is.null(obj2$processes_list$nu_vec)){
      output$processes_list$nu_vec <- rbind(obj1$processes_list$nu_vec,
                                            obj2$processes_list$nu_vec)
      #output$process_nu_vec <- output$processes_list$nu_vec
      output$process_nu_vec <- rbind(obj1$process_nu_vec,obj2$process_nu_vec)
      output$process_nu_var <- rbind(obj1$process_nu_var,obj2$process_nu_var)
    }
    return(output)
  }
   
}

check.convergence <- function(output,std.lim, verbose = TRUE){
  
  N  <- nrow(output$fixed_est_var)
  #check mixed effects
  ind.fixed <- output$fixed_est/sqrt(output$fixed_est_var[N,])
  
  if(N>2){
    n.fixed <- length(output$fixed_est)
    significant <- rep(0,n.fixed)
    for(i in 1:n.fixed){
      significant[i] <- simple.convergence.test(output$fixed_est_vec[(N-2):N,i],output$fixed_est_var[(N-2):N,i])
    }  
    cat("significant: ",significant,"\n")
  }
  
  
  ind.ranef_Sigma <- NULL
  if(!is.null(output$ranef_Sigma))
    ind.ranef_Sigma <- c(output$ranef_Sigma)/sqrt(c(output$ranef_Sigma_var[N,]))

  ind.ranef_nu <- NULL
  if(!is.null(output$ranef_nu) && !is.na(output$ranef_nu))
    ind.ranef_nu <- output$ranef_nu/sqrt(output$ranef_nu_var[N])
  
  ind.ranef_mu <- NULL
  if(!is.null(output$ranef_mu) && !is.na(output$ranef_mu))
    ind.ranef_mu <- output$ranef_mu/sqrt(output$ranef_mu_var[N])
  
  #check error
  ind.meas_error_sigma <- output$meas_error_sigma/sqrt(output$meas_error_sigma_var[N])
  ind.meas_error_nu <- NULL
  if(!is.null(output$meas_error_nu) && !is.na(output$meas_error_nu))
    ind.meas_error_nu <- output$meas_error_nu/sqrt(output$meas_error_nu_var[N]) 
  
  ind.meas_error_mu <- NULL
  if(!is.null(output$meas_error_mu) && !is.na(output$meas_error_mu))
    ind.meas_error_mu <- output$meas_error_mu/sqrt(output$meas_error_mu_var[N]) 
  
  #check operator
  ind.operator_tau <- NULL
  if(!is.null(output$operator_tau) && !is.na(output$operator_tau))
    ind.operator_tau <- output$operator_tau/sqrt(output$operator_tau_var[N]) 
  
  ind.operator_kappa <- NULL
  if(!is.null(output$operator_kappa) && !is.na(output$operator_kappa))
    ind.operator_kappa <- output$operator_kappa/sqrt(output$operator_kappa_var[N]) 
  
  #check process
  ind.process_nu <- NULL
  if(!is.null(output$process_nu) && !is.na(output$process_nu))
    ind.process_nu <- output$process_nu/sqrt(output$process_nu_var[N]) 
  
  ind.process_mu <- NULL
  if(!is.null(output$process_mu) && !is.na(output$process_mu))
    ind.process_mu <- output$process_mu/sqrt(output$process_mu_var[N]) 
  
  converged.fixed <- sum(abs(ind.fixed)<std.lim)
  if(verbose){
    cat("fixed = ", abs(ind.fixed), "\n")
  }
  converged.random <- 0
  if(!is.null(ind.ranef_Sigma)){
    converged.random <- sum(abs(ind.ranef_Sigma)<std.lim)
    if(verbose){
      cat("Sigma = ", abs(ind.ranef_Sigma), "\n")
    }
    if(!is.null(ind.ranef_mu)){
      converged.random <- converged.random + sum(abs(ind.ranef_mu)<std.lim)
      if(verbose){
        cat("mu_mixed = ", abs(ind.ranef_mu), "\n")
      }
    }
      
    if(!is.null(ind.ranef_nu)){
      converged.random <- converged.random + sum(abs(ind.ranef_nu)<std.lim)
      if(verbose){
        cat("nu_mixed = ", abs(ind.ranef_mu), "\n")
      }
    }
  }
  converged.operator <- 0
  if(!is.null(ind.operator_tau)){
    converged.operator <- sum(abs(ind.operator_tau)<std.lim) 
    if(verbose){
      cat("tau = ", abs(ind.operator_tau), "\n")
    }
  }
    
  if(!is.null(ind.operator_kappa)){
    converged.operator <- converged.operator + sum(abs(ind.operator_kappa)<std.lim)
    if(verbose){
      cat("kappa = ", abs(ind.operator_kappa), "\n")
    }
  }
    
  converged.process <-0
  if(!is.null(ind.process_nu)){
    converged.process <- sum(abs(ind.process_nu)<std.lim) 
    if(verbose){
      cat("nu_process = ", abs(ind.process_nu), "\n")
    }
  }
    
  if(!is.null(ind.process_mu)){
    converged.process <- converged.process + sum(abs(ind.process_mu)<std.lim) 
    if(verbose){
      cat("mu_process = ", abs(ind.process_mu), "\n")
    }
  }
    
  
  return(list(fixed = converged.fixed,
              random = converged.random,
              operator = converged.operator,
              process = converged.process,
              fixed.not.converged = output$fixed_est[ind.fixed>std.lim],
              fixed.not.converged.std = sqrt(output$fixed_est_var[N,ind.fixed>std.lim])))
  
}

simple.convergence.test <- function(m,sigma2){
  B <- cbind(c(1,1,1),c(1,2,3))
  Sigma <- diag(sigma2)
  Q <- solve(t(B)%*%solve(Sigma,B))
  beta <- Q%*%(t(B)%*%solve(Sigma,m))
  s2 <- Q[2,2]
  return(abs(beta[2])>2*sqrt(s2))
}
