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
    if(use.random){
      est.merge$fixed_est_vec <- t(matrix(rep(0,length(est.merge$index_fixed)+length(est.merge$index_random))))  
      est.merge$fixed_est_var <- t(matrix(rep(0,length(est.merge$index_fixed)+length(est.merge$index_random))))
    } else {
      est.merge$fixed_est_vec <- t(matrix(rep(0,length(est.merge$index_fixed))))
      est.merge$fixed_est_var <- t(matrix(rep(0,length(est.merge$index_fixed))))
    }
    est.merge$fixed_est_vec[est.merge$index_fixed] <- beta_fixed#betaf_vec
    est.merge$fixed_est[est.merge$index_fixed] <- beta_fixed
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
      est.merge$mixedEffect_list$nu_vec <- matrix(nu_vec)
      est.merge$mixedEffect_list$nu <- nu
      est.merge$ranef_nu_vec <- nu#nu_vec 
      est.merge$ranef_nu <- nu
      est.merge$ranef_nu_var <- nu_var
      est.merge$mixedEffect_list$nu_var <- nu_var
    }
    if(use.mu){
      est.merge$mixedEffect_list$mu_vec <- matrix(mu_vec )
      est.merge$mixedEffect_list$mu <- mu 
      est.merge$mixedEffect_list$mu_var <- mu_var 
      est.merge$ranef_mu_vec <- mu#mu_vec 
      est.merge$ranef_mu <- mu
      est.merge$ranef_mu_var <- mu_var
    }
    #merge error list
    bivariate = FALSE
    if(est.list[[1]]$measurementError_list$noise == "nsNormal"){
      bivariate = TRUE
    }
    if(bivariate){ 
      theta_vec <- est.list[[1]]$measurementError_list$theta_vec/n.cores
      theta_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$measurementError_list$theta)),2,n.cores)
      theta <- apply(theta_v,1,mean)
      theta_var <- apply(theta_v,1,var)/n.cores  
    } else {
      sigma_vec <- est.list[[1]]$measurementError_list$sigma_vec/n.cores
      sigma_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$measurementError_list$sigma)),1,n.cores)
      sigma <- apply(sigma_v,1,mean)
      sigma_var <- apply(sigma_v,1,var)/n.cores  
    }
    
    
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
      if(bivariate){
        theta_vec <- theta_vec + est.list[[i]]$measurementError_list$theta_vec/n.cores
      } else {
        sigma_vec <- sigma_vec + est.list[[i]]$measurementError_list$sigma_vec/n.cores  
      }
      
      if(use.nu)
        nu_vec <- nu_vec + est.list[[i]]$measurementError_list$nu_vec/n.cores
      if(use.mu)
        mu_vec <- mu_vec + est.list[[i]]$measurementError_list$mu_vec/n.cores
    }  
    if(bivariate){
      est.merge$measurementError_list$theta <- theta
      est.merge$measurementError_list$theta_vec <- theta_vec
      est.merge$measurementError_list$theta_var <- theta_var
      est.merge$meas_error_sigma <- exp(theta)
      est.merge$meas_error_sigma_vec <- matrix(theta,1,2)
      est.merge$meas_error_sigma_var <- matrix(theta_var,1,2)
    } else {
      est.merge$measurementError_list$sigma <- sigma
      est.merge$measurementError_list$sigma_vec <- sigma_vec
      est.merge$measurementError_list$sigma_var <- sigma_var
      est.merge$meas_error_sigma <- sigma
      est.merge$meas_error_sigma_var <- sigma_var
      est.merge$meas_error_sigma_vec <- sigma#sigma_vec  
    }
    
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
      if(bivariate){
        tau1_vec <- est.list[[1]]$operator_list$tau1Vec/n.cores
        tau1_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$operator_list$tau1)),1,n.cores)
        tau1 <- apply(tau1_v,1,mean)
        tau1_var <- apply(tau1_v,1,var)/n.cores  

        tau2_vec <- est.list[[1]]$operator_list$tau2Vec/n.cores
        tau2_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$operator_list$tau2)),1,n.cores)
        tau2 <- apply(tau2_v,1,mean)
        tau2_var <- apply(tau2_v,1,var)/n.cores  
      } else {
        tau_vec <- est.list[[1]]$operator_list$tauVec/n.cores
        tau_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$operator_list$tau)),1,n.cores)
        tau <- apply(tau_v,1,mean)
        tau_var <- apply(tau_v,1,var)/n.cores    
      }
      
      use.kappa = FALSE
      if(bivariate){
        if (!is.null(est.list[[1]]$operator_list$kappa1)){
          use.kappa = TRUE
          kappa1_vec <- est.list[[1]]$operator_list$kappa1Vec/n.cores
          kappa1_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$operator_list$kappa1)),1,n.cores)
          kappa1 <- apply(kappa1_v,1,mean)
          kappa1_var <- apply(kappa1_v,1,var)/n.cores  
          
          kappa2_vec <- est.list[[1]]$operator_list$kappa2Vec/n.cores
          kappa2_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$operator_list$kappa2)),1,n.cores)
          kappa2 <- apply(kappa2_v,1,mean)
          kappa2_var <- apply(kappa2_v,1,var)/n.cores  
        }
      } else {
        if (!is.null(est.list[[1]]$operator_list$kappa)){
          use.kappa = TRUE
          kappa_vec <- est.list[[1]]$operator_list$kappaVec/n.cores
          kappa_v <- matrix(unlist(lapply(1:n.cores,function(x) est.list[[x]]$operator_list$kappa)),1,n.cores)
          kappa <- apply(kappa_v,1,mean)
          kappa_var <- apply(kappa_v,1,var)/n.cores  
        }  
      }
      
      for(i in 2:n.cores){
        if(bivariate){
          tau1_vec <- tau1_vec + est.list[[i]]$operator_list$tau1Vec/n.cores
          tau2_vec <- tau2_vec + est.list[[i]]$operator_list$tau2Vec/n.cores
          if(use.kappa){
            kappa1_vec <- kappa1_vec + est.list[[i]]$operator_list$kappa1Vec/n.cores
            kappa2_vec <- kappa2_vec + est.list[[i]]$operator_list$kappa2Vec/n.cores
          }
        } else {
          tau_vec <- tau_vec + est.list[[i]]$operator_list$tauVec/n.cores
          if(use.kappa)
            kappa_vec <- kappa_vec + est.list[[i]]$operator_list$kappaVec/n.cores  
        }
        
      }  
      if(bivariate){
        est.merge$operator_list$tau1 <- tau1
        est.merge$operator_list$tau1Vec <- tau1_vec
        est.merge$operator_list$tau1_var <- tau1_var
        est.merge$operator_list$tau2 <- tau2
        est.merge$operator_list$tau2Vec <- tau2_vec
        est.merge$operator_list$tau2_var <- tau2_var
        
        est.merge$operator_tau <- matrix(c(tau1,tau2),1,2)
        est.merge$operator_tau_var <- matrix(c(tau1_var,tau2_var),1,2)
        est.merge$operator_tau_vec <- matrix(c(tau1,tau2),1,2)
      } else {
        est.merge$operator_list$tau <- tau
        est.merge$operator_list$tauVec <- tau_vec
        est.merge$operator_list$tau_var <- tau_var
        est.merge$operator_tau <- tau
        est.merge$operator_tau_var <- tau_var
        est.merge$operator_tau_vec <- tau#tau_vec  
      }
      
      if(use.kappa){
        if(bivariate){
          est.merge$operator_list$kappa1 <- kappa1
          est.merge$operator_list$kappa1Vec <- kappa1_vec
          est.merge$operator_list$kappa1_var <- kappa1_var
          est.merge$operator_list$kappa2 <- kappa2
          est.merge$operator_list$kappa2Vec <- kappa2_vec
          est.merge$operator_list$kappa2_var <- kappa2_var
          
          est.merge$operator_kappa <- matrix(c(kappa1,kappa2),1,2)
          est.merge$operator_kappa_var <- matrix(c(kappa1_var,kappa2_var),1,2)
          est.merge$operator_kappa_vec <- matrix(c(kappa1,kappa2),1,2)
        } else {
          est.merge$operator_list$kappa <- kappa
          est.merge$operator_list$kappaVec <- kappa_vec
          est.merge$operator_list$kappa_var <- kappa_var
          est.merge$operator_kappa <- kappa
          est.merge$operator_kappa_var <- kappa_var
          est.merge$operator_kappa_vec <- kappa#kappa_vec  
        }
        
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
    bivariate = FALSE
    if(obj1$measurementError_list$noise == "nsNormal"){
      bivariate = TRUE
    }
    if(bivariate){
      output$measurementError_list$theta_vec <- rbind(matrix(obj1$measurementError_list$theta_vec),
                                                      matrix(obj2$measurementError_list$theta_vec))
      output$meas_error_sigma_vec <- rbind(obj1$meas_error_sigma_vec,obj2$meas_error_sigma_vec)
      output$meas_error_sigma_var <- rbind(obj1$meas_error_sigma_var,obj2$meas_error_sigma_var)
      
    } else {
      output$measurementError_list$sigma_vec <- rbind(matrix(obj1$measurementError_list$sigma_vec),
                                                      matrix(obj2$measurementError_list$sigma_vec))
      output$meas_error_sigma_vec <- rbind(obj1$meas_error_sigma_vec,obj2$meas_error_sigma_vec)
      output$meas_error_sigma_var <- rbind(obj1$meas_error_sigma_var,obj2$meas_error_sigma_var)  
    }
    
    
    if(!is.null(obj2$measurementError_list$nu_vec)){
      output$measurementError_list$nu_vec <- rbind(matrix(obj1$measurementError_list$nu_vec),
                                                   matrix(obj2$measurementError_list$nu_vec))
      #output$meas_error_nu_vec <- output$measurementError_list$nu_vec
      output$meas_error_nu_vec <- rbind(obj1$meas_error_nu_vec,obj2$meas_error_nu_vec)
      output$meas_error_nu_var <- rbind(obj1$meas_error_nu_var,obj2$meas_error_nu_var)
      
    }
    if(!is.null(obj2$measurementError_list$mu_vec)){
      output$measurementError_list$mu_vec <- rbind(matrix(obj1$measurementError_list$mu_vec),
                                                   matrix(obj2$measurementError_list$mu_vec))
      #output$meas_error_mu_vec <- output$measurementError_list$mu_vec
      output$meas_error_mu_vec <- rbind(obj1$meas_error_mu_vec,obj2$meas_error_mu_vec)
      output$meas_error_mu_var <- rbind(obj1$meas_error_mu_var,obj2$meas_error_mu_var)
    }
    
    #merge operator list
    if(bivariate){
      if(!is.null(obj2$operator_list$tau1Vec)){
        output$operator_list$tau1Vec <- rbind(matrix(obj1$operator_list$tau1Vec),
                                              matrix(obj2$operator_list$tau1Vec))
        output$operator_list$tau2Vec <- rbind(matrix(obj1$operator_list$tau2Vec),
                                              matrix(obj2$operator_list$tau2Vec))
        
        output$operator_tau_vec <- rbind(obj1$operator_tau_vec,obj2$operator_tau_vec)
        output$operator_tau_var <- rbind(obj1$operator_tau_var,obj2$operator_tau_var)  
      }
    } else {
      if(!is.null(obj2$operator_list$tauVec)){
        output$operator_list$tauVec <- rbind(matrix(obj1$operator_list$tauVec),
                                             matrix(obj2$operator_list$tauVec))
        output$operator_tau_vec <- rbind(obj1$operator_tau_vec,obj2$operator_tau_vec)
        output$operator_tau_var <- rbind(obj1$operator_tau_var,obj2$operator_tau_var)  
      }  
    }
    
    
    if(bivariate){
      if(!is.null(obj2$operator_list$kappa1Vec)){
        output$operator_list$kappa1Vec <- rbind(matrix(obj1$operator_list$kappa1Vec),
                                                matrix(obj2$operator_list$kappa1Vec))
        output$operator_list$kappa2Vec <- rbind(matrix(obj1$operator_list$kappa2Vec),
                                                matrix(obj2$operator_list$kappa2Vec))
        output$operator_kappa_vec <- rbind(obj1$operator_kappa_vec,obj2$operator_kappa_vec)
        output$operator_kappa_var <- rbind(obj1$operator_kappa_var,obj2$operator_kappa_var)
      }  
    } else {
      if(!is.null(obj2$operator_list$kappaVec)){
        output$operator_list$kappaVec <- rbind(matrix(obj1$operator_list$kappaVec),
                                               matrix(obj2$operator_list$kappaVec))
        output$operator_kappa_vec <- rbind(obj1$operator_kappa_vec,obj2$operator_kappa_vec)
        output$operator_kappa_var <- rbind(obj1$operator_kappa_var,obj2$operator_kappa_var)
      }  
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

check.convergence <- function(output,std.lim,silent=FALSE)
  {
  fixed.converged <- simple.convergence.test(output$fixed_est_vec,output$fixed_est_var,std.lim)
  
  ranef_Sigma.converged <- NULL
  if(!is.null(output$ranef_Sigma) && !is.na(output$ranef_Sigma))
      ranef_Sigma.converged <- simple.convergence.test(output$ranef_Sigma_vec,output$ranef_Sigma_var,std.lim)  
    
  ranef_nu.converged <- NULL
  if(!is.null(output$ranef_nu) && !is.na(output$ranef_nu))
    ranef_nu.converged <- simple.convergence.test(output$ranef_nu_vec, output$ranef_nu_var, std.lim)  
  
  ranef_mu.converged <- NULL
  if(!is.null(output$ranef_mu) && !is.na(output$ranef_mu))
    ranef_mu.converged <- simple.convergence.test(output$ranef_mu_vec, output$ranef_mu_var, std.lim)  
  
  random.converged <- c(ranef_Sigma.converged, ranef_nu.converged,ranef_mu.converged)
  
  #check error
  meas_error_sigma.converged <- simple.convergence.test(output$meas_error_sigma_vec, 
                                                        output$meas_error_sigma_var,
                                                        std.lim)  
  meas_error_nu.converged <- NULL
  if(!is.null(output$meas_error_nu) && !is.na(output$meas_error_nu))
    meas_error_nu.converged <- simple.convergence.test(output$meas_error_nu_vec,
                                                       output$meas_error_nu_var,
                                                       std.lim)  
  
  meas_error_mu.converged <- NULL
  if(!is.null(output$meas_error_mu) && !is.na(output$meas_error_mu) && output$measurementError_list$assymetric)
    meas_error_mu.converged <- simple.convergence.test(output$meas_error_mu_vec,
                                                       output$meas_error_mu_var,
                                                       std.lim)  
  meas_error.converged <- c(meas_error_sigma.converged,meas_error_nu.converged,meas_error_mu.converged)
  
  #check operator
  operator_tau.converged <- NULL
  if(!is.null(output$operator_tau) && !is.na(output$operator_tau))
    operator_tau.converged <- simple.convergence.test(output$operator_tau_vec,
                                                      output$operator_tau_var,
                                                      std.lim)  
  
  operator_kappa.converged <- NULL
  if(!is.null(output$operator_kappa) && !is.na(output$operator_kappa))
    operator_kappa.converged <- simple.convergence.test(output$operator_kappa_vec,
                                                        output$operator_kappa_var,
                                                        std.lim)  
  operator.converged <- c(operator_tau.converged,operator_kappa.converged)
  
  #check process
  process_nu.converged <- NULL
  if(!is.null(output$process_nu) && !is.na(output$process_nu))
    process_nu.converged <- simple.convergence.test(output$process_nu_vec,
                                                    output$process_nu_var,
                                                    std.lim)  
  
  process_mu.converged <- NULL
  if(!is.null(output$process_mu) && !is.na(output$process_mu))
    process_mu.converged <- simple.convergence.test(output$process_mu_vec,
                                                    output$process_mu_var,
                                                    std.lim)  
  
  process.converged <- c(process_mu.converged,process_nu.converged)
  
  converged <- TRUE
  if(length(fixed.converged)>0){
    converged <- converged*(sum(!fixed.converged)==0)
    if(!silent){
      cat("Fixed effects converged: ", fixed.converged,"\n")
    }
  }
  if(length(random.converged)>0){
    converged <- converged*(sum(!random.converged)==0)
    if(!silent){
      cat("Random effects converged: ", random.converged,"\n")
    }
  }
  if(length(meas_error.converged)>0){
    converged <- converged*(sum(!meas_error.converged)==0)
    if(!silent){
      cat("Measurement error converged: ", meas_error.converged,"\n")
    }
  }
  if(length(operator.converged)>0){
    converged <- converged*(sum(!operator.converged)==0)
    if(!silent){
      cat("Operator converged: ", operator.converged,"\n")
    }
  }
  if(length(process.converged)>0){
    converged <- converged*(sum(!process.converged)==0)
    if(!silent){
      cat("Process converged: ", process.converged,"\n")
    }
  }
  return(list(fixed = fixed.converged,
              random = random.converged,
              operator = operator.converged,
              process = process.converged,
              error = meas_error.converged,
              converged = converged))
  
}

simple.convergence.test <- function(m,sigma2,std.lim){
  if(!is.null(dim(m))){
    n.test <- dim(m)[2]
    N <- dim(m)[1]
    output <- rep(FALSE,n.test)
    if(N>3){
      n.points <- min(N,4)
      B <- cbind(rep(1,n.points),1:n.points)
      for(i in 1:n.test){
        std.satisfied <- m[N,i]/sqrt(sigma2[N,i])>std.lim
        Sigma <- diag(sigma2[(N-n.points+1):N,i])
        Q <- solve(t(B)%*%solve(Sigma,B))
        beta <- Q%*%(t(B)%*%solve(Sigma,m[(N-n.points+1):N,i]))
        slope.satisfied <- abs(beta[2])<2*sqrt(Q[2,2]) #no significant trend
        output[i] = std.satisfied&slope.satisfied
      }
    }
    return(output)
  } else {
    if(length(m)>3){
      #check that estimate/std is above the limit
      N = length(m)
      std.satisfied <- m[N]/sqrt(sigma2[N])>std.lim
      
      #check if we have a significant trend in the n.points last points
      n.points <- min(length(m),4)
      B <- cbind(rep(1,n.points),1:n.points)
      Sigma <- diag(sigma2[(N-n.points+1):N])
      Q <- solve(t(B)%*%solve(Sigma,B))
      beta <- Q%*%(t(B)%*%solve(Sigma,m[(N-n.points+1):N]))
      slope.satisfied <- abs(beta[2])<2*sqrt(Q[2,2]) #no significant trend
      return(std.satisfied&slope.satisfied)  
    } else {
      return(FALSE)
    }  
  }
}

make.plot <- function(output,est.list,ii,nIter,list_name,vec_name,point_name,point_var_name,title){
  min.y <- min(output[[list_name]][[vec_name]])
  max.y <- max(output[[list_name]][[vec_name]])
  n.cores <- length(est.list)
  for(i in 1:n.cores){
    min.i <- min(est.list[[i]][[list_name]][[vec_name]])
    max.i <- max(est.list[[i]][[list_name]][[vec_name]])
    min.y <- min(min.y,min.i)
    max.y <- max(max.y,max.i)
  }
  plot(output[[list_name]][[vec_name]],type="l", main = title,ylim=c(min.y,max.y),xaxt="n")  
  for(i in 1:n.cores)
    lines((ii-1)*nIter + (1:nIter), est.list[[i]][[list_name]][[vec_name]],col="gray")  
  lines(output[[list_name]][[vec_name]])
  points((1:ii)*nIter,output[[point_name]],col=2)
  points((1:ii)*nIter,output[[point_name]]+2*sqrt(output[[point_var_name]]),col=3)
  points((1:ii)*nIter,output[[point_name]]-2*sqrt(output[[point_var_name]]),col=3)   
  for(i in 1:ii){
    lines(i*nIter*c(1,1),c(output[[point_name]][i]-2*sqrt(output[[point_var_name]][i]),
                           output[[point_name]][i]+2*sqrt(output[[point_var_name]][i])),col=3)
  }  
}

make.plot.k <- function(output,est.list,ii,nIter,list_name,vec_name,point_name,point_var_name,k,k2,title){
  if(k==0){
    min.y <- min(output[[list_name]][[vec_name]])
    max.y <- max(output[[list_name]][[vec_name]])
    n.cores <- length(est.list)
    for(i in 1:n.cores){
      min.i <- min(est.list[[i]][[list_name]][[vec_name]])
      max.i <- max(est.list[[i]][[list_name]][[vec_name]])
      min.y <- min(min.y,min.i)
      max.y <- max(max.y,max.i)
    }
    plot(output[[list_name]][[vec_name]],type="l", main = title,ylim=c(min.y,max.y),xaxt="n")  
    for(i in 1:n.cores)
      lines((ii-1)*nIter + (1:nIter), est.list[[i]][[list_name]][[vec_name]],col="gray")  
    lines(output[[list_name]][[vec_name]])
  } else {
    min.y <- min(output[[list_name]][[vec_name]][,k])
    max.y <- max(output[[list_name]][[vec_name]][,k])
    n.cores <- length(est.list)
    for(i in 1:n.cores){
      min.i <- min(est.list[[i]][[list_name]][[vec_name]][,k])
      max.i <- max(est.list[[i]][[list_name]][[vec_name]][,k])
      min.y <- min(min.y,min.i)
      max.y <- max(max.y,max.i)
    }
    plot(output[[list_name]][[vec_name]][,k],type="l", main = title,ylim=c(min.y,max.y),xaxt="n")  
    for(i in 1:n.cores)
      lines((ii-1)*nIter + (1:nIter), est.list[[i]][[list_name]][[vec_name]][,k],col="gray")  
    lines(output[[list_name]][[vec_name]][,k]) 
  }
  points((1:ii)*nIter,output[[point_name]][,k2],col=2)
  points((1:ii)*nIter,output[[point_name]][,k2]+2*sqrt(output[[point_var_name]][,k2]),col=3)
  points((1:ii)*nIter,output[[point_name]][,k2]-2*sqrt(output[[point_var_name]][,k2]),col=3)   
  for(i in 1:ii){
    lines(i*nIter*c(1,1),c(output[[point_name]][i,k2]-2*sqrt(output[[point_var_name]][i,k2]),
                           output[[point_name]][i,k2]+2*sqrt(output[[point_var_name]][i,k2])),col=3)
  }  
}

plot.output <- function(output,est.list,ii,nIter,plot.type){
  
  if(plot.type=="Fixed" || plot.type==TRUE || plot.type=="All"){
    n.cores = length(est.list)
    #check which parameters we have 
    n.fixed = dim(output$mixedEffect_list$betaf_vec)[2]
    if(is.null(output$mixedEffect_list$betar_vec)){
      n.random = 0
    } else {
      n.random = dim(output$mixedEffect_list$betar_vec)[2]  
    }
    
    n.effects <- n.fixed + n.random
    
    n.random.sigma = n.random*(n.random+1)/2
    n.random.nu = !(is.null(output$ranef_nu)|is.na(output$ranef_nu))
    n.random.mu = !(is.null(output$ranef_mu)|is.na(output$ranef_mu))
    n.random.dist <- n.random.sigma + n.random.mu + n.random.nu
    
    n.operator.kappa = 0
    if(!(is.null(output$operator_kappa)|is.na(output$operator_kappa[1])))
      n.operator.kappa = length(output$operator_kappa)  
    
    n.operator.tau = 0
    if(!(is.null(output$operator_tau)|is.na(output$operator_tau[1])))
      n.operator.tau = length(output$operator_tau)
    
    n.process.nu = 0
    if(!(is.null(output$process_nu)|is.na(output$process_nu)))
      n.process.nu = length(output$process_nu)
    
    n.process.mu = 0
    if(!(is.null(output$process_mu)|is.na(output$process_mu)))
      n.process.mu = length(output$process_mu)
    
    n.process <- n.operator.tau + n.operator.kappa + n.process.mu + n.process.nu
    
    n.meas.nu = 0
    if(!(is.null(output$meas_error_nu)|is.na(output$meas_error_nu)))
      n.meas.nu = length(output$meas_error_nu)
    
    n.meas.mu = 0
    if((!is.null(output$measurementError_list$assymetric)) && (output$measurementError_list$assymetric==1))
      n.meas.mu = length(output$meas_error_mu)
    
    n.error = length(output$meas_error_sigma) + n.meas.nu + n.meas.mu
    n.tot = n.effects + n.random.dist + n.process + n.error 
    
    bivariate = FALSE
    if(output$measurementError_list$noise == "nsNormal"){
      bivariate = TRUE
    }
    
    #Plot the fixed effects
    
    if(plot.type=="Fixed" || plot.type==TRUE){
      n.p = min(ceiling(sqrt(n.effects)),4)
      if(n.p*(n.p-1)>=n.effects)
        par(mfcol = c(n.p-1, n.p), mai = c(0.05, 0.3, 0.1, 0.05))
      else
        par(mfcol = c(n.p, n.p), mai = c(0.05, 0.3, 0.1, 0.05))
    } else if(plot.type=="All"){
      n.p = min(ceiling(sqrt(n.tot)),4)
      if(n.p*(n.p-1)>=n.tot)
        par(mfcol = c(n.p-1, n.p), mai = c(0.05, 0.3, 0.1, 0.05))
      else
        par(mfcol = c(n.p, n.p), mai = c(0.05, 0.3, 0.1, 0.05))
    }
    
    total.plotted = 0
    for(k in 1:min(n.fixed,16)){
      total.plotted = total.plotted + 1
      make.plot.k(output,est.list,ii,nIter,
                  "mixedEffect_list","betaf_vec","fixed_est_vec","fixed_est_var",
                  k,output$index_fixed[k],"fixed")
    }
    if(n.random>0){
      for(k in 1:n.random){
        if(total.plotted < 16){
          total.plotted = total.plotted + 1 
          make.plot.k(output,est.list,ii,nIter,
                      "mixedEffect_list","betar_vec","fixed_est_vec","fixed_est_var",
                      k,output$index_random[k],"random")
        }
      }  
    }
    
    if(plot.type=="All"){
      if(n.random>0){
        kk = 0
        for(i.r in 1:n.random){
          for(j.r in 1:n.random){
            kk=kk+1
            if(j.r>=i.r && total.plotted < 16){
              total.plotted = total.plotted + 1
              make.plot.k(output,est.list,ii,nIter,
                          "mixedEffect_list","Sigma_vec","ranef_Sigma_vec","ranef_Sigma_var",
                          kk,kk,"Sigma random")
            }
          }
        }  
      }
      if(n.random.mu>0 & total.plotted < 16){
        total.plotted = total.plotted + 1
        make.plot(output,est.list,ii,nIter,
                  "mixedEffect_list","mu_vec","ranef_mu_vec","ranef_mu_var","mu random")
      }
      if(n.random.nu>0 & total.plotted < 16){
        total.plotted = total.plotted + 1
        make.plot(output,est.list,ii,nIter,
                  "mixedEffect_list","nu_vec","ranef_nu_vec","ranef_nu_var","nu random")
      }
      
      if(total.plotted < 16){
        if(bivariate){
          for(k in 1:2){
            if(total.plotted < 16){
              total.plotted = total.plotted + 1
              make.plot.k(output,est.list,ii,nIter,
                          "measurementError_list","theta_vec","meas_error_sigma_vec",
                          "meas_error_sigma_var",k,k,"sigma error")
            }
          }
        } else {
          total.plotted = total.plotted + 1
          make.plot(output,est.list,ii,nIter,
                    "measurementError_list","sigma_vec","meas_error_sigma_vec","meas_error_sigma_var",
                    "sigma error")
        }
      }
      
      if(total.plotted < 16 & n.meas.nu>0){
        total.plotted = total.plotted + 1
        make.plot(output,est.list,ii,nIter,
                  "measurementError_list","nu_vec","meas_error_nu_vec","meas_error_nu_var",
                  "nu error")
      }
      
      if(total.plotted < 16 & n.meas.mu>0){
        total.plotted = total.plotted + 1
        make.plot(output,est.list,ii,nIter,
                  "measurementError_list","mu_vec","meas_error_mu_vec","meas_error_mu_var",
                  "mu error")
      }
      
      if(n.operator.tau>0 & total.plotted < 16){
        if(bivariate){
          total.plotted = total.plotted + 1
          make.plot.k(output,est.list,ii,nIter,
                      "operator_list","tau1Vec","operator_tau_vec","operator_tau_var",0,1,"tau1")
          if(total.plotted < 16){
            total.plotted = total.plotted + 1
            make.plot.k(output,est.list,ii,nIter,
                        "operator_list","tau2Vec","operator_tau_vec","operator_tau_var",0,2,"tau2")  
          }
        } else {
          total.plotted = total.plotted + 1
          make.plot(output,est.list,ii,nIter,
                    "operator_list","tauVec","operator_tau_vec","operator_tau_var","tau")  
        }
      } 
      if(n.operator.kappa>0 & total.plotted < 16){
        if(bivariate){
          total.plotted = total.plotted + 1
          make.plot.k(output,est.list,ii,nIter,
                      "operator_list","kappa1Vec","operator_kappa_vec","operator_kappa_var",0,1,"kappa1")
          if(total.plotted < 16){
            total.plotted = total.plotted + 1
            make.plot.k(output,est.list,ii,nIter,
                        "operator_list","kappa2Vec","operator_kappa_vec","operator_kappa_var",0,2,"kappa2")
          }
        } else {
          total.plotted = total.plotted + 1
          make.plot(output,est.list,ii,nIter,
                    "operator_list","kappaVec","operator_kappa_vec","operator_kappa_var","kappa")
        }
      } 
      if(n.process.nu>0 & total.plotted<16){
        total.plotted = total.plotted + 1
        make.plot(output,est.list,ii,nIter,
                  "processes_list","nu_vec","process_nu_vec","process_nu_var","nu process")
      }
      if(n.process.mu>0 & total.plotted<16){
        total.plotted = total.plotted + 1
        make.plot(output,est.list,ii,nIter,
                  "processes_list","mu_vec","process_mu_vec","process_mu_var","mu process")
      }
    }
  }
}