
nglda <- function(fixed, random, data = NULL, 
                  reffects = "NIG",  error = "NIG", 
                  initials = list(), alpha = 0.1, step0 = 0.33, 
                  Niter = 1000, nSim = 2, silent = TRUE){

  #fixed: two sided formula for fixed effects
  #random: one sided formula for random effects (id will be included in this formula)
  #timeVar: time variable
  #id: id of the patients, just extract from the data
  #data: the data set
  #reffects: "nig" = normal-inverse Gaussian, "gaussian" = for Gaussian 
  #process: "nig" = normal-inverse Gaussian, "gaussian" = for Gaussian
  #error: "nig" = normal-inverse Gaussian, "gaussian" = for Gaussian
  #initials: a list for initial values
  #tolerance: tolerance value to stop the stochastic gradient algorithm - i am not sure what is the stopping rule at the moment
  #Niter: number of iterations for the Stochastic gradient algorithm
  #nSim: number of Gibbs sample in each iteration
  #silent: if FALSE, details of the stochastic gradient are printed when the algorithm is running
  
  #nGibbs: Niter in JW code
  #reffects: noise in JW code, alternative is Normal
  #error: meas_noise in JW code, alternative is Normal
 


  # response matrix and fixed effects design matrix 
  mf_fixed <- model.frame(formula = fixed, data = data)
  y        <- as.matrix(model.extract(mf_fixed, "response"))
  x_fixed  <- as.matrix(model.matrix(attr(mf_fixed, "terms"), data = mf_fixed))
  colnames(x_fixed)[1] <- gsub("[[:punct:]]", "", colnames(x_fixed)[1])
  
  #random effects design matrix and id variable
  idname <- rev(unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE)))[1]
  id <- data[, idname]  
  random_names             <- unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE))
  random_names_id_excluded <- random_names[!(random_names %in% idname)]  
  random_formula           <- as.formula(paste("~", paste(random_names_id_excluded, collapse = "+")))
  
  mf_random <- model.frame(formula = random_formula, data = data)
  x_random  <- as.matrix(model.matrix(attr(mf_random, "terms"), data = mf_random))
  colnames(x_random)[1] <- gsub("[[:punct:]]", "", colnames(x_random)[1])
  
  idlist <- unique(id)     # unique id's - patids in JW code
  #nsubj  <- length(idlist) # number of subjects
  #ntotal <- nrow(y)        # total number of observations

  # converting to lists: fixed effects design matrix, random effects design matrix, response matrix, time variable
  data_fixed  <- data.frame(cbind(id, x_fixed))
  DM_fixed    <- split(data_fixed[, -1], data_fixed[,1]) #B_list in JW code
  DM_fixed    <- lapply(DM_fixed, function(x) as.matrix(x))
  
  data_random <- data.frame(cbind(id, x_random))
  DM_random   <- split(data_random[, -1], data_random[,1])
  DM_random   <- lapply(DM_random, function(x) as.matrix(x))
  
  YM    <- tapply(y, id, function(x) x) #Y_list in JW code
  nobs  <- tapply(id, id, function(x) length(x)) # number of observations per patient
  #Time  <- tapply(timeVar, id, function(x) x) # probably to be used when stochastic process is added
  
  #Vin <- tapply(rep(1, nrow(y)), id, function(x) x) # to be used in meas_list below for the initials
 
  
  ### decomposing initials 
  Sigma <- matrix(Sigma, ncol = sqrt(length(initials$Sigma)))
  beta_random <- initials$beta_random
  sigma_eps <- initials$sigma_eps
  meas_list <-   list(Vs=Vin, sigma = initials$meas_list[[1]], nu = initials$meas_list[[2]])
  nu = initials$nu
  mu = initials$mu
  meas_list$sigma_eps <- sigma_eps
  meas_list$noise <- error
  mixedEffect_list <- list(B_fixed     = DM_fixed, 
                           B_random    = DM_random, 
                           Sigma       = Sigma,
                           beta_random = beta_random,
                           noise       = reffects,
                           nu          = nu,
                           mu          = mu
                           )
  input <- list(Y = YM, 
              Niter       = Niter,
              nSim        = nSim,
              meas_noise  = error,
              alpha       = alpha,
              step0       = step0,
              mixedEffect_list        = mixedEffect_list,
              measurementError_list   = meas_list)
  res_m <- estimateME(input)
  return(res_m)


}
