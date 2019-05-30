#'
#' @title Standardise covariates.
#'
#' @description A function to standardise covariates to avoid numerical issues.
#' @param B.list A numerical list that contains covariate matrices.
#'
#' @details This function is supplementary and internally used.
#'
#' @return A list of standardised covariate matrices.
#'
#' @seealso \code{\link{estimateLong}}
#'
#' @examples
#'   \dontrun{
#'   standardize.covariates(...)
#'   }
#'

standardize.covariates <- function(B.list)
{
  if(!is.null(B.list)){
    n.cov <- dim(B.list[[1]])[2]
    if(n.cov>1){
      BB = lapply(B.list,function(x) t(x)%*%x)
      Bs = Reduce("+",BB)/length(B.list)
      R = chol(Bs)
      Ri = solve(R)
      B.out = lapply(B.list,function(x) x%*%Ri)
    } else {
      return(list(B  =B.list, M = matrix(1), Minv = matrix(1)))
    }

    return(list(B  =B.out, M = R, Minv = Ri))

  } else {
    return(NULL)
  }
}

#' @title Scale fixed effects coefficients.
#'
#' @description A function to scale the fixed effects coefficients.
#'
#' @param beta A numerical matrix of fixed effects estimates (in original scale).
#' @param B.list A numerical list of re-scaled covariate matrices.
#' @param inc A logical variable;
#'    \code{"TRUE"} indicates the use of inverse re-scaling
#'    (meaningful to use if \code{"sigma"} is given in transformed scale),
#'    \code{"FALSE"} no use of inverse re-scaling.
#'
#' @details This function is supplementary and internally used.
#'
#' @return A matrix of scaled fixed effects coefficients.
#'
#' @seealso \code{\link{estimateLong}}
#'
#' @examples
#'   \dontrun{
#'   scale.beta(...)
#'   }
#'

scale.beta <- function(beta, B.list, inv = FALSE)
{
  if(!is.matrix(beta))
    beta = as.matrix(beta)

  nc = dim(B.list$M)[1]
  if(inv){
    if(dim(beta)[1]==nc)
      return(B.list$Minv%*%beta)
    else
      return(t(B.list$Minv%*%t(beta)))
  } else {
    if(dim(beta)[1]==nc)
      return(B.list$M%*%beta)
    else
      return(t(B.list$M%*%t(beta)))
  }
}

#' @title Scale covariance matrix.
#'
#' @description A function to scale covariance matrix.
#'
#' @param sigma A numerical matrix of covariance matrix (in original scale).
#' @inheritParams scale.beta
#'
#' @details This function is supplementary and internally used.
#'
#' @return A list of scaled covariance matrix.
#'
#' @seealso \code{\link{estimateLong}}
#'
#' @examples
#'   \dontrun{
#'   scale.sigma(...)
#'   }
#'

scale.sigma <- function(sigma, B.list, inv = FALSE)
{
  if(inv){
    return(B.list$Minv%*%sigma%*%t(B.list$Minv))
  } else {
    return(B.list$M%*%sigma%*%t(B.list$M))
  }
}

#' @title Group individuals for grouped sub-sampler.
#'
#' @description A function to group individuals based on fixed effect matrix B.
#'
#' @param B A numeric list that contains fixed effects covariate matrices
#'
#' @return Returns a list of groups and a vector free.
#'
#' @details This function is supplementary and internally used.
#' Each element in groups contains a vector of elements sufficient to make t(B)*B full rank.
#  The elments in free can be sampled freely given that one of the vectors in groups is sampled
#'
#' @seealso \code{\link{estimateLong}}
#'
#' @examples
#'   \dontrun{
#'   scale.fixed(...)
#'   }
#'
group.fixed <- function(B)
{
  I = 1:length(B)
  groups <- list()
  free <- rep(1,0)
  k = 1
  while(length(I)>0) {
    Gi = create.group(I,B)
    Bi = 0
    for(i in 1:length(Gi)){
      Bi = Bi + t(B[[1+Gi[i]]])%*%B[[1+Gi[i]]]
    }
    if(rankMatrix(Bi) == dim(B[[1]])[2]){
      groups[[k]] = Gi
      I = setdiff(I,1+ groups[[k]])
      k = k+1
    } else {
      free = I-1
      I = NULL
    }
  }
  if(length(groups) == 0){
    stop("Could not create one group with full rank, must add further subjects to fit the model.")
  }
  return(list(groups = groups,
              free = free))
}

#internal function for group creation
create.group <- function(I,B)
{
  G <- I[1]
  Bi = t(B[[I[1]]])%*%B[[I[1]]]
  I <- I[-1]
  while(length(I)>0 && rankMatrix(Bi)[1] < dim(Bi)[2]){
    if(rankMatrix(Bi + t(B[[I[1]]])%*%B[[I[1]]])[1] > rankMatrix(Bi)[1]){
      G <- c(G,I[1])
      Bi = Bi + t(B[[I[1]]])%*%B[[I[1]]]
    }
    I = I[-1]
  }
  return(G-1)
}

#Crop all lists in an ngme object to only those in ind.
crop.lists <- function(object,ind){
  object$measurementError_list$Vs <- object$measurementError_list$Vs[ind]

  object$mixedEffect_list$B_random <- object$mixedEffect_list$B_random[ind]
  object$mixedEffect_list$B_fixed <- object$mixedEffect_list$B_fixed[ind]
  object$mixedEffect_list$U <- t(as.matrix(object$mixedEffect_list$U[ind]))

  object$operator_list$loc <-  object$operator_list$loc[ind]
  object$operator_list$h <-  object$operator_list$h[ind]
  if(object$operator_list$type == "matern"){
    object$operator_list$C <-  object$operator_list$C[ind]
    object$operator_list$Ci <-  object$operator_list$Ci[ind]
    object$operator_list$G <-  object$operator_list$G[ind]
    object$operator_list$Ce <-  object$operator_list$Ce[ind]
  } else if(object$operator_list$type == "fd2"){
    object$operator_list$Q <-  object$operator_list$Q[ind]
  }

  if(object$use_process){
    object$processes_list$X <- object$processes_list$X[ind]
    object$processes_list$V <- object$processes_list$V[ind]
  }

  object$Y <- object$Y[ind]
  object$locs <- object$locs[ind]
  return(object)
}

extract.effects <- function(data = NULL, fixed = NULL, random = NULL, idname = idname,
                            na.action = NULL){
  # response matrix and fixed effects design matrix
  mf_fixed <- model.frame(formula = fixed, data = data,  na.action= na.action)
  y        <- as.matrix(model.extract(mf_fixed, "response"))
  x_fixed_f  <- as.matrix(model.matrix(attr(mf_fixed, "terms"), data = mf_fixed))
  colnames(x_fixed_f)[1] <- gsub("[[:punct:]]", "", colnames(x_fixed_f)[1])  
  id <- as.factor(data[,idname])
  idlist <- unique(id)
  # excluding the intercept and the covariates that are specified in random
  if(!is.null(random)){
    cov_list_fixed  <- attr(terms(fixed), "term.labels")
    cov_list_random <- unlist(strsplit(attr(terms(random), "term.labels"), " | ", fixed = TRUE))
    cov_list_random <- c(strsplit(cov_list_random[1], " + ", fixed=TRUE)[[1]], cov_list_random[2])
    cov_list_random <- cov_list_random[-length(cov_list_random)]  
    to_del_x_fixed <- c("Intercept", cov_list_fixed[(cov_list_fixed %in% cov_list_random)])
    x_fixed <- x_fixed_f[, !(colnames(x_fixed_f) %in% to_del_x_fixed), drop = FALSE]
    #random effects design matrix
    random_names             <- unlist(strsplit(as.character(random)[-1], " | ", fixed = TRUE))
    random_names_id_excluded <- random_names[!(random_names %in% idname)]
    random_formula           <- as.formula(paste("~", paste(random_names_id_excluded, collapse = "+")))
    
    mf_random <- model.frame(formula = random_formula, data = data)
    x_random  <- as.matrix(model.matrix(attr(mf_random, "terms"), data = mf_random))
    colnames(x_random)[1] <- gsub("[[:punct:]]", "", colnames(x_random)[1])
    
    data_random <- data.frame(cbind(id, x_random))
    B_random    <- split(data_random[, -1], data_random[,1])
    B_random    <- lapply(B_random, function(x) as.matrix(x))
    
    Y <- tapply(y, id, function(x) x)
    data_fixed <- data.frame(cbind(id, x_fixed))
    B_fixed    <- split(data_fixed[, -1], data_fixed[,1])
    B_fixed    <- lapply(B_fixed, function(x) as.matrix(x))  
  } else {
    x_fixed <- x_fixed_f
    x_random <- NA
    to_del_x_fixed = NA
    if(is.null(idname)){ # no replicates
      Y <- list(y)
      B_fixed = list(x_fixed)
    } else {
      Y <- tapply(y, id, function(x) x)
      data_fixed <- data.frame(cbind(id, x_fixed))
      B_fixed    <- split(data_fixed[, -1], data_fixed[,1])
      B_fixed    <- lapply(B_fixed, function(x) as.matrix(x))
    }
  }
  return(list(Y=Y,
              B_fixed = B_fixed,
              x_fixed = x_fixed,
              x_random = x_random,
              x_fixed_f = x_fixed_f,
              to_del_x_fixed = to_del_x_fixed))
}

