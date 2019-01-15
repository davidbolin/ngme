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

#' @title Group individuals.
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
  debug = 0
  ncov = dim(B[[1]])[2] #number of covariates
  if(debug)
    cat("Desired rank = ",ncov, "\n")
  nper = length(B) #number of subjects
  BB = list()
  grouped = rep(0,nper) #indicator of which subjects who have been assined to groups, will contain the group number for each subject
  
  #compute matrices for each subject
  for(i in 1:nper){
    BB[[i]] = t(B[[i]])%*%B[[i]]
  }
  
  groups <- Bgroup <- list()
  
  #start creating the groups
  groupnumber = 1
  while(sum(abs(grouped)>0)<nper) #while there are subjects left ungrouped
  {
    ng = which(grouped==0) #subjects not grouped (and not having full ranks)
    gi = ng[1] #the first ungrouped subject
    
    Bi = BB[[ng[1]]]
    crank = rankMatrix(Bi)[1] #rank of the subject
    if(debug)
      cat("adding subject ",ng[1], "to group", groupnumber, "rank = ", crank,"\n")
    grouped[ng[1]] = groupnumber #create new group for the subject
    
    if(length(ng)==1){ #if there are no further subjects leave the current ungrouped and stop
      if(crank < ncov)
        grouped[ng[1]] = 0
      break
    } else { #otherwise remove the patient from the ungrouped
      ng = ng[2:length(ng)] 
    }
    
    while(crank<ncov) #while the rank of the group is too small
    {
      Bp = Bi + BB[[ng[1]]] #check the rank if the first ungrouped subject is added to the group
      prank = rankMatrix(Bp)[1]
      if(prank>crank){ #if it is larger than the current rank, add the subject to the group
        gi = c(gi, ng[1])
        grouped[ng[1]] = groupnumber
        Bi = Bp
        crank = prank
        if(debug)
          cat("adding subject ",ng[1], "to group", groupnumber, "rank = ", crank,"\n")
      }
      if(length(ng)==1){ #if there are no further subjects, stop
        break
      } else {
        ng = ng[2:length(ng)] #otherwise remove the subject from the ungrouped and continue
      }
    }
    #at this point, the procedure has stopped because there are no further subjects or because the group is complete.
    
    #if the procedure stopped because there are no further subjects stop the group formation.
    #if the group does not have full rank, ungroup all subjects in the group first
    if(length(ng)==1 && crank < ncov){ 
      if(debug)
        cat("ungrouping group ", groupnumber, "rank = ", crank,"\n")
      grouped[grouped==groupnumber] = 0
      break 
    }
    #if we are here, there are further subjects to group and the group is complete
    if(debug){
      cat("creating group ", groupnumber,"rank =",crank, "\n")
      cat("subjects ", gi, "\n")  
    }
    groups[[groupnumber]] = gi-1 #create the group 
    Bgroup[[groupnumber]] = Bi #save the matrix of the group
    groupnumber = groupnumber +1 #increate group number and start the next group
  }
  free.samples = which(grouped<=0) #create the group G0 out of all ungrouped subjects.
  if(debug)
    cat("creating G0 with samples", free.samples, "\n")
  if(groupnumber == 1){
    Bi = 0
    for(i in 1:nper){
      Bi = Bi + BB[[i]]
    }
    stop("Could not create one group with full rank, must add further subjects to fit the model.")
  }
  return(list(groups = groups,
              free = free.samples-1))
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
