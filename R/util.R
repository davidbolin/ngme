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
  ncov = dim(B[[1]])[2]
  nper = length(B)
  BB = list()
  grouped = rep(0,nper)
  for(i in 1:nper){
    BB[[i]] = t(B[[i]])%*%B[[i]]
    if(rankMatrix(BB[[i]])[1] == ncov)
      grouped[i] = -1
  }
  groups <- Bgroup <- list()
  groupnumber = 1
  while(sum(abs(grouped)>0)<nper)
  {
    ng = which(grouped==0)
    gi = ng[1]
    Bi = BB[[ng[1]]]
    crank = rankMatrix(Bi)[1]

    grouped[ng[1]] = groupnumber

    if(length(ng)==1){
      if(crank < ncov)
        grouped[ng[1]] = 0
      break
    } else {
      ng = ng[2:length(ng)]
    }

    while(crank<ncov)
    {
      Bp = Bi + BB[[ng[1]]]
      prank = rankMatrix(Bp)[1]
      if(prank>crank){
        gi = c(gi, ng[1])
        grouped[ng[1]] = groupnumber
        Bi = Bp
        crank = prank
      }
      if(length(ng)==1){
        break
      } else {
        ng = ng[2:length(ng)]
      }
    }
    if(length(ng)==1){
      if(crank < ncov)
        grouped[grouped==groupnumber] = 0
      break
    }
    groups[[groupnumber]] = gi-1
    Bgroup[[groupnumber]] = Bi
    groupnumber = groupnumber +1
  }
  free.samples = which(grouped<=0)
  return(list(groups = groups,
              free = free.samples-1))
}
