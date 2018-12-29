#' @title Residual values.
#'
#' @description Calculates residuals.
#'
#' @param object A fitted object returned by the \code{"summary.ngme"} function.
#' @param type A character string for the type of fitted values;
#'   \code{"marginal"} for fixed effects.
#' @param ... Additional arguments; none used currently.
#'
#' @details Subject-specific residaul values are not available currently.
#' @seealso \code{\link{ngme}}, \code{\link{fitted.ngme}}
#'
#'
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   fitted(fit)
#'   }
#'
residuals.ngme <- function(object, type = "marginal", ...){

  if(type == "marginal"){
    resid <- as.numeric(unlist(object$Y) - fitted(object, type = "marginal", ...))
  } else {
    resid <- list()
    n.sub <- length(object$Y)
    for(i in 1:n.sub){
      i.fix <- object$mixedEffect_list$B_fixed[[i]]%*%as.vector(object$mixedEffect_list$beta_fixed)
      i.ran <- object$mixedEffect_list$B_random[[i]]%*%object$mixedEffect_list$U[,i]
      i.pro <- object$processes_list$X[[i]]
      A.i <- build.A.matrix(object$operator_list,object$locs,i)
      resid[[i]] <- as.double(object$Y[[i]] - (i.fix + i.ran + A.i%*%i.pro))
    }
    resid <- as.numeric(unlist(resid))
  }
  return(resid)
}
