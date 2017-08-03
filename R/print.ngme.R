
#' @title Print function for \code{"ngme"} objects.
#'
#' @description A function to print results stored in the 
#'    objects returned by the \code{"ngme"} function.
#'
#' @param object A fitted object by calling \code{"ngme"}.
#' @param ... Additional arguments; none used currently.
#'
#' @seealso \code{\link{ngme}}, \code{\link{print}}
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   fit
#'   }

print.ngme <- function(object, ...){

  cat("Call:\n")
  print(object$call)

  cat("\n")
  cat("Fixed effects parameter estimates:\n")
  print(object$fixed_est)

  cat("\n")
  cat("Estimate of random effects Sigma matrix:\n")
  print(object$ranef_Sigma)

  if(object$random_distr %in% c("NIG")){
    cat("\n")
    cat("Estimate of random effects mu:", object$ranef_mu, "\n")
    cat("Estimate of random effects nu:", object$ranef_nu, "\n")
  }

  if(object$use_process == TRUE){
    cat("\n")
    if(object$process_distr %in% c("NIG", "GAL")){
      cat("Estimate of process mu:", object$process_mu, "\n")
      cat("Estimate of process nu:", object$process_nu, "\n")
    }
    cat("Estimate of tau (for operator):", object$operator_tau, "\n")
    if(object$operator_type %in% c("matern")){
      cat("Estimate of kappa (for operator):", object$operator_kappa, "\n")
    }
    }

  cat("\n")
  cat("Estimate of the measurement error sigma:", object$meas_error_sigma, "\n")
  if(object$error_distr %in% c("NIG", "tdist")){
  cat("Estimate of the measurement error nu:", object$meas_error_nu, "\n")
  }
}

