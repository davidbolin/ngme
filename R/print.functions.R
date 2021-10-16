
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
#' @export
#' @method print ngme
#' 
print.ngme <- function(object, ...){

  cat("Call:\n")
  print(object$call)

  cat("\n")
  cat("Fixed effects parameter estimates:\n")
  print(object$fixed_est)

  cat("\n")
  cat("Estimate of random effects Sigma matrix:\n")
  print(object$ranef_Sigma)

  if(object$random_distr %in% c("NIG", "tdist")){
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

#' @title Print function for the function \code{intervals}.
#'
#' @description Prints the output produced by \code{intervals} function concisely.
#'
#' @param object A fitted object by calling \code{"intervals"}.
#' @param ... Additional arguments; none used currently.

#' @seealso \code{\link{intervals}}
#'
#' @examples
#'   \dontrun{
#'   fit <- nglda_est(...)
#'   intervals(fit)
#'   }
#' @export
#' @method print intervals

print.intervals <- function(object, ...){
  
  if(object$type == "fixed"){
    cat(paste0(object$level*100, "% confidence intervals for the fixed effects"))
    
    cat("\n"); cat("\n")
    
    print(object$result)
  }
  
}


#' @title Print function for \code{"predict.ngme"} objects.
#'
#' @description Prints the output contained in the \code{"predict.ngme"} objects concisely.
#'
#' @param object A fitted object by calling \code{"predict.ngme"}.
#' @param ... Additional arguments; none used currently.
#'
#' @details MAE, RMSE, CRPS, Interval coverage, Interval width
#'
#' @seealso \code{\link{predict.ngme}}, \code{\link{print}}
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   pred <- predict.ngme(fit, ...)
#'   pred
#'   }
#' @export
#' @method print predict.ngme
#' 
print.predict.ngme <- function(object, ...){
  
  cat("Call:\n")
  print(object$call)
  
  cat("\n")
  
  if(object$type == "Filter"){
    cat("Filtering for subjects with the following ID(s) (", length(object$id), " of them in total)", ":\n", sep = "")
    if(length(object$id) > 10){
      cat(head(object$id, 10), "...\n")
    }else{
      cat(object$id, "\n")
    }
  }else{
    cat("Smoothing for subjects with the following ID(s) (", length(object$id), " of them in total)", ":\n", sep = "")
    if(length(object$id) > 10){
      cat(head(object$id, 10), "...\n")
    }else{
      cat(object$id, "\n")
    }
  }
  
  cat("\n")
  
  cat("Mean predictor MAE (SD) = ", object$mean.mae.mean.predictor," (", object$std.mae.mean.predictor, ")\n", sep = "")
  cat("Median predictor MAE (SD) = ", object$mean.mae.median.predictor," (", object$std.mae.median.predictor, ")\n", sep = "")
  cat("Mean predictor RMSE (SD) = ",  object$mean.rmse.mean.predictor," (", object$std.rmse.mean.predictor,")\n", sep = "")
  cat("Median predictor RMSE (SD) = ",  object$mean.rmse.median.predictor," (", object$std.rmse.median.predictor,")\n", sep = "")
  cat("CRPS (SD) = ", object$mean.crps," (", object$std.crps, ")\n", sep = "")
  cat("Interval coverage (SD) = ", object$coverage.mean, " (", object$coverage.std,")\n", sep = "")
  cat("Interval width (SD) = ", object$mean.int.width, " (", object$std.int.width,")\n", sep = "")
  
}


#' @title Print function for \code{"summary.ngme"} function.
#'
#' @description A function to print \code{"summary.ngme"} objects concisely.
#'
#' @param object A fitted object returned by the \code{"summary.ngme"} function.
#' @param ... Additional arguments; none used currently.
#'
#' @seealso \code{\link{ngme}}, \code{\link{summary.ngme}},
#'   \code{\link{print}}, \code{\link{summary}}
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   summary(fit)
#'   }
#'
#' @export
#' @method print summary.ngme
#' 
print.summary.ngme <- function(object, ...){
  
  cat("Call:\n")
  print(object$call)
  
  if(object$estimate_fisher){#fisher matrix is estimated
    cat("\n"); cat("\n")
    cat("Fixed effects:\n")
    cat("\n")
    printCoefmat(object$fixed_results, P.valueS = T, has.Pvalue = T)
    
    cat("\n"); cat("\n")
    cat("Random effects:\n")
    cat("\n")
    cat("Sigma matrix:\n")
    print(object$random_results$Sigma_est)
    cat("\n")
    cat("Standard errors for Sigma matrix:\n")
    print(object$random_results$Sigma_est_se)
    
    if(object$random_distr %in% c("NIG", "tdist")){
      mu_results <-
        matrix(
          c(object$random_results$mu_est,
            object$random_results$mu_est_se),
          ncol = 2, byrow = F,
          dimnames = list(paste0("mu_", 1:length(object$random_results$mu_est)), c("Estimate", "SE"))
        )
      
      nu_results <-
        matrix(c(object$random_results$nu_est,
                 object$random_results$nu_est_se), nrow = 1,
               dimnames = list(c("nu"), c("Estimate", "SE")))
      
      cat("\n")
      print(mu_results)
      
      cat("\n")
      print(nu_results)
    }
    
    if(object$use_process){
      cat("\n"); cat("\n")
      cat("Operator:\n")
      cat("\n")
      print(object$operator_results)
      
      if(object$process_distr %in% c("NIG", "GAL")){
        cat("\n"); cat("\n")
        cat("Process:\n")
        cat("\n")
        print(object$process_results)
      }
    }
    
    cat("\n"); cat("\n")
    cat("Measurement error:\n")
    cat("\n")
    print(object$meas_error_results)
    
  }else{#fisher is not estimated
    
    cat("\n"); cat("\n")
    cat("Fixed effects:\n")
    cat("\n")
    print(object$fixed_results)
    
    cat("\n"); cat("\n")
    cat("Random effects:\n")
    cat("\n")
    cat("Sigma matrix:\n")
    print(object$random_results$Sigma_est)
    
    if(object$random_distr %in% c("NIG", "tdist")){
      mu_results <-
        matrix(
          object$random_results$mu_est,
          ncol = 1,
          dimnames = list(paste0("mu_", 1:length(object$random_results$mu_est)), c("Estimate"))
        )
      
      nu_results <-
        matrix(object$random_results$nu_est,
               nrow = 1,
               dimnames = list(c("nu"), c("Estimate")))
      
      cat("\n")
      print(mu_results)
      
      cat("\n")
      print(nu_results)
    }
    
    if(object$use_process){
      
      cat("\n"); cat("\n")
      cat("Operator:\n")
      cat("\n")
      print(object$operator_results)
      
      if(object$process_distr %in% c("NIG", "GAL")){
        cat("\n"); cat("\n")
        cat("Process:\n")
        cat("\n")
        print(object$process_results)
      }
    }
    
    cat("\n"); cat("\n")
    cat("Measurement error:\n")
    cat("\n")
    print(object$meas_error_results)
    
  }
  
}
