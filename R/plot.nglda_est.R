#' @title Trace plots.
#'
#' @description Plots the parameter estimates at successive iterations
#'    of the stochastic algorithm.
#'
#' @param object A fitted object by calling \code{"nglda_est"}.
#' @param param A character string for which parameter set the traces
#'   are to be plotted;
#'   \code{"fixed"} for fixed effects,
#'   \code{"random"} for time-invariant random effects,
#'   \code{"process"} for stochastic process,
#'   \code{"error"} for error term.
#' @param ... Additional arguments; none used currently.
#'
#' @seealso \code{\link{nglda_est}}
#'
#' @examples
#'   \dontrun{
#'   fit <- nglda_est(...)
#'   fitted(fit)
#'   }
#'

plot.nglda_est <- function(object, param = "fixed", ...){

  # if process was not specified in the model
  if(param == "process" && object$use_process == FALSE){
    stop("No process was included in the model fit")
  }


  if(param == "fixed"){

    plot.ts(object$fixed_est_vec, main = "Fixed effects", xlab = "Iteration")#, ...)

  }else if(param == "random"){

    # preparing the ranef_Sigma_vec: excluding duplicated columns and adding column names
    ncol_Sigma        <- ncol(summary(object)$random_results$Sigma_est)
    cols_omit         <- unlist(lapply(1:(ncol_Sigma-1), function(i) (ncol_Sigma * i + 1):(ncol_Sigma * i + i)))
    ranef_results_vec <- object$ranef_Sigma_vec
    ranef_results_vec <- ranef_results_vec[, -cols_omit]

    colnames_Sigma_est          <- colnames(summary(object)$random_results$Sigma_est)
    #colnames_Sigma_est          <- substr(colnames_Sigma_est, 1, 3)
    colnames_Sigma_est_ext      <- paste(rep(colnames_Sigma_est,each = ncol_Sigma), colnames_Sigma_est, sep = "_")
    colnames_Sigma_est_ext      <- colnames_Sigma_est_ext[-cols_omit]
    colnames(ranef_results_vec) <- colnames_Sigma_est_ext

    if(summary(object)$random_distr %in% ("Normal")){
      plot.ts(ranef_results_vec, main = "Random effects", xlab = "Iteration")#, ...)
    }else if(summary(object)$random_distr %in% ("NIG")){
      ranef_results_vec <- cbind(ranef_results_vec, object$ranef_mu_vec, object$ranef_nu_vec)
      colnames(ranef_results_vec)[(ncol(ranef_results_vec)-ncol_Sigma ) : (ncol(ranef_results_vec))] <-
        c(paste("mu", colnames_Sigma_est, sep = "_"), "nu")
      plot.ts(ranef_results_vec, main = "Random effects", xlab = "Iteration")#, ...)

    }

  }else if(param == "process"){

    if(object$process_distr %in% c("Normal", "CH")){
      plot(object$operator_tau_est_vec, main = "Process",
           xlab = "Iteration", ylab = "tau", type = "l")#, ...)
    }else if(object$process_distr %in% c("NIG", "GAL")){
      process_results_vec <- cbind(object$operator_tau_est_vec,
                                   object$process_mu_vec,
                                   object$process_nu_vec)
      colnames(process_results_vec) <- c("tau", "mu", "nu")
      plot.ts(process_results_vec, main = "Process", xlab = "Iteration")#, ...)
    }


  }else if(param == "error"){

    if(object$error_distr == "Normal"){
      plot(object$meas_error_sigma_vec, main = "Measurement error",
           xlab = "Iteration", ylab = "sigma", type = "l")#, ...)
    }else if(object$error_distr %in% c("NIG", "tdist")){
      meas_error_result_vec <- cbind(object$meas_error_sigma_vec,
                                   object$meas_error_nu_vec)
      colnames(meas_error_result_vec) <- c("sigma", "nu")
      plot.ts(meas_error_result_vec, main = "Measurement error", xlab = "Iteration")#, ...)
    }
  }
}


