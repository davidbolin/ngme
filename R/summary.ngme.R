
#' @title Summary function for \code{"ngme"} objects.
#'
#' @description A function to summarise the output contained in 
#'    \code{"ngme"} objects.
#'
#' @param object A fitted object by calling \code{"ngme"}.
#' @param ... Additional arguments; none used currently.
#'
#' @seealso \code{\link{ngme}}, \code{\link{print.summary.ngme}},
#'   \code{\link{print}}, \code{\link{summary}}
#' @examples
#'   \dontrun{
#'   fit <- ngme(...)
#'   summary(fit)
#'   }

summary.ngme <- function(object, ...){

  if(object$estimate_fisher){ # if fisher is estimated
    se_est       <- sqrt(diag(solve(object$fisher_est)))

    fixed_est    <- object$fixed_est
    fixed_se_est <- se_est[names(se_est) %in% names(fixed_est)]
    fixed_se_est <- fixed_se_est[names(fixed_est)]
    fixed_z      <- fixed_est/fixed_se_est
    fixed_p      <- pnorm(abs(fixed_z), lower.tail = FALSE) * 2

    fixed_results <- cbind(fixed_est, fixed_se_est, fixed_z, fixed_p)
    colnames(fixed_results) <- c("Estimate", "SE", "Z", "p-value")

    if(object$random_distr == "Normal"){
      random_results <- list(Sigma_est    = object$ranef_Sigma,
                             Sigma_est_se = se_est[grep("Sigma_random", names(se_est))])
    }else if(object$random_distr %in% c("NIG", "tdist")){
      random_results <- list(mu_est       = as.numeric(object$ranef_mu),
                             mu_est_se    = se_est[grep("mu_random", names(se_est))],
                             Sigma_est    = object$ranef_Sigma,
                             Sigma_est_se = se_est[grep("Sigma_random", names(se_est))],
                             nu_est       = object$ranef_nu,
                             nu_est_se    = se_est[names(se_est) == "nu"]
                             )
    }

    if(object$use_process){#if use.process = TRUE
      operator_tau_est <- object$operator_tau
      operator_tau_est_se  <- se_est[names(se_est) == "tau_operator"]
      operator_results <- cbind(operator_tau_est, operator_tau_est_se)
      colnames(operator_results) <- c("Estimate", "SE")
      rownames(operator_results) <- "operator_tau"

      if(object$operator_type %in% c("matern")){
        operator_kappa_est <- object$operator_kappa
        operator_kappa_est_se <- se_est[names(se_est) == "kappa_operator"]
        operator_results <- rbind(operator_results, c(operator_kappa_est, operator_kappa_est_se))
        rownames(operator_results)[2] <- "operator_kappa"
      }

      if(object$process_distr %in% c("NIG", "GAL")){
        process_mu_est    <- object$process_mu
        process_mu_est_se <- se_est[names(se_est) == "mu_process"]
        process_nu_est    <-  object$process_nu
        process_nu_est_se <- se_est[names(se_est) == "nu_process"]

        process_results <- matrix(c(process_mu_est,
                                    process_mu_est_se,
                                    process_nu_est,
                                    process_nu_est_se), ncol = 2, byrow = T)
        rownames(process_results) <- c("mu", "nu")
        colnames(process_results) <- c("Estimate", "SE")
      }
    }

    if(object$error_distr == "Normal"){
      meas_error_sigma_est    <- object$meas_error_sigma
      meas_error_sigma_est_se <- se_est[names(se_est) == "sigma_error"]

      meas_error_results           <- cbind(meas_error_sigma_est, meas_error_sigma_est_se)
      colnames(meas_error_results) <- c("Estimate", "SE")
      rownames(meas_error_results) <- "sigma"
    }else if(object$error_distr %in% c("NIG", "tdist")){
      meas_error_sigma_est    <- object$meas_error_sigma
      meas_error_sigma_est_se <- se_est[names(se_est) == "sigma_error"]
      meas_error_nu_est       <- object$meas_error_nu
      meas_error_nu_est_se    <- se_est[names(se_est) == "nu_error"]

      meas_error_results <- matrix(c(meas_error_sigma_est,
                                     meas_error_sigma_est_se,
                                     meas_error_nu_est,
                                     meas_error_nu_est_se), ncol = 2, byrow = T)
      rownames(meas_error_results) <- c("sigma", "nu")
      colnames(meas_error_results) <- c("Estimate", "SE")

    }

  }else{ # if fisher is not estimated
    fixed_est    <- object$fixed_est

    fixed_results <- matrix(fixed_est, ncol = 1,
                            dimnames = list(names(fixed_est),"Estimate"))

    if(object$random_distr == "Normal"){
      random_results <- list(Sigma_est    = object$ranef_Sigma)
    }else if(object$random_distr %in% c("NIG", "tdist")){
      random_results <- list(mu_est    = as.numeric(object$ranef_mu),
                             Sigma_est = object$ranef_Sigma,
                             nu_est    = object$ranef_nu)
    }

    if(object$use_process){ # if use.process = TRUE
      operator_results <- matrix(object$operator_tau, ncol = 1)
      colnames(operator_results) <- "Estimate"
      rownames(operator_results) <- "operator_tau"

      if(object$operator_type %in% c("matern")){
        operator_results <- rbind(operator_results, object$operator_kappa)
        rownames(operator_results)[2] <- "operator_kappa"
      }

      if(object$process_distr %in% c("NIG", "GAL")){
        process_mu_est    <- object$process_mu
        process_nu_est    <- object$process_nu

        process_results <- matrix(c(process_mu_est,
                                    process_nu_est), ncol = 1)
        rownames(process_results) <- c("mu", "nu")
        colnames(process_results) <- c("Estimate")
      }
    }

    if(object$error_distr == "Normal"){

      meas_error_results <- matrix(c(object$meas_error_sigma))
      rownames(meas_error_results) <- c("sigma")
      colnames(meas_error_results) <- c("Estimate")

    }else if(object$error_distr %in% c("NIG", "tdist")){
      meas_error_sigma_est    <- object$meas_error_sigma
      meas_error_nu_est       <- object$meas_error_nu

      meas_error_results <- matrix(c(meas_error_sigma_est,
                                     meas_error_nu_est), ncol = 1)
      rownames(meas_error_results) <- c("sigma", "nu")
      colnames(meas_error_results) <- c("Estimate")
    }
  }

  # preparing the output
  if(object$use_process){#if use.process = TRUE
    if(object$process_distr %in% c("NIG", "GAL")){
      out <- list(call = object$call,
                  fixed_results = fixed_results,
                  random_results = random_results,
                  operator_results = operator_results,
                  process_results = process_results,
                  meas_error_results = meas_error_results,
                  use_process = object$use_process,
                  random_distr = object$random_distr,
                  process_distr = object$process_distr,
                  operator_type = object$operator_type,
                  error_distr = object$error_distr,
                  estimate_fisher = object$estimate_fisher)

    }else if(object$process_distr %in% c("Normal", "CH")){
      out <- list(call = object$call,
                  fixed_results = fixed_results,
                  random_results = random_results,
                  operator_results = operator_results,
                  meas_error_results = meas_error_results,
                  use_process = object$use_process,
                  random_distr = object$random_distr,
                  process_distr = object$process_distr,
                  operator_type = object$operator_type,
                  error_distr = object$error_distr,
                  estimate_fisher = object$estimate_fisher)
    }

  }else{#if use.process = FALSE
    out <- list(call = object$call,
                fixed_results = fixed_results,
                random_results = random_results,
                meas_error_results = meas_error_results,
                use_process = object$use_process,
                random_distr = object$random_distr,
                process_distr = object$process_distr,
                operator_type = object$operator_type,
                error_distr = object$error_distr,
                estimate_fisher = object$estimate_fisher)
  }

    class(out) <- "summary.ngme"
    out

}
