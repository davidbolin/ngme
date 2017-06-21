#' @title Confidence intervals.
#'
#' @description Calculates confidence intervals for model parameters.
#'
#' @param object A fitted object by calling \code{"nglda_est"}.
#' @param type A character string for the type of parameters;
#'   \code{"fixed"} for fixed effects.
#' @param level A numerical value indicating the confidence level.
#'
#' @details Parameter estimates for random-effects and noise parameters
#'    are not available currently.
#' @seealso \code{\link{nglda_est}}
#'
#' @examples
#'   \dontrun{
#'   fit <- nglda_est(...)
#'   fitted(fit)
#'   }
#'
  intervals <- function(object, type = "fixed", level = 0.95){

    if(!summary(object)$estimate_fisher) stop("Fisher matrix was not estimated")

    if(type == "fixed"){

      fixed_results <- summary(object)$fixed_results
      fixed_est     <- fixed_results[, 1]
      fixed_est_se  <- fixed_results[, 2]

      quant_norm <- abs(qnorm((1 - level)/2))

      llim <- fixed_est - quant_norm * fixed_est_se
      ulim <- fixed_est + quant_norm * fixed_est_se

      result           <- cbind(llim, fixed_est, ulim)
      rownames(result) <- names(fixed_est)
      colnames(result) <- c("lower", "estimate", "upper")

      out <- list(result = result,
                  type = type,
                  level = level)
      class(out) <- "intervals"
      out
    }

  }
