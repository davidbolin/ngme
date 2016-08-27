
#' estimation function    
#' @param input_list -     
#' @param   Y           - list with the observations
#' @param   Niter       - number of iteration of the stochastic gradient
#' @param   meas_noise  - the aviable noise classes: Normal or NIG
#' @param   noise       - the distribution of the mixed effect
#' @param   B_random    - list for the random effect covariates (needs to be matrix, can be NULL)
#' @param   B_fixed     - list for the fixed  effect covariates (needs to be matrix, can be NULL) 
#' @param   beta_random - initial parameters of the random effect (mean parameter) (if not specified set to zero)
#' @param   beta_fixed  - initial parameters of the fixed  effect (if not specified set to zero)
#' @param   Sigma       - initial parameters of the covariance of random effect (if not specified set to I )
#' @param   nu          - shape parameter for noise (NIG only)
#' @param   mu          - shift parameter for noise (NIG only)
#' @param   U           - (list) inital guess of the random effect
#' @param   V           - (list) inital guess of the variance effect
#' @param   meas_list   - list for measurement error:
#' @param     sigma       - measurement noise variance 
#' @param     nu          - shape parameter for noise (NIG only)
#' @param     Vs          - (list) inital guess for the noise of the measurement
#' output::
#' @param output_list - list
#' 
#' @param   Y_list       - the observations
#' @param   mu           -  tracjorty of the shift parameter (of the mixed effect)
#' @param   nuVec        -  tracjorty of the shape parameter (of the mixed effect)
#' @param   nu_measerror -  tracjorty of the shape parameter (of the measurement error)
#' @param   betaVec      -  tracjorty of the random effect mean parameter
#' @param   sigma_eps    -  tracjorty of the standard devation parameter (of measurement error)
#' @param    mixedeffect - list the output from mixed effects
#' @param       noise       - string of distribution type
#' @param       B_random    - list for the random effect covariates
#' @param       B_fixed     - list for the fixed  effect covariates
#' @param       beta_random - the last iteration of beta_random (mean of random effect) 
#' @param       beta_fixed  - the last iteration of beta_fixed 
#' @param       Sigma       - the last iteration of Sigma ("covaraince of random effect")
#' @param       nu          - the last iteration of shape parameter
#' @param       mu          - the last iteration of shift parameter
#' @param       U           - the last sample of random effect
#' @param       V           - the last sample of the variance coeffient of random effect
#' 
#' @param    measerror - list the output from measurement error
#' @param       noise       - string of distribution type of the measurement error
#' @param       sigma       - the last iteration of standard devation parameter
#' @param       nu          - the last iteration of shape parameter
#' @param       Vs           - the last sample of the random variance of the measurement error
estiamteMixedGH <- function(input_list)
{
  return(estimateME(input_list))
} 

#' measurement error input
#' TODO: generlize the class so that one put in only theta
#'       and then theta arguments should be passed on depending on noise class
#'       
#' @param Y     - list with the observations
#' @param noise - the aviable noise classes
#' @param sigma - measurement variance
#' @param nu    - shape parameter for noise (NIG only)
MeasurementErrorInit <- function(Y,
                                 noise = c("Normal","NIG"),
                                 sigma = 1.,
                                 nu    = 1.){
  
  noise <- match.arg(noise)
  output <- list(noise = noise)
  output$sigma = sigma
  
  if(noise == "NIG"){
    output$nu = nu
    output$Vin <- list()
    for(i in 1:length(Y))
    {
      output$Vin[[i]] <- rep(1, length(Y[[i]]))
    }
  }
  return(output)
}


#' Mixed effect input
#' TODO: generlize the class so that one put in only theta
#'       and then theta arguments should be passed on depending on noise class
#'       
#' @param noise - the aviable noise classes
#' @param B_fixed     - fixed effect
#' @param B_random    - random effect
MixedInit <- function(noise = c("Normal","NIG"),
                      B_fixed  = NULL,
                      B_random = NULL){
  
  noise <- match.arg(noise)
  output <- list(noise = noise)
  if(is.null(B_random) == 0)
    output$B_random = B_random 
  
  if(is.null(B_fixed) == 0)
    output$B_fixed = B_fixed
  
  return(output)
}