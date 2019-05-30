
#' @title Parameter estimation.
#'
#' @description Maximum likelihood model estimation using parallel runs of 
#'   stochastic gradient estimation. See \code{\link{ngme}} and \code{\link{ngme.spatial}} for explanation of the model specification.
#'
#' @param n.cores Number of cores, and the number of parallel chains, to use. Default is 4. 
#' @param std.lim Parameter for convergence criterium. The estimation is stopped when the rate of change 
#' in each parameter is small enough, and when the estimate divided by the estimated Monte Carlo standard
#' deviation for that parameter is larger than std.lim. Default is 10. 
#' @param max.rep The total number of iterations that is run is given by \code{nIter*max.rep}. Convergence is checked
#' after every nIter iterations, and max.rep thus sets how many batches of nIter iterations that should be run at most. 
#' Default is 10. 
#' @param nIter The number of iterations per batch of runs. Default is 1000. 
#' @param plot.type Set to "All" to get parameter trajectories of all estimated parameters. However, at most 16 parameters
#' are plotted at once. Set to "TRUE" or "Fixed" to get plots of only the fixed effects. 
#' @param ... Other parameter needed by \code{\link{ngme}} 
#' @return A list of outputs.
#' @details The function calls \code{\link{ngme}} or \code{\link{ngme.spatial}} internally. See these functions for further information on the actual
#' model specification. When plots of parameter trajectories are shown, the gray lines show the trajectories of the individual
#' runs and the black curve is the estimate obtained by averaging the individual trajectories. The green lines show approximate
#' 95% confidence bands for the estimate.  
#' 
#' @seealso \code{\link{ngme}}, \code{\link{ngme.spatial}} 
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @examples
#'   \dontrun{
#'   data(srft_data)
#'
#'   # transform pwl to decrese it correlation with bage
#'   # then center all the covariates for better convergence
#'   srft_data$pwl2 <- with(srft_data, pwl - bage/1.5)
#'   srft_data[, c("sex_cent", "bage_cent", "fu_cent", "pwl2_cent")] <-
#'     scale(srft_data[, c("sex", "bage", "fu", "pwl2")], scale = FALSE)
#'
#'   # fit the model with normal assumption for random effects, process and error
#'   # covariance function is integrated random walk
#'   set.seed(123)
#'   fit <- ngme.par(n.cores = 5, std.lim = 100,max.rep = 20,
#'                   fixed = log(egfr) ~ sex_cent + bage_cent + fu_cent + pwl2_cent,
#'                   random = ~ 1|id,
#'                   data = srft_data,
#'                   reffects = "Normal",
#'                   process = c("Normal", "fd2"),
#'                   error = "Normal",
#'                   timeVar = "fu",
#'                   use.process = TRUE,
#'                   silent = FALSE,
#'                   mesh = list(cutoff = 1/365,max.dist = 1/12,extend = 0.01),
#'                   controls = list(pSubsample = 0.025,subsample.type = 1))
#'}

ngme.par <- function(n.cores = 4, 
                     std.lim = 10,
                     max.rep = 10,
                     controls = NULL, 
                     controls.init = NULL,
                     nIter = 1000,
                     timevar = NULL,
                     location.names = NULL,
                     init.fit = NULL,
                     silent = FALSE,
                     plot.type="All",
                     ...)
{
  if(n.cores < 2){
    stop("Must use at least 2 cores. For non-parallel estimation, use ngme for temporal models or ngme.spatial for spatial models.")
  } 
  temporal.model = FALSE
  if(is.null(timevar) && is.null(location.names)){
    stop("You must specify timevar for temporal models or location.names for spatial models")
  } else if(!is.null(timevar) && !is.null(location.names)){
    stop("You can only specify either timevar (if a temporal model is estimiated) or location.names (if a spatial model is estimated).")
  } else if(!is.null(timevar)){
    temporal.model = TRUE
  } 
  call = match.call()
  output <- NULL
  if(is.null(controls)){
    step0  <- 1
    alpha <- 0.6
  } else {
    step0 <- controls$step0
    if(is.null(step0))
      step0 <- 1  
    alpha <- controls$alpha
    if(is.null(alpha))
      alpha = 0.6
  }
  
  for(ii in 1:max.rep){
    cl <- makeCluster(n.cores)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = n.cores, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    parallel::clusterExport(cl, varlist = c(), envir = environment())
    if(is.null(controls)){
      controls <- list()
    }
    controls$alpha = ii*alpha/max.rep
    if(ii>1){
      controls$iter.start <- (ii-1)*nIter
      controls$nBurnin = 5
      
    }
      
    est.list <- foreach(i = 1:n.cores, .options.snow = opts,.packages = c("ngme","Matrix")) %dopar%
    {
      if(ii==1){
        if(temporal.model){
          est <- ngme::ngme(controls=controls,
                            controls.init = controls.init,
                            timevar = timevar,
                            nIter = nIter,
                            init.fit = init.fit,
                            silent = TRUE,...)  
        } else {
          est <- ngme::ngme.spatial(controls=controls,
                            controls.init = controls.init,
                            location.names = location.names,
                            nIter = nIter,
                            init.fit = init.fit,
                            silent = TRUE,...)  
        }
      } else {
        if(temporal.model){
          est <- ngme::ngme(init.fit = est.list.old[[i]],
                            controls=controls,
                            timevar = timevar,
                            nIter = nIter,
                            silent=TRUE,...)  
        } else {
          est <- ngme::ngme.spatial(init.fit = est.list.old[[i]],
                                    controls=controls,
                                    location.names = location.names,
                                    nIter = nIter,
                                    silent=TRUE,...)  
        }
      }
      return(est)  
    }
    close(pb)
    stopCluster(cl)
    est.list.old <- est.list
    est.merge <- merge.ngme.outputs(est.list)
    output <- attach.ngme.output(output,est.merge)
    
    plot.output(output,est.list,ii,nIter,plot.type=plot.type)
    
    converged <- check.convergence(output,std.lim,silent)
    
    if(converged$converged)
          break
  }
  output$call <- call
  return(output)
}

