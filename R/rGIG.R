#' @title Simulate generalised inverse-Gaussian (GIG) random variables.
#'
#' @description A function to simulate random numbers 
#'     from GIG distribution -- GIG(p, a, b).
#' @param p A numeric vector for p parameters.
#' @param a A numeric vector for a parameters.
#' @param b A numeric vector for b parameters.
#' @param seed A numeric value for setting the seed of random number generation.
#' @details This function is a user-friendly wrapper that calls 
#'     the \code{rGIG_cpp} function. WE NEED PDF OF GIG(p, a, b).
#' @return A list of outputs.
#' @examples
#'   \dontrun{
#'   rGIG(...)
#'   }

rGIG <- function (p, a, b, seed = 0) 
{
  if (length(a) != length(p)) 
    stop("vector a does not have same length as vector p")
  if (length(a) != length(b)) 
    stop("vector a does not have same length as vector b")
  if (min(a) < 0) 
    stop("vector a must be  pos")
  if (min(b) < 0) 
    stop("vector b must be  pos")
  return(rGIG_cpp(p, a, b, seed))
}
