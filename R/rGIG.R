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