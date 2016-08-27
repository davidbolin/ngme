# simulate data from the prior model
# @param Y only used to get size of objects
simulateLongPrior <- function( Y,
                               locs,
                               mixedEffect_list,
                               measurment_list,
                               processes_list,
                               operator_list)
{
  obs_list <- list()
  for(i in 1:length(locs))
    obs_list[[i]] <- list(A = spde.A(locs[[i]],
                                     operator_list$loc,
                                     right.boundary = operator_list$right.boundary,
                                     left.boundary = operator_list$left.boundary),
                      Y=Y[[i]],
                      locs = locs[[i]])


  input <- list( obs_list = obs_list,
                 operator_list = operator_list,
                 measurment_list = measurment_list,
                 mixedEffect_list = mixedEffect_list,
                 processes_list = processes_list)

  output <- simulateLongGH_cpp(input)
  return(output)
}