simulate.ngme <- function(object,id=NULL)
{

  if(!is.null(id)){
    id_list <- as.numeric(names(object$Y))
    ind <- which(id_list %in% id)
    object <- crop.lists(object,ind)

  }

  sim_object <- simulateLongPrior(Y                 = object$Y,
                                  locs              = object$locs,
                                  mixedEffect_list  = object$mixedEffect_list,
                                  measurment_list   = object$measurementError_list,
                                  processes_list    = object$processes_list,
                                  operator_list     = object$operator_list)

  a <- attributes(object$Y)
  for(i in 1:length(sim_object)){
    if(names(sim_object)[i]!="U"){
      attributes(sim_object[[names(sim_object)[i]]]) <- a
    }
  }

  return(sim_object)
}
#' @title Simulates prior model using processes and operator
#' 
#' @description 
#' 
#' @param n number of simulations
#' @param process_list procsess list
#' @param operator_list operator list
#' @return A list of simulation
simulate.process <- function(n, operator_list, process_list){
  
  process_list$V <- list()
  process_list$X <- list()
  for(i in 1:length(operator_list$h)){
    process_list$X[[i]] <- rep(0, length(operator_list$h[[i]]))
    process_list$V[[i]] <- operator_list$h[[i]]
  }
  simulateLongProcesses_cpp(list(nsim = n,
                                 operator_list=operator_list,
                                 processes_list=process_list))
}
