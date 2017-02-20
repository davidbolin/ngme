standardize.covariates <- function(B.list)
{
  if(!is.null(B.list)){
    if(0){
      n <- length(B.list)
      n.cov <- dim(B.list[[1]])[2]
      n.obs <- sum(unlist(lapply(B.list,length))/n.cov)
      m <- lapply(B.list,colMeans)
      #means <- rowMeans(matrix(unlist(m,use.names = FALSE),n.cov,n))
      #ss = lapply(B.list,function(x) colSums(t((t(x)-means)^2)))
      ss = lapply(B.list,function(x) colSums(x^2))
      stds <- sqrt(rowSums(matrix(unlist(ss,use.names = FALSE),n.cov,n))/n.obs)

      B.out <- B.list
      ss <- rep(0,n.cov)
      for(i in 1:n){
        B.out[[i]] <- t(t(B.list[[i]])/stds)
      }
      return(list(B = B.out,M = diag(stds), Minv = diag(1/stds)))
    } else {
      n.cov <- dim(B.list[[1]])[2]
      if(n.cov>1){
        BB = lapply(B.list,function(x) t(x)%*%x)
        Bs = Reduce("+",BB)/length(B.list)
        R = chol(Bs)
        Ri = solve(R)
        B.out = lapply(B.list,function(x) x%*%Ri)
      } else {
        return(list(B  =B.list, M = matrix(1), Minv = matrix(1)))
      }

      return(list(B  =B.out, M = R, Minv = Ri))
    }
  } else {
    return(NULL)
  }
}
scale.beta <- function(beta,B.list,inv = FALSE)
{
  if(!is.matrix(beta))
    beta = as.matrix(beta)

  nc = dim(B.list$M)[1]
  if(inv){
    if(dim(beta)[1]==nc)
      return(B.list$Minv%*%beta)
    else
      return(t(B.list$Minv%*%t(beta)))
  } else {
    if(dim(beta)[1]==nc)
      return(B.list$M%*%beta)
    else
      return(t(B.list$M%*%t(beta)))
  }
}

scale.sigma <- function(sigma,B.list,inv = FALSE)
{
  if(inv){
    return(B.list$Minv%*%sigma%*%t(B.list$Minv))
  } else {
    return(B.list$M%*%sigma%*%t(B.list$M))
  }
}

#group individuals based on Fixed effect matrix B
#function returns a list groups and a vector free.
#each element in groups contains a vector of elements sufficient to make t(B)*B full rank
#the elments in free can be sampled freely given that one of the vectors in groups is sampled
group.fixed <- function(B)
{
  ncov = dim(B[[1]])[2]
  nper = length(B)
  BB = list()
  grouped = rep(0,nper)
  for(i in 1:nper){
    BB[[i]] = t(B[[i]])%*%B[[i]]
    if(rankMatrix(BB[[i]])[1] == ncov)
      grouped[i] = -1
  }
  groups <- Bgroup <- list()
  groupnumber = 1
  while(sum(abs(grouped)>0)<nper)
  {
    ng = which(grouped==0)
    gi = ng[1]
    Bi = BB[[ng[1]]]
    crank = rankMatrix(Bi)[1]

    grouped[ng[1]] = groupnumber

    if(length(ng)==1){
      if(crank < ncov)
        grouped[ng[1]] = 0
      break
    } else {
      ng = ng[2:length(ng)]
    }

    while(crank<ncov)
    {
      Bp = Bi + BB[[ng[1]]]
      prank = rankMatrix(Bp)[1]
      if(prank>crank){
        gi = c(gi, ng[1])
        grouped[ng[1]] = groupnumber
        Bi = Bp
        crank = prank
      }
      if(length(ng)==1){
        break
      } else {
        ng = ng[2:length(ng)]
      }
    }
    if(length(ng)==1){
      if(crank < ncov)
        grouped[grouped==groupnumber] = 0
      break
    }
    groups[[groupnumber]] = gi-1
    Bgroup[[groupnumber]] = Bi
    groupnumber = groupnumber +1
  }
  free.samples = which(grouped<=0)
  return(list(groups = groups,
              free = free.samples-1))
}
