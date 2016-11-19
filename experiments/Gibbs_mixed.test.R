V_0 <- LDMod::rGIG(-1.5, nu.mixed, nu.mixed, seed ) 
V1    <- matrix(0, nrow = sim + 1, ncol = 1)
V1[1] <- V_0
U1    <- matrix(0, nrow= sim, ncol=  2)
p_GIG = 0.5 * (-1 - dim(B_random)[2])
eSigma <- eigen(Sigma) 
Sigma_inv <- eSigma$vectors%*%diag(1/eSigma$values)%*%t(eSigma$vectors)
R         <- eSigma$vectors%*%diag(sqrt(1/eSigma$values))%*%t(eSigma$vectors)
R_sigma   <- eSigma$vectors%*%diag(sqrt(eSigma$values))%*%t(eSigma$vectors)

b_U <- (t(B_random)%*%(Y - B_random%*%beta_random))/sigma_eps^2
Q_U <- (t(B_random)%*%B_random)/sigma_eps^2

###
# Gibbs part
#
###
sampleU <- function(V)
{
  Q   =  Sigma_inv/V
  #b_U <- (t(B_random)%*%Y)/sigma_eps^2
  #Q_U <- (t(B_random)%*%B_random)/sigma_eps^2
  b   = b_U + Q%*%(-mu.mixed + V*mu.mixed)
  return(rmvnorm(n = 1, mean = solve(Q_U + Q, b), sigma = solve(Q + Q_U)))
}

sampleV <- function(U)
{
  
  #sample V
  U_ <- U + mu.mixed #- beta_random
  b = t(U_)%*%Sigma_inv%*%(U_ )
  b = b + nu.mixed
  a_GIG = t(mu.mixed)%*%Sigma_inv%*%mu.mixed + nu.mixed
  return(LDMod::rGIG(p_GIG, a_GIG, b, ceiling(runif(1, 0, 10^12))))  
}