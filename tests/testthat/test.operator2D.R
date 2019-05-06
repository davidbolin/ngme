library(testthat)
library(ngme)
library(INLA)
context("2Doperator")


set.operator <- function(theta,rho,tau1,tau2,kappa1,kappa2,G,C){
  D <- matrix(0,2,2)
  D[1,1] = cos(theta) + rho*sin(theta)
  D[1,2] = -sin(theta)*sqrt(1+rho^2)
  D[2,1] = sin(theta) - rho*cos(theta)
  D[2,2] = cos(theta)*sqrt(1+rho^2)
  
  K1 = (tau1/kappa1)*operator_list$G[[1]] + tau1*kappa1*operator_list$C[[1]]
  K2 = (tau2/kappa2)*operator_list$G[[1]] + tau2*kappa2*operator_list$C[[1]]
  
  return(rbind(cbind(D[1,1]*K1,D[1,2]*K2),
               cbind(D[2,1]*K1,D[2,2]*K2)))
  
}

like <- function(K,X,iV){
  KX = K%*%matrix(X,length(X),1)
  return(log(det(K)) - 0.5*t(KX)%*%iV%*%KX)
}

test_that("2Doperator_gradient", {

noise="Gaussian"
kappa1 = 1
kappa2 = 1
tau1 = 5
tau2 = 5
rho = 1
theta = 1
sigma.e = c(0.01,0.1)
beta.fixed = c(0,0)

n.lattice = 3
n.rep = 1 #number of replicates

#create mesh 
x=seq(from=0,to=10,length.out=n.lattice)
lattice=inla.mesh.lattice(x=x,y=x)
mesh=inla.mesh.create(lattice=lattice, extend=FALSE, refine=FALSE)


operator_list <- create_operator_matern2Dbivariate(mesh)

operator_list$kappa1 <- kappa1
operator_list$kappa2 <- kappa2
operator_list$tau1   <- tau1
operator_list$tau2   <- tau2
operator_list$rho   <- rho
operator_list$theta   <- theta

processes_list = list(noise = "Normal", V <- list())
for(i in 1:n.rep){
  processes_list$V[[i]] <- c(operator_list$h[[1]],operator_list$h[[1]])  
  processes_list$X[[i]] <- (1:length(operator_list$h[[1]]))
}

Ci = diag(1/operator_list$h[[1]])
grad = test_2Doperator(processes_list, operator_list)

K <- set.operator(theta,rho,tau1,tau2,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])

l <- like(K,processes_list$X[[1]],Ci)
e <- 0.00001

K.eps <- set.operator(theta,rho+e,tau1,tau2,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le <- like(K.eps,processes_list$X[[1]],Ci)
drho = (le-l)/e
expect_equal(as.vector(drho),grad$drho,tolerance=0.1)

K.eps <- set.operator(theta,rho-e,tau1,tau2,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le2 <- like(K.eps,processes_list$X[[1]],Ci)
d2rho = (le-2*l+le2)/e^2
expect_equal(as.vector(d2rho),grad$d2rho,tolerance=0.1)

K.eps <- set.operator(theta,rho,tau1+e,tau2,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le <- like(K.eps,processes_list$X[[1]],Ci)
dtau1 = (le-l)/e
expect_equal(as.vector(dtau1),grad$dtau1,tolerance=0.1)


K.eps <- set.operator(theta,rho,tau1-e,tau2,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le2 <- like(K.eps,processes_list$X[[1]],Ci)
d2tau1 = (le-2*l+le2)/e^2
expect_equal(as.vector(d2tau1),grad$d2tau1,tolerance=0.1)

K.eps <- set.operator(theta,rho,tau1,tau2+e,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le <- like(K.eps,processes_list$X[[1]],Ci)
dtau2 = (le-l)/e
expect_equal(as.vector(dtau2),grad$dtau2,tolerance=0.1)

K.eps <- set.operator(theta,rho,tau1,tau2-e,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le2 <- like(K.eps,processes_list$X[[1]],Ci)
d2tau2 = (le-2*l+le2)/e^2
expect_equal(as.vector(d2tau2),grad$d2tau2,tolerance=0.1)

K.eps <- set.operator(theta,rho,tau1,tau2,kappa1+e,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le <- like(K.eps,processes_list$X[[1]],Ci)
dkappa1 = (le-l)/e
expect_equal(as.vector(dkappa1),grad$dkappa1,tolerance=0.1)

K.eps <- set.operator(theta,rho,tau1,tau2,kappa1-e,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le2 <- like(K.eps,processes_list$X[[1]],Ci)
d2kappa1 = (le-2*l+le2)/e^2
expect_equal(as.vector(d2kappa1),grad$d2kappa1,tolerance=0.1)

K.eps <- set.operator(theta,rho,tau1,tau2,kappa1,kappa2+e,operator_list$G[[1]],operator_list$C[[1]])
le <- like(K.eps,processes_list$X[[1]],Ci)
dkappa2 = (le-l)/e

expect_equal(as.vector(dkappa2),grad$dkappa2,tolerance=0.1)

K.eps <- set.operator(theta,rho,tau1,tau2,kappa1,kappa2-e,operator_list$G[[1]],operator_list$C[[1]])
le2 <- like(K.eps,processes_list$X[[1]],Ci)
d2kappa2 = (le-2*l+le2)/e^2
expect_equal(as.vector(d2kappa2),grad$d2kappa2,tolerance=0.1)


K.eps <- set.operator(theta+e,rho,tau1,tau2,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le <- like(K.eps,processes_list$X[[1]],Ci)
dtheta = (le-l)/e

expect_equal(as.vector(dtheta),grad$dtheta,tolerance=0.1)

K.eps <- set.operator(theta-e,rho,tau1,tau2,kappa1,kappa2,operator_list$G[[1]],operator_list$C[[1]])
le2 <- like(K.eps,processes_list$X[[1]],Ci)
d2theta = (le-2*l+le2)/e^2
expect_equal(as.vector(d2theta),grad$d2theta,tolerance=0.1)

})