library(INLA)
library(ngme)
library(Matrix)
graphics.off()
n.lattice <- 30

##
# parameters
## 
theta <- list()

#non-G parameters
theta$type = "NIG" #NIG - normal inverse gaussian, 
                   #GAL - Generalized lapalace distribution
                   #Normal - Gaussian
 

theta$nu <- 1 # shape parameter (not applicable if Normal)
theta$mu <- 1 # assymetrical parameter (not applicable if Normal)

#cov parameters
theta$tau  =  0.2 # precision of field
theta$kappa = 1.5 # range of field

# nugget
theta$sigma_e <- 0.1
##
#create mesh for latent field
##

x=seq(from=0,to=10,length.out=n.lattice)
lattice=inla.mesh.lattice(x=x,y=x)
mesh=inla.mesh.create(lattice=lattice, extend=FALSE, refine=FALSE)


##
# build operator
##
operator_list <- create_operator_matern2D(mesh)
K = theta$tau * (operator_list$G[[1]] + theta$kappa*operator_list$C[[1]])
n <- dim(K)[1]

##
# simulate latent random field
#
##"
if(theta$type == "NIG"){
  V =          rGIG(rep(-0.5, n),
                   rep( theta$nu, n),
                   operator_list$h[[1]]^2 * theta$nu, sample.int(10^6, 1))
}else if(theta$type == "GAL"){
  V = rgamma(n, operator_list$h[[1]] * theta$nu, rep(theta$nu, n)) + 10e-14
}else{
  V =   operator_list$h[[1]]
}
Z <- (- operator_list$h[[1]]  + V) * theta$mu + sqrt(V) * rnorm(n)
X <- as.vector(solve(K, Z)) #the latent field


##
# create random location for observations
# here we shrink the location to remove boundary effect
##
n.obs<- 20
dx <- diff(range(x))
shrink <- 0.3
obs.loc = cbind(runif(n.obs)*dx*(1-2*shrink)  + min(x) + shrink*dx,
                runif(n.obs)*dx*(1-2*shrink)  + min(x) + shrink*dx)


A <-  build.A.matrix(operator_list, locs = obs.loc)

Y <- A%*%X + theta$sigma_e * rnorm(n.obs)

D <- dist(obs.loc)
##
# covariance of Y (without nugget)
##
Ci = Matrix::Diagonal(n,1/operator_list$h[[1]])
Sigma_Y <- as.matrix(A%*%solve(K%*%Ci%*%K,t(A)))
plot(c(as.matrix(D)),c(Sigma_Y),ylab='cov',xlab='dist')
