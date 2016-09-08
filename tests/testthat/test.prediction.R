context("Prediction")
test_that("Subset prediction", {
   n.obs  <- 5
  n <- 10
  pred.type <- "Filter"
  nSim <- 10000
  Y <- locs <- B_random <- B_fixed  <- X <- V <- list()
  for(i in 1:3)
  {
    B_random[[i]] <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
    Y[[i]] <- rep(0,n.obs)
    Y[[i]][i] = 1
    locs[[i]] <- (1:n.obs)*i
    B_fixed[[i]]  <- as.matrix(locs[[i]])
    X[[i]] <- i*rep(0, n)
    V[[i]] <- i*rep(1, n)
  }

  mError_list <- list(noise = "Normal", sigma = 0.1, nu = 1)
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = as.matrix(c(2,-1)),
                            beta_fixed  = as.matrix(c(.1)),
                            Sigma = diag(c(0.1, 0.2)),
                            noise = "Normal",
                            Sigma_epsilon=1)

  operator_list <- create_operator(locs, n, name = "Matern")
  operator_list$kappa <- 2
  operator_list$tau   <- 3


  processes_list = list(noise = "Normal")
  processes_list$X = X
  processes_list$V = V



  res <- predictLong( Y = Y,locs = locs, Brandom.pred = B_random, Bfixed.pred = B_fixed, type = pred.type,
                      nSim = nSim, mixedEffect_list = mixedEffect_list, measurment_list  = mError_list,
                      processes_list = processes_list, operator_list = operator_list,silent = TRUE)

  res2 <- predictLong(pInd = 2, Y = Y,locs = locs, Brandom.pred = B_random, Bfixed.pred = B_fixed, type = pred.type,
                      nSim = nSim, mixedEffect_list = mixedEffect_list, measurment_list  = mError_list,
                      processes_list = processes_list, operator_list = operator_list,silent = TRUE)
  v = sum(abs(res$X.summary[[2]]$Mean - res2$X.summary[[1]]$Mean))
  expect_equal(v,0,tolerance=1e-2)
})


test_that("Subset prediction NIG", {

  n.obs  <- 5
  n <- 10
  pred.type <- "Filter"
  nSim <- 10000
  Y <- locs <- B_random <- B_fixed  <- X <- V <- list()
  for(i in 1:3)
  {
    B_random[[i]] <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
    Y[[i]] <- rep(0,n.obs)
    Y[[i]][i] = 1
    locs[[i]] <- (1:n.obs)*i
    B_fixed[[i]]  <- as.matrix(locs[[i]])
    X[[i]] <- rep(0, n)
    V[[i]] <- rep(1, n)
  }

  mError_list <- list(noise = "Normal", sigma = 0.1, nu = 1)
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = as.matrix(c(2,-1)),
                            beta_fixed  = as.matrix(c(.1)),
                            Sigma = diag(c(0.1, 0.2)),
                            noise = "Normal",
                            Sigma_epsilon=1)

  operator_list <- create_operator(locs, n, name = "Matern",common.grid = FALSE)
  operator_list$kappa <- 2
  operator_list$tau   <- 3



  processes_list = list(noise = "NIG")
  processes_list$X = X
  processes_list$V = V



  res <- predictLong( Y = Y,locs = locs, Brandom.pred = B_random, Bfixed.pred = B_fixed, type = pred.type,
                      nSim = nSim, mixedEffect_list = mixedEffect_list, measurment_list  = mError_list,
                      processes_list = processes_list, operator_list = operator_list,silent=TRUE)

  res2 <- predictLong(pInd = 2, Y = Y,locs = locs, Brandom.pred = B_random, Bfixed.pred = B_fixed, type = pred.type,
                      nSim = nSim, mixedEffect_list = mixedEffect_list, measurment_list  = mError_list,
                      processes_list = processes_list, operator_list = operator_list,silent=TRUE)
  v = sum(abs(res$X.summary[[2]]$Mean - res2$X.summary[[1]]$Mean))
  expect_equal(v,0,tolerance=1e-2)
})


test_that("Normal coverage", {
  n.obs  <- 500
  n <- 100
  n.rep <- 100
  pred.type <- "Filter"
  nSim <- 500
  Y <- locs <- B_random <- B_fixed  <- X <- V <- list()
  for(i in 1:n.rep)
  {
    B_random[[i]] <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
    Y[[i]] <- rep(0,n.obs)
    Y[[i]][i] = 1
    locs[[i]] <- (1:n.obs)*i
    B_fixed[[i]]  <- as.matrix(locs[[i]])
    X[[i]] <- i*rep(0, n)
    V[[i]] <- i*rep(1, n)
  }

  mError_list <- list(noise = "Normal", sigma = 0.1, nu = 1)
  mixedEffect_list  <- list(B_random = B_random,
                            B_fixed  = B_fixed,
                            beta_random = as.matrix(c(2,-1)),
                            beta_fixed  = as.matrix(c(.1)),
                            Sigma = diag(c(0.1, 0.2)),
                            noise = "Normal",
                            Sigma_epsilon=1)

  operator_list <- create_operator(locs, n, name = "Matern")
  operator_list$kappa <- 2
  operator_list$tau   <- 3


  processes_list = list(noise = "Normal")
  processes_list$X = X
  processes_list$V = V


  sim_res <- simulateLongPrior( Y                 = Y,
                                locs              = locs,
                                mixedEffect_list  = mixedEffect_list,
                                measurment_list   = mError_list,
                                processes_list    = processes_list,
                                operator_list     = operator_list)

  processes_list$X <- sim_res$X
  Y = sim_res$Y

  res <- predictLong( Y = Y,locs = locs, Brandom.pred = B_random, Bfixed.pred = B_fixed, type = pred.type,
                      nSim = nSim, mixedEffect_list = mixedEffect_list, measurment_list  = mError_list,
                      processes_list = processes_list, operator_list = operator_list,
                      quantiles = c(0.025,0.975), silent = TRUE)

  covered <- NULL
  for(i in 1:length(locs)){
    covered <- c(covered,(res$Y.summary[[i]]$quantiles[[1]]$field < Y[[i]]) & (res$Y.summary[[i]]$quantiles[[2]]$field> Y[[i]]))
  }

  coverage.mean <- 100*mean(covered)
  coverage.std <- 100*sqrt(var(covered)/sum(n.obs))
  expect_equal(v,0,tolerance=1e-2)
})
