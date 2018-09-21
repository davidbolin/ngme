library(ngme)
data(Orthodont)
set.seed(123)
NormalMVD <- ngme(                   fixed       = distance ~ age,
                                     random      = ~ 1|id,
                                     data        = Orthodont,
                                     reffects    = "NIG",
                                     use.process = F,
                                     silent      = F,
                                     nIter       = 2,
                                     controls    = list(estimate.fisher = FALSE,
                                                        subsample.type  = 0,
                                                        nSim  =2,
                                                        nBurnin = 2))
print(NormalMVD$fixed_est)