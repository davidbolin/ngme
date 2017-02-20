library("LDMod")
if(Sys.info()[7] == "jonaswallin"){
  setwd("~/Dropbox/articles/NonGaussianLDA/Thorax data set/Thorax_data_analysis/")
}else{
  setwd('~/Dropbox/research/nongaussian_fields/NonGaussianLDA/Thorax data set/Thorax_data_analysis/')
}
source("estimate.parameters.R")
source("run2case.R")
thorax_data <- read.table("Thorax_data.txt", header = T)
thorax_data$cohort_relev <- relevel(thorax_data$cohort, ref = "1968-01-01")
thorax_data$pib_01       <- ifelse(thorax_data$PIb == 2, 0, 1)

thorax_data$sex <- thorax_data$sex-1
mf <- model.frame(fev1 ~ as.factor(DM2) + age3 + age2 * (cohort_relev + as.factor(pib_01)), thorax_data)
x  <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
y  <- as.matrix(model.extract(mf, "response"))

DM   <- x
RM   <- y
ID   <- thorax_data$id
Time <- matrix(thorax_data$age2)
n.pers <- length(unique(ID))
patids <- unique(ID)
count <- 0
obs_list <- Vin <- Y <- locs <- B_random <- B_fixed <- list()
IDs <- list()
for(i in 1:n.pers){
  index <- ID == patids[i]
  # 285 gigantisk hopp, 10 stort hopp

  if(sum(index) > 1){
    count             <- count + 1
    IDs[[count]]      <- patids[i]
    locs[[count]]     <- Time[index]
    B_random[[count]] <- as.matrix(rep(1, sum(index)))
    B_fixed[[count]]  <- as.matrix(x[index ,2:dim(x)[2]])
    Vin[[count]]      <- rep(1, sum(index))
    Y[[count]]        <- y[index]
  }
}

groups <- group.fixed(B_fixed)

dvars = rep(FALSE,dim(B_fixed[[1]])[2])
for(i in 1:length(B_fixed))
  dvars = dvars + apply(B_fixed[[i]], 2, var)
dvars = dvars==0

ncov = dim(B_fixed[[1]])[2]
covi = matrix(0,length(B_fixed),sum(dvars))
for(i in 1:length(B_fixed)){
  covi[i,] = B_fixed[[i]][1,dvars]
}

clist = list()
for(i in 1:length(groups$groups)){
  clist[[i]] = covi[groups$groups[[i]],]
}
nsamp = 10;
ng = length(groups$groups)
nfree = length(groups$free)
ranks = rep(0,1000)
for(i in 1:1000){
  ind = groups$groups[[sample(ng,1)]] #sample group
  nf = nsamp - length(ind)
  if(nf>0)
    ind = c(ind,sample(nfree,nf)) #sample rest
  #constrct B
  B = 0
  for(j in ind){
    B = B + t(B_fixed[[j]])%*%B_fixed[[j]]
  }
  ranks[i] = rankMatrix(B)[1]
}


res <- estimate.parameters("Normal","Normal","Normal",operator.type = "matern", nIter = 30,
                                pSubsample = 1, subsample.type = 4, save.results = FALSE,standMix=FALSE)
