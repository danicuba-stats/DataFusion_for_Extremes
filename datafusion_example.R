################################################################################
### Data fusion for extreme prediction over the Greater London Area
### Author: Daniela Cuba (dcuba@agricarbon.co.uk)
### Date: September 20th, 2024
################################################################################

# Load libraries
library(evd)
library(psych)
library(fda)
library(MASS)
library(mvtnorm)
library(Rcpp)
library(POT)
library(NSD)
library(dplyr)
library(scales)
library(fBasics)
library(boot)

# Source functions
setwd("D:/DataFusion_Extremes")
source("./run_ENSDmodel_ddGPD_xiXxiY_p_lambda_v2.R")
sourceCpp("./ExNSD_cpp_code_lambda_v2.cpp")
load("./london_data.RData")

# Prep EAC4 data - matrix form with every location a different column and censored at the 80th quantile
coords.Pred <- expand.grid(pm2.5mean_all$longs,pm2.5mean_all$lats)
x.pred.All <- matrix(0, ncol=220, nrow=365)
for(i in 1:220){
  row <- ((i-1)%/%11)
  cols <- i-(row*11)
  tmp <- pm2.5mean_all$data[(row+1),cols,((3652-364):(3652))]
  tmp.q80 <- quantile(tmp,0.8,na.rm=T)
  tmp <- tmp - tmp.q80
  tmp[which(tmp<0)] <- 0
  x.pred.All[,i] <- tmp
}

# Prep indices
m <- 60
nData <- 12
nPred <- 10
r <- 365
aurn.dates <- 1:365
d.pred.ind <- (6 + 4 * m * nData + 2 * m * nPred + 4 * nData):(5 + 4 * m * nData + 3 * m * nPred + 4 * nData)
c.pred.ind <- (6 + 4 * m * nData + 3 * m * nPred + 4 * nData):(5 + 4 * m * nData + 4 * m * nPred + 4 * nData)
xiX.pred.ind <- (6 + 4 * m * nData + 4 * m * nPred + 4 * nData)
xiY.pred.ind <- (7 + 4 * m * nData + 4 * m * nPred + 4 * nData)
lam.pred.ind <- (8 + 4 * m * nData + 4 * m * nPred + 4 * nData): (7 + 4 * m * nData + 4 * m * nPred + 4 * nData + 4 * nPred)
y.pred.ind <-  (8 + 4 * m * nData + 4 * m * nPred + 4 * nData + 4 * nPred):(7 + 4 * m * nData + 4 * m * nPred + r*nPred + 4 * nData + 4 * nPred)

# Set up initial values (and others)
laminits <- lapply(1:12, function(x){
  xcov <- ifelse(RSdata[,x]>0,1,0)
  tmp.df <- data.frame(int=rep(1,365),
                       tm1=c(NA,xcov[1:364]),
                       t=c(xcov[1:365]),
                       tp1=c(xcov[2:365],NA),
                       ISind=ifelse(ISdata[,x]>0,1,0))
  tmp <- glm(ISind ~ tm1 + t + tp1, data = tmp.df, family = "binomial")
  return(tmp$coefficients)
})

lamMus <- logit(c(0.05, 0.7,0.9,0.7)) # Rough mean of coefficients across locations
laminits <- do.call(rbind,laminits)
aurn.dates <- 1:365
pred.ints <- seq(1,220,by=10)

# Save predictions only
save_preds_list <- list()
# Run model (save every 10)
for(i in 1:length(pred.ints)){
  i.preds <- pred.ints[i] + 0:9
  start <- Sys.time()
  model_i <- run.ExNSDmodelXi.R(nIter = 500000, yData = ISdata, xData = RSdata,
                                xPred = x.pred.All[,i.preds],
                                coordsData = as.matrix(coords.IS),
                                coordsPred = as.matrix(coords.Pred[i.preds,]),
                                times.yData = aurn.dates,
                                times.xData = aurn.dates,
                                times.xPred = aurn.dates,
                                times.yPred = aurn.dates,
                                basis.type = "bspline", basis.dim = m,
                                aAlpha=2,bAlpha=10,aBeta=2,bBeta=10, # improves when bs are 20
                                aC=2,bC=20,
                                sigmaCTune= 0.85, sigmaDTune = 0.5,
                                sigmaATune = 0.3, sigmaBTune = 0.1,
                                xiXInit=-0.3, xiYInit=-0.2,
                                uniqueXi=0, xiXTune = 0.5, xiYTune = 0.5,
                                xiXb = 0.05, xiXMu = 0, xiYb = 0.05, xiYMu = 0,
                                phiAlpha=1.6,phiBeta=1.6, nChains=1,nThin =20, nBurnIn=300000,
                                lambdaInit = laminits, lambdaTune=c(2,rep(1,3)),
                                lambdaMu = lamMus, lambdaSig = c(1,rep(1,3)))
  end <- Sys.time()
  print(paste0("Finished with i = ",i, " at ",Sys.time(),". It took ",end - start))
  save_preds_list[[i]] <- model_i[[1]][,y.pred.ind]
  save(model_i,
       file=paste0("./output/model_",i,".RData"))
  
}

### Saving only the prediction of each location for mapping
save(save_preds_list, file="./output/map_data.RData")

