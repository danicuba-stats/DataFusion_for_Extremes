

run.ExNSDmodelXi.R<-function(nIter=1000,nBurnIn=100,nChains=2,nThin=10,
                         yData,xData,xPred,coordsData,coordsPred,
                        times.yData=NULL,
                        times.xData=NULL,
                        times.xPred=NULL,
                        times.yPred=NULL,
                        basis.type="bspline",
                        basis.dim=5,
                        period=1,
                        phiAlpha=0.1,phiBeta=0.1,
                        aAlpha=2,bAlpha=1,aBeta=2,bBeta=1,aC=2,bC=1,xiY=NULL,xiX=NULL,
                        muD=NULL,SigmaD=NULL,sigmaCTune=10,sigmaDTune = 10,
                        sigmaAlphaPrecInit=1,sigmaBetaPrecInit=1,sigmaYPrecInit=1,sigmaCPrecInit=1,
                        alphaInit=NULL,betaInit=NULL,cInit=NULL,sigmaXPrecInit=1,dInit=NULL,
                        #muAlphaInit=NULL,muBetaInit=NULL,sigmaMuAlpha=1,sigmaMuBeta=1,
                        xiXInit=NULL, xiYInit = NULL, uniqueXi=F, xiXTune = 0.1, xiYTune = 0.1 ,# Initial values for xiX_i and xiY_i. If xiXInit is a single value, assume statX = T (stationary X)
                        xiXb = 2, xiXMu = 0, xiYb = 2, xiYMu = 0,
                        sigmaATune = 1, sigmaBTune = 1, 
                        # lambds.run){
                        # things needed for the reg
                        lambdaInit, lambdaTune=rep(0.01,4), lambdaMu, lambdaSig){

  coordsPred<-matrix(coordsPred,ncol = 2)
  if(!is.list(xPred)&!is.matrix(xPred)&!is.data.frame(xPred)){xPred<-matrix(xPred,ncol=1)}
  if(is.list(yData)){
    if(length(yData)!=nrow(coordsData)){
      stop("'length(yData)' and 'nrow(coordsData)' do not match")
    }
  }else{
    if(is.matrix(yData)|is.data.frame(yData)){
      if(ncol(yData)!=nrow(coordsData)){
        stop("'ncol(yData)' and 'nrow(coordsData)' do not match")
      }
    }else{
      stop("'yData' must be a list, matrix or data frame")
    }
  }
  if(is.list(xData)){
    if(length(xData)!=nrow(coordsData)){
      stop("'length(xData)' and 'nrow(coordsData)' do not match")
    }
  }else{
    if(is.matrix(xData)|is.data.frame(xData)){
      if(ncol(xData)!=nrow(coordsData)){
        stop("'ncol(xData)' and 'nrow(coordsData)' do not match")
      }
    }else{
      stop("'xData' must be a list, matrix or data frame")
    }
  }
  if(is.list(xPred)){
    if(length(xPred)!=nrow(coordsPred)){
      stop("'length(xPred)' and 'nrow(coordsPred)' do not match")
    }
  }else{
    if(is.matrix(xPred)|is.data.frame(xPred)){
      if(ncol(xPred)!=nrow(coordsPred)){
        stop("'ncol(xPred)' and 'nrow(coordsPred)' do not match")
      }
    }else{
      stop("'xPred' must be a list, matrix or data frame")
    }
  }
  
  if(!is.list(times.yData)&inherits(times.yData,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yData' has been provided in a Date format, but will be converted to a numeric vector");times.yData<-lubridate::decimal_date(times.yData)}
  if(!is.list(times.xData)&inherits(times.xData,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xData' has been provided in a Date format, but will be converted to a numeric vector");times.xData<-lubridate::decimal_date(times.xData)}
  if(!is.list(times.xPred)&inherits(times.xPred,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xPred' has been provided in a Date format, but will be converted to a numeric vector");times.xPred<-lubridate::decimal_date(times.xPred)}
  if(!is.list(times.yPred)&inherits(times.yPred,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yPred' has been provided in a Date format, but will be converted to a numeric vector");times.yPred<-lubridate::decimal_date(times.yPred)}
  if(is.list(times.yData)&inherits(times.yData[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yData' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.yData)){times.yData[[i]]<-lubridate::decimal_date(times.yData[[i]])}}
  if(is.list(times.xData)&inherits(times.xData[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xData' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.xData)){times.xData[[i]]<-lubridate::decimal_date(times.xData[[i]])}}
  if(is.list(times.xPred)&inherits(times.xPred[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xPred' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.xPred)){times.xPred[[i]]<-lubridate::decimal_date(times.xPred[[i]])}}
  if(is.list(times.yPred)&inherits(times.yPred[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yPred' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.yPred)){times.yPred[[i]]<-lubridate::decimal_date(times.yPred[[i]])}}
  if(is.null(times.yData)){
    warning("'times.yData' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(yData)|is.data.frame(yData)){
      times.yData<-1:nrow(yData)
    }else{
      if(is.list(yData)){
        times.yData<-vector("list",length(yData))
        for(i in 1:length(yData)){
          times.yData[[i]]<-1:length(yData[[i]])
        }
      }
    }
  }
 
  if(is.null(times.xData)){
    warning("'times.xData' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(xData)|is.data.frame(xData)){
      times.xData<-1:nrow(xData)
    }else{
      if(is.list(xData)){
        times.xData<-vector("list",length(xData))
        for(i in 1:length(xData)){
          times.xData[[i]]<-1:length(xData[[i]])
        }
      }
    }
  }
  
  if(is.null(times.xPred)){
    warning("'times.xPred' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(xPred)|is.data.frame(xPred)){
      times.xPred<-1:nrow(xPred)
    }else{
      if(is.list(xPred)){
        times.xPred<-vector("list",length(xPred))
        for(i in 1:length(xPred)){
          times.xPred[[i]]<-1:length(xPred[[i]])
        }
      }
    }
  }
  if(is.null(times.yPred)){
    if(is.numeric(times.yData)){
      warning("'times.yPred' has not been provided. Defaulting to predictions at 'times.yData'")
      times.yPred<-times.yData
    }else{
      if(is.list(times.yData)){
        times.yPred<-times.yData[[1]]
        warning("'times.yPred' has not been provided. Defaulting to predictions at 'times.yData[[1]]'")
      }
    }
  }
  if(is.matrix(yData)|is.data.frame(yData)){
    yData.mat<-yData
    yData<-matrix.to.list(yData)
  }
  if(is.matrix(xData)|is.data.frame(xData)){
    xData.mat<-xData
    xData<-matrix.to.list(xData)
  }
  if(is.matrix(xPred)|is.data.frame(xPred)){
    xPred.mat<-xPred
    xPred<-matrix.to.list(xPred)
  }

  if(!is.null(times.yData)){
    if(is.matrix(times.yData)|is.data.frame(times.yData)){
      times.yData<-matrix.to.list(times.yData)
    }
    if(is.vector(times.yData)){
      if(exists("yData.mat")){
        times.yData<-times.to.list(yData.mat,times.yData)
      }else{
        stop("'yData' has been provided as a list, while 'times.yData' has been provided as a vector")
      }
    }
  }
  if(!is.null(times.xData)){
    if(is.matrix(times.xData)|is.data.frame(times.xData)){
      times.xData<-matrix.to.list(times.xData)
    }
    if(is.vector(times.xData)){
      if(exists("xData.mat")){
        times.xData<-times.to.list(xData.mat,times.xData)
      }else{
        stop("'xData' has been provided as a list, while 'times.xData' has been provided as a vector")
      }
    }
  }
  if(!is.null(times.xPred)){
    if(is.matrix(times.xPred)|is.data.frame(times.xPred)){
      times.xPred<-matrix.to.list(times.xPred)
    }
    if(is.vector(times.xPred)){
      if(exists("xPred.mat")){
        times.xPred<-times.to.list(xPred.mat,times.xPred)
      }else{
        stop("'xPred' has been provided as a list, while 'times.xPred' has been provided as a vector")
      }
    }
  }

  if(is.matrix(times.yPred)|is.data.frame(times.yPred)|is.list(times.yPred)){
    stop("times.yPred must be a numeric vector")
  }
  if(!is.null(times.yData)&!is.null(times.xData)&!is.null(times.xPred)&!is.null(times.yPred)){
    if(basis.type=="bspline"){
      basis<-fda::create.bspline.basis(rangeval = c(min(c(unlist(times.yData),unlist(times.xData),unlist(times.xPred),times.yPred)),
                                               max(c(unlist(times.yData),unlist(times.xData),unlist(times.xPred),times.yPred))),
                                  nbasis = basis.dim)
    }else{
      if(basis.type=="fourier"){
        basis<-fda::create.fourier.basis(rangeval = c(min(c(unlist(times.yData),unlist(times.xData),unlist(times.xPred),times.yPred)),
                                                 max(c(unlist(times.yData),unlist(times.xData),unlist(times.xPred),times.yPred))),
                                    nbasis = basis.dim,
                                    period = period)
      }else{
        stop("'basis.type' must be either 'bspline' or 'fourier'")
      }
    }
    
    By<-vector("list",length(yData))
    for(i in 1:length(yData)){
      By[[i]]<-fda::eval.basis(times.yData[[i]],basis)
    }
    # browser()
    # browser()
    Bx<-vector("list",length(xData))
    for(i in 1:length(xData)){
      Bx[[i]]<-fda::eval.basis(times.xData[[i]],basis)
    }
    # browser()
    ByPred<-fda::eval.basis(times.yPred,basis)
    BxPred<-vector("list",length(xPred))
    for(i in 1:length(xPred)){
      BxPred[[i]]<-fda::eval.basis(times.xPred[[i]],basis)
    }
  }

  # if(is.null(muAlphaInit)){
  #   muAlphaInit<-matrix(1,nrow=basis.dim,ncol=length(yData))
  # }
  # if(is.null(muBetaInit)){
  #   muBetaInit<-matrix(1,nrow=basis.dim,ncol=length(yData))
  # }
  if(is.null(alphaInit)){
    alphaInit<-matrix(1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(alphaInit)){
    alphaInit<-matrix(1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(betaInit)){
    betaInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(cInit)){
    cInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(dInit)){
    dInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(muD)){
    muD<-rep(0,basis.dim)
  }
  if(is.null(SigmaD)){
    SigmaD<-10*diag(basis.dim)
  }
  if(is.null(xiXInit)){
    xiXInit<-0
  }
  if(is.null(xiYInit)){
    xiYInit<-0
  }
  if(is.list(ByPred)){
    ByPred<-matrix(unlist(ByPred),nrow=nrow(ByPred[[1]]),ncol=ncol(ByPred[[1]]))
  }
  # browser()
  # Other things for the regression part
  yData.ind <- lapply(1:length(yData), function(x){
    return(ifelse(yData.mat[,x]>0,1,0))
  })
  xData.cov <- lapply(xData, function(x){
    xcov <- x
    xcov <- ifelse(xcov>0,1,0)
    tmp <- data.frame(int=rep(1,365),
                     tm1=c(NA,xcov[1:364]),
                     t=c(xcov[1:365]),
                     tp1=c(xcov[2:365],NA))
    # tmp <- data.frame(int = rep(1,length(x)),
    #                   tm1 = c(NA,x[1:(length(x)-1)]), 
    #                   t = x,
    #                   tp1 = c(x[2:length(x)],NA))
    tmp <- as.matrix(tmp)
    return(tmp)
  })
  # xData.cov <- lapply(1:length(xData.cov),function(x){
  #   tmp1 <- xData.cov[[x]]
  #   tmp2 <- yData.mat[,x]
  #   return(tmp1[which(!is.na(tmp2)),])
  # })
  xPredData.cov <- lapply(xPred, function(x){
    xcov <- x
    xcov <- ifelse(xcov>0,1,0)
    tmp <- data.frame(int=rep(1,365),
                      tm1=c(NA,xcov[1:364]),
                      t=c(xcov[1:365]),
                      tp1=c(xcov[2:365],NA))
    tmp <- as.matrix(tmp)
    return(tmp)
  })
  
  ##
  xBasis.ind <- lapply(Bx, function(x){
    apply(x,1,function(m){which(m>0)})
  })
  xBasis.ind <- lapply(xBasis.ind, 
                       function(y){
                         lapply(1:basis.dim, function(x) {
                           lapply(y,function(m) {any(m==x)}) %>% unlist %>% which})})
  
  xBasisPred.ind <- lapply(BxPred, function(x){
    apply(x,1,function(m){which(m>0)})
  })
  xBasisPred.ind <- lapply(xBasisPred.ind, 
                       function(y){
                         lapply(1:basis.dim, function(x) {
                           lapply(y,function(m) {any(m==x)}) %>% unlist %>% which})})
  
  ##
  yBasis.ind <- lapply(By, function(x){
    apply(x,1,function(m){which(m>0)})
  })
  yBasis.ind <- lapply(yBasis.ind, 
                       function(y){
                         lapply(1:basis.dim, function(x) {
                           lapply(y,function(m) {any(m==x)}) %>% unlist %>% which})})
  
  
  m<-ncol(Bx[[1]])
  nData<-nrow(coordsData)
  nPred<-nrow(coordsPred)
  r<-nrow(ByPred)
  output<-vector("list",nChains)
  n.lam <- 4
  for(i in 1:nChains){
# browser()
    output[[i]]<-ExNSDmodel(nIter=nIter, nThin=nThin, yData=yData, xData=xData, 
                            xPred=xPred,coordsData=coordsData,coordsPred=coordsPred,By=By,Bx=Bx, #9
                            ByPred,BxPred,phiAlpha,phiBeta,aAlpha,bAlpha,aBeta,bBeta,#8
                            aC = aC,bC=bC, muD=muD,SigmaD=SigmaD,
                            sigmaAlphaPrecInit = sigmaAlphaPrecInit,sigmaBetaPrecInit = sigmaBetaPrecInit,sigmaCPrecInit = sigmaCPrecInit,#6
                            alphaInit = alphaInit,betaInit = betaInit,cInit=cInit,dInit=dInit,#5
                            sigmaCTune=sigmaCTune, sigmaDTune=sigmaDTune,#2
                            # muAlphaInit=muAlphaInit, muBetaInit=muBetaInit,#2
                            # sigmaMuAlpha=sigmaMuAlpha,sigmaMuBeta=sigmaMuBeta,#2
                            xiXInit = xiXInit, xiYInit = xiYInit, 
                            uniqueXi=uniqueXi, xiXTune = xiXTune, 
                            xiYTune = xiYTune ,# Initial values for xiX_i and xiY_i. If xiXInit is a single value, assume statX = T (stationary X)
                            xiXb = xiXb, xiXMu = xiXMu, xiYb = xiYb, xiYMu = xiYMu,
                            sigmaATune = sigmaATune, sigmaBTune = sigmaBTune,
                            # lambds_run = lambds.run,
                            lambdaInit = lambdaInit, lambdaTune=lambdaTune,
                            lambdaMu = lambdaMu, lambdaSig = lambdaSig, yData_ind = yData.ind,
                            xData_cov = xData.cov, xPredData_cov = xPredData.cov,
                            yBasis_ind=yBasis.ind, xBasis_ind=xBasis.ind, 
                            xBasisPred_ind = xBasisPred.ind)

     if(nBurnIn>0){
      output[[i]]<-output[[i]][-(1:(ceiling(nBurnIn/nThin)+1)),]
    }else{
      output[[i]]<-output[[i]][-1,]
    }
 
    colnames(output[[i]])<-c("sigmaAlphaPrec","sigmaBetaPrec","sigmaCPrec",
                             paste0("alpha[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),
                             paste0("beta[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),
                             paste0("c[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),
                             # "sigmaXPrec",
                             paste0("d[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),
                             # paste0("muAlpha[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),
                             # paste0("muBeta[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),
                             paste0("xiX[",1,"]"),
                             paste0("xiY[",1,"]"),
                             paste0("lambda[",rep(1:n.lam,nData),",",rep(1:nData,each=n.lam),"]"),
                             paste0("alphaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),
                             paste0("betaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),
                             paste0("cPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),
                             paste0("dPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),
                             paste0("lambdaPred[",rep(1:n.lam,nPred),",",rep(1:nPred,each=n.lam),"]"),
                             # paste0("muAlphaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),
                             # paste0("muBetaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),
                             paste0("xiPredX[",1,"]"),
                             paste0("xiPredY[",1,"]"),
                             paste0("yPred[",rep(1:r,nPred),",",rep(1:nPred,each=r),"]"))
    output[[i]]<-coda::mcmc(output[[i]],end = nIter,thin = nThin)
  }
  output<-coda::mcmc.list(output)
  return(output)
}

