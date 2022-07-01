library(PReMiuM)
#library(fields)
library(mvtnorm)


mlgamma=function(alpha,p){
  res=pi^(p*(p-1)/4)
  for (j in 1:p){
    res=res*gamma(alpha+(1-j)/2)
  }
  return(log(res))
}

member_proba=function(runInfoObj,clusObj=NULL, clustering=NULL, X=NULL, file_prefix='output', nskip=9, na.number=-999)
{
  #na.string: when we read the covariate table runInfoObj$xMat the NA may be transformed in number or string, na.string 
  #is the character string corresponding to NAs
  #we add the computation of criterium Cz* which indicates how separated are the clusters in Z*
  
  xModel=NULL
  yModel=NULL
  varSelect=NULL
  nSubjects=NULL
  xMat=NULL
  directoryPath=NULL
  fileStem=NULL
  includeResponse=NULL
  nFixedEffects=NULL
  reportBurnIn=NULL
  nBurn=NULL
  nFilter=NULL
  nSweeps=NULL
  nProgress=NULL
  dPitmanYor=NULL
  nCovariates=NULL
  useNormInvWishPrior=NULL
  
  for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])
  
  
  #clusters
  if (is.null(clustering)){
    diss=clusObj$clustering 
    nk=clusObj$clusterSizes 
  }else{
    diss=clustering
    nk=table(clustering)
  }
  K=length(nk) #nb of Clusters
  
  #sweeps
  if (reportBurnIn==TRUE){
    nskipini=floor(nBurn/nFilter)
  }else{nskipini=0}
  nSweeps=floor(nSweeps/nFilter)
  
  #Data
  if (includeResponse==TRUE){Y=runInfoObj$yMat[,1]}
  if (includeResponse==TRUE){ if ((yModel=='Poisson' | yModel=='Binomial')){YT=runInfoObj$yMat[,2]}} #offset for Poisson, number of trials for binomial
  if (nCovariates>=1){
    Xdata=as.matrix(runInfoObj$xMat)
    Xdata[which(Xdata==na.number)]=NA
    if(is.null(X)){
      X=Xdata
    }else if((dim(X)[2]!=nCovariates)){
      print('WARNING: the vector X has not the dimensions required, it will be ignored')
      X=Xdata
    }
  }
  if (nFixedEffects>=1){W=as.matrix(runInfoObj$wMat)}
  
  #Only usable for Normal mixture with NIWP prior
  if(xModel!='Normal'|useNormInvWishPrior==FALSE)stop('ERROR: This function is only implemented for Normal mixture with normal inverse Wishart prior.') 
  
  #Read Prior hyperparams
  runData<-readLines(file.path(directoryPath,paste(fileStem,'_log.txt',sep='')))
  # set alpha, whether it has been estimated or it is fixed
  # find out if it was fixed or estimated 
  alphaUpdate<-runData[grep('Update alpha',runData)]
  alphaUpdate<-substr(alphaUpdate,regexpr(':',alphaUpdate)+1,nchar(alphaUpdate))
  alphaUpdate<-gsub(' ','',alphaUpdate)
  if (alphaUpdate=="False") {
    alpha<-runData[grep('Fixed alpha: ',runData)]
    alpha<-substr(alpha,regexpr(':',alpha)+1,nchar(alpha))
    alpha<-as.integer(alpha)
  } else {
    # if alpha wasn't fixed, take median value of chain
    firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
    skipLines<-ifelse(reportBurnIn,nBurn/nFilter+1,0)
    lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter		
    alphaFileName <- file(file.path(directoryPath,paste(fileStem,'_alpha.txt',sep='')))
    open(alphaFileName)
    alphaValues<-vector()
    alphaValues[1]<-scan(alphaFileName,what=double(),skip=skipLines,nlines=1,quiet=T)
    for (i in (firstLine+1):lastLine){
      alphaValues[i-firstLine]<-scan(alphaFileName,what=double(),skip=0,nlines=1,quiet=T)
    }
    close(alphaFileName)
    alpha<-median(alphaValues)
  }
  runInfoObj$alphaMPP <- alpha
  if (xModel=="Normal"){
    nu0<-runData[grep('nu0',runData)]
    nu0<-substr(nu0,regexpr(':',nu0)+1,nchar(nu0))
    nu0<-gsub(' ','',nu0)
    nu0<-gsub('\t','',nu0)
    nu0<-as.numeric(nu0)
    
    kappa0<-runData[grep('kappa0',runData)]
    kappa0<-substr(kappa0,regexpr(':',kappa0)+1,nchar(kappa0))
    kappa0<-gsub(' ','',kappa0)
    kappa0<-gsub('\t','',kappa0)
    kappa0<-as.integer(kappa0)
    
    rd<-runData[grep('R0',runData)+1:nCovariates]
    R0=c()
    for (j in 1:nCovariates){
      sub=as.numeric(strsplit(rd[j],split=" ")[[1]])
      sub=sub[which(is.na(sub)==FALSE)]
      R0=c(R0,sub)
    }  
    R0=matrix(R0,ncol=nCovariates,byrow=TRUE)
    
    rd<-runData[grep('mu0',runData)+1:nCovariates]
    mu0=c()
    for (j in 1:nCovariates){
      mu0=c(mu0, as.numeric(rd[j]))
    }
    
    hyperParams<-list(mu0=mu0, nu0=nu0,kappa0=kappa0,R0=R0)
  }
  
  runInfoObj$hyperParams <- hyperParams
  ldetR0=log(det(R0))
  mlgammaKappa0=mlgamma(kappa0/2,nCovariates)
  #MargZ=0.0
   
  #return
  pZ=matrix(0,ncol=K,nrow=length(X[,1]))
  #MARGZ=rep(NA,floor(nSweeps/(nskip+1)))
  #Cz=0
  #Open files
  zFile<-file(file.path(directoryPath,paste(fileStem,'_z.txt',sep='')))
  open(zFile)
  nClustersFile<-file(file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep='')))
  open(nClustersFile)
  PsiFile=file(file.path(directoryPath,paste(fileStem, '_psi.txt',sep='')))
  open(PsiFile)
  if (xModel=="Normal"){muFile<-file(file.path(directoryPath,paste(fileStem,  '_mu.txt',sep='')))
                        open(muFile)
                          SigmaFile<-file(file.path(directoryPath,paste(fileStem, '_Sigma.txt',sep='')))
                        open(SigmaFile)
  }
  if (xModel=="Discrete"){phiFile<-file(file.path(directoryPath,paste(fileStem, '_phi.txt',sep='')))
                          open(phiFile)
  }
  if (includeResponse==TRUE){thetaFile<-file(file.path(directoryPath,paste(fileStem, '_theta.txt',sep='')))
                             open(thetaFile)
                if (nFixedEffects>=1){betaFile<-file(file.path(directoryPath,paste(fileStem, '_beta.txt',sep='')))
                                      open(betaFile)
                }}
  
  #Read files
  currNClusters<-scan(nClustersFile,skip=nskipini,n=1,what=integer(),quiet=T)
  currZ<-scan(zFile,what=integer(),skip=nskipini,n=nSubjects,quiet=T)
  currZ<-1+currZ
  isOpen(zFile)
  currPsi=scan(PsiFile, what=double(), skip=nskipini, n=currNClusters, quiet=T)
  if (xModel=="Normal"){currMu<-scan(muFile,what=double(),skip=nskipini,n=currNClusters*nCovariates,quiet=T)
                          currMu<-array(currMu,dim=c(currNClusters,nCovariates))
                          sigma=scan(SigmaFile,what=double(),skip=nskipini,n=currNClusters*nCovariates*nCovariates,quiet=T)
                          sigma<-array(sigma,dim=c(currNClusters,nCovariates,nCovariates))}
  if (xModel=="Discrete"){currPhi=scan(phiFile,what=double(),skip=nskipini,n=currNClusters*maxNcategories*nCovariates,quiet=T)                     
                          currPhi[which(currPhi==na.number)]=NA
                          currPhi=array(currPhi, dim=c(currNClusters,maxNcategories,nCovariates))}
  if (includeResponse==TRUE){currTheta<-scan(thetaFile,what=double(),skip=nskipini,n=currNClusters,quiet=T)
                if (nFE>=1){currBeta<-scan(betaFile,what=double(),skip=nskipini,n=nFE,quiet=T)}}  
  
  
  if(floor(nSweeps/(nskip+1))<1){print('too large nskip')
  }else{
  for(iit in 1:floor(nSweeps/(nskip+1)))
  {
    print(iit*(nskip+1))

    if (xModel=="Normal"){
      clus=unique(currZ)
      k=length(clus)
      n=length(currZ)
      #MargZ=k*kappa0/2*ldetR0-k*mlgamma(kappa0/2,p=nCovariates) 
      ftilde=matrix(0,ncol=K,nrow=length(X[,1]))
      currftilde=matrix(0,ncol=K,nrow=length(X[,1]))
      for (ic in 1:k){
        #Compute the multivariate t-distribution for each cluster
        c=clus[ic]
        qui=which(currZ==c)
        nc=length(qui)
        nuc=nu0+nc
        df=kappa0+nc-nCovariates+1
        if (nc==1){xc=Xdata[qui,]}else{if (nCovariates==1){xc=mean(Xdata[qui])}else{xc=colMeans(Xdata[qui,])}}
        muc=(nu0*mu0+nc*xc)/(nu0+nc)
        Xsic=R0+nc*nu0/(nc+nu0)*matrix(xc-mu0,ncol=1)%*%matrix(xc-mu0,nrow=1)
        
        Mxc=matrix(xc[1],ncol=1,nrow=nc)
        if(nCovariates>=2){for (j in 2:nCovariates) {Mxc=cbind(Mxc,rep(xc[j],nc))}}
        if(nCovariates>=2){SS=rowSums(apply(Xdata[qui,]-Mxc,MARGIN=1,
                                            FUN=function(x){return(matrix(x,ncol=1)%*%matrix(x,nrow=1))}))
        }else{SS=sum((Xdata[qui,]-Mxc)^2)}
        
        Xsic=Xsic+matrix(SS,byrow=T,ncol=nCovariates)
        Sigma=(nu0+nc+1)/((nu0+nc)*(kappa0+nc-nCovariates+1))*Xsic
        fctilde=matrix(apply(X,MARGIN=1,FUN=dmvt,delta=muc, sigma=Sigma,df=df, log=FALSE),ncol=1)
        for (ik in 1:K){
          ncl=length(which(diss[qui]==ik))
          if (nk[ik]>0){
            currftilde[,ik]=currftilde[,ik]+ncl*fctilde/nk[ik]
          }
        }
        #MargZ=MargZ-nc*nCovariates/2*log(pi)+kappa0/2*ldetR0-(kappa0+nc)/2*log(det(Xsic))+mlgamma((kappa0+nc)/2,nCovariates)-mlgammaKappa0+nCovariates/2*log(nu0/(nu0+nc))
      }
      #Add prior P(Z) (Ewens formula)
      #MargZ=MargZ+k*log(alpha)-sum(log(alpha+seq(0,nSubjects-1)))
      #for (ic in 1:k){
      #  c=clus[ic]
      #  qui=which(currZ==c)
      #  nc=length(qui)
      #  MargZ=MargZ+lfactorial(nc-1)
      #}
      #ftilde=exp(MargZ)*currftilde
      pZ[,1:K]=pZ[,1:K]+currftilde
      #MARGZ[iit]=MargZ
    }
    
    #Compute the criterium Cz*
    #for (i in unique(currZ)){
    #  temp=which(currZ==i)
    #  if (length(unique(diss[temp]))>1){Cz=Cz+length(temp)}
    #}
    
    currNClusters<-scan(nClustersFile,skip=nskip,n=1,what=integer(),quiet=T)
    currZ<-scan(zFile,what=integer(),skip=nskip,n=nSubjects,quiet=T)
    currZ<-1+currZ
    currPsi=scan(PsiFile, what=double(), skip=nskipini, n=currNClusters, quiet=T)
    if (xModel=="Normal"){currMu<-scan(muFile,what=double(),skip=nskip,n=currNClusters*nCovariates,quiet=T)
                            currMu<-array(currMu,dim=c(currNClusters,nCovariates))
                            sigma=scan(SigmaFile,what=double(),skip=nskip,n=currNClusters*nCovariates*nCovariates,quiet=T)
                            sigma<-array(sigma,dim=c(currNClusters,nCovariates,nCovariates))}
    if (xModel=="Discrete"){currPhi=scan(phiFile,what=double(),skip=nskip,n=currNClusters*maxNcategories*nCovariates,quiet=T)                     
                           currPhi[which(currPhi==na.number)]=NA
                           currPhi=array(currPhi, dim=c(currNClusters,maxNcategories,nCovariates))}
    if (includeResponse==TRUE){currTheta<-scan(thetaFile,what=double(),skip=nskip,n=currNClusters,quiet=T)
                  if (nFixedEffects>=1){currBeta<-scan(betaFile,what=double(),skip=nskip,n=nFE,quiet=T)}}
    }
  
  close(zFile)
  close(nClustersFile)
  close(PsiFile)
  if (xModel=="Normal"){close(muFile)
                          close(SigmaFile)}
  if (xModel=="Discrete"){close(phiFile)}
  if (includeResponse==TRUE){close(thetaFile)
  if (nFixedEffects>0){close(betaFile)}}
  
  pZ=pZ/floor(nSweeps/(nskip+1))
  #pZ=pZ/sum(exp(MARGZ))
  fprior=matrix(apply(X,MARGIN=1,FUN=dmvt,delta=mu0, 
                      sigma=(1+nu0)*R0/((kappa0-nCovariates+1)*nu0),
                      df=kappa0-nCovariates+1, log=FALSE),ncol=1)
  pZ[,1:K] = matrix((alpha/(alpha+nSubjects))*fprior,ncol=1)%*%
    matrix(1,ncol=K,nrow=1)+(nSubjects/(alpha+nSubjects))*pZ[,1:K]
  pZ[,1:K] = (matrix(1,ncol=1,nrow=length(X[,1]))%*%matrix(nk/nSubjects,nrow=1))*pZ[,1:K]
  
  #pZ[,K+1]=alpha*fprior/(alpha+nSubjects)
  #pZ[,1:K]=pZ[,1:K]/(alpha+nSubjects)
  #Cz=Cz/(floor(nSweeps/(nskip+1))*(alpha+n))
  
  }
  return(pZ)#list(pZ=pZ,margZ=MARGZ,Cz=Cz))
}#end function
