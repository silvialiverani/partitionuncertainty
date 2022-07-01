
Ratio_Rate=function(pz, clusObj=NULL, clustering=NULL){
  n=length(pz[,1])
  p=length(pz[1,])
  if (is.null(clustering)){cl=clusObj$clustering
  }else{cl=clustering}
  scl=cumsum(table(cl))
  scl=n-scl
  RES=data.frame(LRmean=rep(NA,n),LRmax=rep(NA,n),LRmin=rep(NA,n),RATE=rep(NA,n))
  
  for (ip in 1:p){
    qui=which(cl==ip)
    if (length(qui)==1){
      RES$LRmax[qui]=pz[qui,ip]/max(pz[qui,-ip])
      RES$LRmin[qui]=pz[qui,ip]/min(pz[qui,-ip])
      RES$LRmean[qui]=pz[qui,ip]/mean(pz[qui,-ip])
      RES$RATE[qui]=pz[qui,ip]
    }else{
      RES$LRmax[qui]=apply(pz[qui,], MARGIN=1, FUN=function(x){x[ip]/max(x[-ip])})
      RES$LRmin[qui]=apply(pz[qui,], MARGIN=1, FUN=function(x){x[ip]/min(x[-ip])})
      RES$LRmean[qui]=apply(pz[qui,], MARGIN=1, FUN=function(x){x[ip]/mean(x[-ip])})
      RES$RATE[qui]=pz[qui,ip]
    }
  }
  return(RES)    
}


PlotClustersUncertainty=function(pz, clusObj=NULL, clustering=NULL, ratio=1.5){
  #ratio is the ratio between height and width, the width is fxed to 1, by default height=1.5
  n=length(pz[,1])
  p=length(pz[1,])
  if (is.null(clustering)){cl=clusObj$clustering
  }else{cl=clustering}
  scl=cumsum(table(cl))
  scl=n-scl
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  layout(matrix(c(0,2,4,3,1,4),2,3,byrow = TRUE), width=c(0.08,1,0.3), heights=c(0.15,ratio),respect=TRUE)
  #each rect figuring an individual/cluster probability is of size ratio/n times 1/p
  #The individual in cluster 1 with higher probability of belonging to this cluster is located at the first row
  #Then we consider cluster 2, 3 and so on. 
  #Bold black lines are added to figured out the different clusters 
  #and names of the clusters are added in the margin of the table
  
  palette<-colorRampPalette(c("Gray 92","Gray 16"))
   
  #central plot
  par(mar=c(2,0,0,0))
  plot.new()
  plot.window(xlim=c(0,p),ylim=c(0,n),xaxs="i",yaxs='i')
  co=1
  for (ip in 1:(p)){
    qui=which(cl==ip)
    
    if (length(qui)==1){
      PZ=pz[qui,]
      col=palette(100)[findInterval(PZ,seq(0,1,0.01),all.inside=TRUE)]
      rect(xleft=seq(0,p-1,1),ybottom=rep(n-co,p),xright=seq(1,p,1),ytop=rep(n-co+1,p-1),col=col,border=col)
      co=co+1
    }else{
      PZ=pz[qui[order(pz[qui,ip],decreasing=TRUE)],]
      lp=length(PZ[,1])
      for ( i in 1:lp){
        col=palette(100)[findInterval(PZ[i,],seq(0,1,0.01),all.inside=TRUE)]
        rect(xleft=seq(0,p-1,1),ybottom=rep(n-co,p),xright=seq(1,p,1),ytop=rep(n-co+1,p-1),col=col,border=col)
        co=co+1
      }
    }
  }
  for (ip in 1:p){
    lines(c(0,p),c(scl[ip],scl[ip]))
    lines(c(ip,ip),c(0,n))
  }
  rect(xleft=0,xright=p,ybottom=0,ytop=n)
  
  #top margin
  par(mar=c(0,0,2,0))
  plot.new()
  plot.window(xlim=c(0,p),ylim=c(0,1),xaxs="i",yaxs='i')
  rect(xleft=seq(0,p-1,1),ybottom=rep(0,p),xright=seq(1,p,1),ytop=rep(1,p),border=rep(1,p))
  text(x=0.5+0:(p-1),y=rep(0.5,p),labels=c(paste('cluster',seq(1,p))))

  #left margin
  par(mar=c(2,0,0,0))
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,n),xaxs="i",yaxs='i')
  rect(xleft=rep(0,p),ybottom=scl,xright=rep(1,p),ytop=c(scl[2:p],n),border=rep(1,p))
  text(x=0.5,y=(scl+c(n,scl[1:(p)-1]))/2,labels=paste('c.',seq(1,p)))
  
  #legend scale
  par(mar=c(5,2,5,2),las=1)
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,100),xaxs="i",yaxs='i')
  axis(side=4, at=seq(0,100,10),labels=c('0',NA,'0.2',NA,'0.4',NA,'0.6',NA,'0.8',NA,'1'))
  
  for( i in 1:100){
    rect(xleft=0,ybottom=i-1,xright=1,ytop=i,col=palette(100)[i],border=palette(100)[i])
  }
  rect(xleft=0,ybottom=0,ytop=100,xright=1)

  
  par(def.par)  #- reset to default
  return('Plot done')
  
}