PreCompDistL1IsoParIID<-function(modg,R,y,matz,wthr=0){
  #################################################################################
  # A version of PreCompDistPar when y and z are iid samples to be compared
  # Case L = 1 
  # Attempt to implement an approximation of GLLiM for iid data
  # y and z are now of length DxR
  # When R=1 the  PreCompDistPar results are recovered  (Validated on SMA2)
  # modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
  # same structure but the mixtures parameters are modified
  # to account for the R observations
  #
  # Compute both the MW2 and L2 distances between the GLLiM posterior for y and 
  # the Gllim posterior forall z in matz 
  # The specificity is to exploit the specific structure of the GLLiM mixtures
  # to avoid repeating the same computations and accelerate computation using 
  # cores and foreach package. Requires to run before the following cmd: 
  # numCores <- detectCores() #numCores = 8
  # registerDoParallel(numCores)
  # 
  # Rem: For faster computation when K is large, the mixture components with very small 
  # weights can be removed as they are unlikely to impact the distance values
  #################################################################################
  ### Input:
  #  - modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
  #  - y and matz are DR x 1 and DR x M observations 
  #  - covstarR is 
  #  - logdetVR is 
  #  - R is the number of iid replications
  #  - wthr= weight threshold eg 0.001 recommended: we keep the mixture weights 
  #    higher than wthr. Default is 0 meaning no threshold
  ### Output: a list with the MW2 and L2 distance values
  #################################################################################
  
  #### step1 : compute posteriors parameters for y and z (mixtures)
  # notational shortcut: gllim direct parameters
  modeleg<-modg$mod
  Aa<-modeleg$A
  Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  pia<-modeleg$pi
  
  invGammaa<-modg$invGamma
  covstarR<-modg$covstarR
  logdetVR<-modg$logdetVR
  
  
  D= dim(ba)[1]
  L=dim(ca)[1]
  K<-dim(ba)[2]
  
  M<-dim(matz)[2]
  
  # y is a DR-dim vector turned into a D x R matrix
  yDR<-matrix(y, ncol=R)
  # matz is a DRxM matrix turned into an array (D,R,M)
  matzDR<-array(matz,c(D,R,M))
  
  # Compute post mixture weights for y and all the z: Piy dim is 1 x K 
  # and Pimatz M x K
  logPiy<-NULL
  logPimatz<-matrix(0,M,K)
  
  tmp<-matrix(0,D,K)
  tmpM<-array(0,c(D,D,K))
  # t(t()) trick for L=1 case
  for(k in 1:K){
    tmpk<-Aa[,,k]%*%t(t(ca[,k]))+ba[,k]
    tmpyDR<-t(rowSums(yDR-tmpk[,1]))
    # tmpM pourrait etre calculer avant
    tmpMk<-Aa[,,k]%*%t(t(covstarR[,,k]))%*%t(Aa[,,k])/(Sigmaa[1,1,k]^2)
    logpiyk<--0.5*sum((yDR-tmpk[,1])^2/Sigmaa[1,1,k])+0.5*quad.form(tmpMk,t(tmpyDR))-0.5*logdetVR[k]
    logPiy<-c(logPiy, log(pia[k])+logpiyk)
    tmp[,k]<-tmpk
    tmpM[,,k]<-tmpMk
  }
  
  # added for very small values? 
  #logPiy<-logPiy-min(logPiy)
  #postweightyz<-gllimpredyz$alpha
  
  
  den=mylogsumexp(logPiy);
  logPiy= logPiy-den 
  
  Piy=exp(logPiy); 
  
  # in case of very low weights
  seuil<-min(wthr,sort(Piy, decreasing=TRUE)[3])
  leftky<-seq(1,K)[Piy>seuil]
  dimay<-sum(Piy>seuil)
  Piy<-Piy[Piy>seuil]
  #Piy<-Piy/sum(Piy)
  
  
  #for(i in 1:M){
  # Pimatz is  M x K 
  Pimatz<- foreach(i=1:M, .combine=rbind) %dopar% {
    logPiz<-NULL
    # trick just for R=1! 
    matzDRi<-t(t(matzDR[,,i]))
    for (k in 1:K){
      tmpzDR<-t(rowSums(matzDRi-tmp[,k]))
      logpizk<--0.5*sum((matzDRi-tmp[,k])^2/Sigmaa[1,1,k])+0.5*quad.form(tmpM[,,k],t(tmpzDR))-0.5*logdetVR[k]
      logPiz<-c(logPiz, log(pia[k])+logpizk)
      #logPimatz[i,]<-logPiz
    }
    # logPiz 
    den=mylogsumexp(logPiz);
    logPiz= logPiz-den
    exp(logPiz)
  }
  
  # End compute posterir weights
  
  postcovy<-array(0,c(L,L,dimay))
  postmeany<-matrix(0,L,dimay)
  
  postcovy[,,1:dimay]=covstarR[,,leftky]
  
  #L=1 case
  # Pre-computation (faster when L=1) of the trace part in the Wassertein distance, does not
  # depend on z or y and could even be pre-computed out of the function like
  # covstarR (todo)
  # if computed here depends on y only via leftky 
  # tracecost is leftky x K 
  tracecost<-matrix(0,dimay,K)
  for (ii in 1:dimay) {
    sigma1<-postcovy[,,ii]
    for (jj in 1:K){
      sigma2<-covstarR[,,jj]
      tracecost[ii,jj] <-  sigma1+sigma2 - 2*sqrt(sigma1*sigma2)
    }
  }
  #
  
  
  Asybs<-NULL
  for (k in leftky){
    # ISOTROPIC case
    Asybs<-cbind(Asybs,t(t(covstarR[,,k]))%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%t(t(rowSums(yDR-ba[,k]))) + invGammaa[,,k]%*%t(t(ca[,k]))) )
  }
  # L x dimay  (y)
  postmeany=Asybs
  # problem with networkflow when K=1
  # to avoid we duplicate the mixture in this case
  
  if (dimay==1) {
    postmeany<-cbind(postmeany,postmeany)
    temp=array(0,c(L,L,2))
    temp[,,1]<-postcovy
    temp[,,2]<-postcovy
    postcovy<-temp
    Piy=c(0.5,0.5)
    dimay<-2
  } 
  
  #
  # Pre computation of the first quadratic form (dep on y) in the L2 distance
  # could be improved probably (todo)
  gramy <- matrix(NA,nrow=dimay,ncol=dimay) # symmetic
  for (ii in 1:dimay) {
    for (jj in ii:dimay) {
      gramy[ii,jj] <- L2scal2normalL1(postmeany[,ii],
                                    postmeany[,jj],
                                    postcovy[,,ii],
                                    postcovy[,,jj]);
      gramy[jj,ii] <- gramy[ii,jj]
    }
  }
  qgramy<-quad.form(gramy,Piy)
  
  
  ### for all other z 
  listdist<-NULL
  # foreach is used to use the cores of the computer
  
  Piz=Pimatz[1,]
  seuil<-min(wthr,sort(Piz, decreasing=TRUE)[3])
  leftkz<-seq(1,K)[Piz>seuil]
  dimaz<-sum(Piz>seuil)
  
  listdist<- foreach (i=1:M, .combine=cbind) %dopar% {
    
    Piz=Pimatz[i,]
    seuil<-min(wthr,sort(Piz, decreasing=TRUE)[3])
    leftkz<-seq(1,K)[Piz>seuil]
    dimaz<-sum(Piz>seuil)
    Piz<-Piz[Piz>seuil]
    
    postcovz<-array(0,c(L,L,dimaz))
    postmeanz<-matrix(0,L,dimaz)
    
    postcovz[,,1:dimaz]=covstarR[,,leftkz]
    
    Asybs<-NULL
    for (k in leftkz){
      # ISOTROPIC case
      Asybs<-cbind(Asybs,t(t(covstarR[,,k]))%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%t(t(rowSums(t(t(matzDR[,,i]))-ba[,k]))) + invGammaa[,,k]%*%t(t(ca[,k])) ) )
    }
    # L x dimaz  (y)
    postmeanz=Asybs
    
    if (dimaz==1) {
      postmeanz<-cbind(postmeanz,postmeanz)
      temp=array(0,c(L,L,2))
      temp[,,1]<-postcovz
      temp[,,2]<-postcovz
      postcovz<-temp
      Piz=c(0.5,0.5)
      dimaz<-2
      leftkz<-c(leftkz,leftkz)
    } 
    
    mixy<-list("Mu"=postmeany, "S"=postcovy, "Pi"=Piy)
    mixz<-list("Mu"=postmeanz, "S"=postcovz, "Pi"=Piz)
    
    # tracecost has to be reduced to leftky x leftkz before the call
    #tracecostz<-matrix(0,dimay,dimaz)
    tracecostz<-tracecost[,leftkz]
    
    MW2dist <- was2mixPreComp(mixy,mixz,tracecostz)
    L2dist <- L2normalPreCompL1(mixy,mixz,qgramy)
    
    c(MW2dist,L2dist) 
  }
}
