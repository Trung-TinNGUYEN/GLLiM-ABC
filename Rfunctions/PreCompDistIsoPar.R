PreCompDistIsoPar<-function(modg,y,matz,wthr=0){
  #################################################################################
  # Case L >1
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
  #  - y and matz are D x 1 and D x M observations 
  #  - wthr= weight threshold eg 0.001 recommended: we keep the mixture weights 
  #    higher than wthr. Default is 0 meaning no threshold
  ### Output: a list with the MW2 and L2 distance values
  #################################################################################
  
  
  #### step1 : compute posteriors parameters for y and z (mixtures)
  # notational shortcut: gllim direct parameters
  modeleg<-modg$mod
  #Gammaa<-modeleg$Gamma
  Aa<-modeleg$A
  Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  covstara<-modg$covstar
  invGammaa<-modg$invGamma
  
  L<-dim(invGammaa)[1]
  
  
  # compute inverse model for y and z : a D x M+1 matrix
  matyz<-cbind(y,matz)
  # xllim inverse model
  postweightyz<-gllim_inverse_map(matyz,modeleg)$alpha
  dimMplus1<-dim(postweightyz)[1] #M+1
  
  Piy=t(postweightyz)[,1]
  Kdim<-length(Piy)
  # in case of very low weights
  seuil<-min(wthr,sort(Piy, decreasing=TRUE)[3])
  leftky<-seq(1,Kdim)[Piy>seuil]
  #leftky<-seq(1,Kdim)[Piy>wthr]
  # dimay: number of Gaussians whose weight is above the threshold
  dimay<-sum(Piy>seuil)
  Piy<-Piy[Piy>seuil]
  
  postcovy<-array(0,c(L,L,dimay))
  postmeany<-matrix(0,L,dimay)
  
  postcovy[,,1:dimay]=covstara[,,leftky]
  
  Asybs<-NULL
  for (k in leftky){
    # ISOTROPIC case
    Asybs<-cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(matyz[,1]-ba[,k]) +invGammaa[,,k]%*%ca[,k]) )
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
    leftky<-c(leftky,leftky)
  } 
  
  # Pre-computation of the trace part in the Wassertein distance, does not
  # depend on z or y and could even be pre-computed out of the function like
  # covstara (todo)
  # if computed here depends on y only via leftky 
  # tracecost is leftky x K 
  tracecost<-matrix(0,dimay,Kdim)
  
  for (ii in 1:dimay) {
    sigma1<-postcovy[,,ii]
    for (jj in 1:Kdim){
      sigma2<-covstara[,,jj]
      E2 <- eigen(sigma2)
      V2 <- E2$vectors
      U2 <- solve(V2)
      D2 <- diag(E2$values) 
      sqrt2 <- V2 %*% D2^(1/2) %*% U2
      E <- eigen(sqrt2 %*%sigma1 %*% sqrt2)
      V <- E$vectors
      U <- solve(V)
      D <- diag(E$values) 
      sqrtout <- V %*% D^(1/2) %*% U
      tracecost[ii,jj] <-  sum(diag(sigma1+sigma2 - 2*sqrtout))
    }
  }
  #
  
  #
  # Pre computation of the first quadratic form (dep on y) in the L2 distance
  # could be improved probably (todo)
  gramy <- matrix(NA,nrow=dimay,ncol=dimay) # symmetic
  for (ii in 1:dimay) {
    for (jj in ii:dimay) {
      gramy[ii,jj] <- L2scal2normal(postmeany[,ii],
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
  listdist<- foreach (ny = 2:dimMplus1, .combine=cbind) %dopar% {
  
    Piz=t(postweightyz)[,ny]
    # in case of very low weights
    seuil<-min(wthr,sort(Piz, decreasing=TRUE)[3])
    leftkz<-seq(1,Kdim)[Piz>seuil]
    #leftkz<-seq(1,Kdim)[Piz>wthr]
    dimaz<-sum(Piz>seuil)
    Piz<-Piz[Piz>seuil]
    
    postcovz<-array(0,c(L,L,dimaz))
    postmeanz<-matrix(0,L,dimaz)
    
    postcovz[,,1:dimaz]=covstara[,,leftkz]
    
    Asybs<-NULL
    for (k in leftkz){
      # ISOTROPIC case
      Asybs<-cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(matyz[,ny]-ba[,k]) +invGammaa[,,k]%*%ca[,k]) )
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
    
    L2dist <- L2normalPreComp(mixy,mixz,qgramy)
    
     c(MW2dist,L2dist) 
  }
}
