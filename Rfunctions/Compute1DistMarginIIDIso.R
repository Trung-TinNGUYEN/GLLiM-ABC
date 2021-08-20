Compute1DistMarginIIDIso<-function(modg, covstarR,logdetVR,R, y,z,imarg=0,wthr=0){
##############################################################################
# Attempt to implement an approximation of GLLiM for iid data
# y and z are now of length DxR
# When R=1 the Compute1DistMargin results are recovered  (Validated on SMA2)
# modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
# same structure as Compute1DistMargin but the mixtures parameters are modified
# to account for the R observations
# TO BE CONTINUED and further tested: results do not seem to improve for now
##############################################################################
 
  # y is a DR-dim vector turned into a D x R matrix
  yDR<-matrix(y, ncol=R)
  zDR<-matrix(z, ncol=R)
  
  #### step1 : compute posteriors parameters for y and z (mixtures)
  # notational shortcut: gllim direct parameters
  modeleg<-modg$mod
  #Gammaa<-modeleg$Gamma
  Aa<-modeleg$A
  Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  pia<-modeleg$pi
  #covstara<-covstarR
  invGammaa<-modg$invGamma
  #detVa<-modg$detV
  
  # D= dim(ba)[1]
  L=dim(ca)[1]
  if(imarg==0) imarg<-seq(1,L)
  
  # post mixture weights: dim is 2 x K 
  
  logPiy<-NULL
  logPiz<-NULL
  
  for(k in 1:K){
    tmp<-Aa[,,k]%*%ca[,k]+ba[,k]
    tmpyDR<-t(rowSums(yDR-tmp[,1]))
    
    # tmpM pourrait etre calculer avant
    tmpM<-Aa[,,k]%*%covstarR[,,k]%*%t(Aa[,,k])/(Sigmaa[1,1,k]^2)
    
    logpiyk<--0.5*sum((yDR-tmp[,1])^2/Sigmaa[1,1,k])+0.5*quad.form(tmpM,t(tmpyDR))-0.5*logdetVR[k]
    logPiy<-c(logPiy, log(pia[k])+logpiyk)
    
    tmpzDR<-t(rowSums(zDR-tmp[,1]))
    logpizk<--0.5*sum((zDR-tmp[,1])^2/Sigmaa[1,1,k])+0.5*quad.form(tmpM,t(tmpzDR))-0.5*logdetVR[k]
    logPiz<-c(logPiz, log(pia[k])+logpizk)
  }
  
  
  # added for very small values? 
  #logPiy<-logPiy-min(logPiy)
  #postweightyz<-gllimpredyz$alpha
  den=mylogsumexp(logPiy);
  logPiy= logPiy-den 
  Piy=exp(logPiy); 
  leftky<-seq(1,K)[Piy>wthr]
  dimay<-sum(Piy>wthr)
  Piy<-Piy[Piy>wthr]
  #Piy<-Piy/sum(Piy)
  
  den=mylogsumexp(logPiz);
  logPiz= logPiz-den
  Piz=exp(logPiz); 
  leftkz<-seq(1,K)[Piz>wthr]
  dimaz<-sum(Piz>wthr)
  Piz<-Piz[Piz>wthr]
  # Piz<-Piz/sum(Piz)
  
  # quantites that depend on y (z)
  #Aks y + bks: 2 because y and z 
  postcovy<-array(0,c(L,L,dimay))
  postmeany<-matrix(0,L,dimay)
  postcovz<-array(0,c(L,L,dimaz))
  postmeanz<-matrix(0,L,dimaz)
  
  
  # cov dep pas de y or z
  postcovy[,,1:dimay]=covstarR[,,leftky]
  postcovz[,,1:dimaz]=covstarR[,,leftkz]
  # dim L x K : the mixtures means
  # for y 
  Asybs<-NULL
  for (k in leftky){
    # ISOTROPIC case
    Asybs<-cbind(Asybs,covstarR[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%t(t(rowSums(yDR-ba[,k]))) + invGammaa[,,k]%*%ca[,k]) )
  }
  
# L x dimay  (y)
  postmeany=Asybs
  
  # for z 
  Asybs<-NULL
  for (k in leftkz){
    # ISOTROPIC case
    Asybs <- cbind(Asybs,covstarR[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%t(t(rowSums(zDR-ba[,k]))) +invGammaa[,,k]%*%ca[,k])) 
  }
  # L x dimaz  ( z)
  postmeanz=Asybs
  
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
  if (dimaz==1) {
    postmeanz<-cbind(postmeanz,postmeanz)
    temp=array(0,c(L,L,2))
    temp[,,1]<-postcovz
    temp[,,2]<-postcovz
    postcovz<-temp
    Piz=c(0.5,0.5)
    dimaz<-2
  } 
  
  ## here we compute the mixture for margin in imarg
  MW2dist<-NULL
  L2dist<-NULL
  for(im in imarg){
    postmeanymarg<-postmeany[im,]
    postcovymarg<-postcovy[im,im,]
    postmeanzmarg<-postmeanz[im,]
    postcovzmarg<-postcovz[im,im,]
    
    # trick for dim 1
    postcovymarg<-array(postcovymarg, dim = c(1,1,dimay))
    postcovzmarg<-array(postcovzmarg, dim = c(1,1,dimaz))
    postmeanymarg<-array(postmeanymarg, dim = c(1,dimay))
    postmeanzmarg<-array(postmeanzmarg, dim = c(1,dimaz))
    
    mixy<-list("Mu"=postmeanymarg, "S"=postcovymarg, "Pi"=Piy)
    mixz<-list("Mu"=postmeanzmarg, "S"=postcovzmarg, "Pi"=Piz)
    
    #### step2 : Compute MW2 and L2 between imarg posterior margin given y and post given z
    MW2dist <- c(MW2dist, was2mixL1(mixy,mixz))
    L2dist <- c(L2dist, L2normalL1(mixy,mixz))
  }
  
  # end for
  c(MW2dist,L2dist)
}
