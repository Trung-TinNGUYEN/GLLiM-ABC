Compute1DistL1Iso<-function(modg,y,z,wthr=0){
  #################################################################################
  # Case L = 1 
  # and one distance at a time    (for use with doPar or parallelisation)
  # Compute both the MW2 and L2 distances between the GLLiM posterior for y and 
  # the Gllim posterior for z (1 distance at a time in contrast to other
  # ComputeDist function.
  # For faster computation when K is large, the mixture components with very small 
  # weights can be removed as they are unlikely to impact the distance values
  #################################################################################
  ### Input:
  #  - modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
  #  - y and z are each D x 1 observations (vector in one-column matrix format)
  #  - wthr= weight threshold eg 0.001 recommended: we keep the mixture weights 
  #    higher than wthr. Default is 0 meaning no threshold
  ### Output: a list with the MW2 and L2 distance values
  #################################################################################
  
  L=1
  #### step1 : compute posteriors parameters for y and z (mixtures)
  # notational shortcut: gllim direct parameters
  modeleg<-modg$mod
  Gammaa<-modeleg$Gamma
  Aa<-modeleg$A
  Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  covstara<-modg$covstar
  invGammaa<-modg$invGamma
  
  # compute inverse model for y and z : a D x 2 matrix
  matyz<-cbind(y,z)
  # xllim inverse model
  gllimpredyz<-gllim_inverse_map(matyz,modeleg)
  # post mixture weights: dim is 2 x K 
  postweightyz<-gllimpredyz$alpha
  Piy=t(postweightyz)[,1]
  leftky<-seq(1,K)[Piy>wthr]
  dimay<-sum(Piy>wthr)
  Piy<-Piy[Piy>wthr]
  
  Piz=t(postweightyz)[,2]
  leftkz<-seq(1,K)[Piz>wthr]
  dimaz<-sum(Piz>wthr)
  Piz<-Piz[Piz>wthr]
  
  # quantites that depend on y (z)
  #Aks y + bks: 2 because y and z 
  postcovy<-array(0,c(L,L,dimay))
  postmeany<-matrix(0,L,dimay)
  postcovz<-array(0,c(L,L,dimaz))
  postmeanz<-matrix(0,L,dimaz)
  
  
  # cov dep pas de y or z
  postcovy[,,1:dimay]=covstara[,,leftky]
  postcovz[,,1:dimaz]=covstara[,,leftkz]
  # dim L x K : the mixtures means
  # for y 
  Asybs<-NULL
  for (k in leftky){
    # ISOTROPIC case
    Asybs<-cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(matyz[,1]-ba[,k]) +invGammaa[,,k]%*%ca[,k]) )
  }
  # L x dimay  (y)
  postmeany=Asybs
  
  # for z 
  Asybs<-NULL
  for (k in leftkz){
    # ISOTROPIC case
    Asybs <- cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(matyz[,2]-ba[,k]) +invGammaa[,,k]%*%ca[,k])) 
  }
  # L x dimaz  ( z)
  postmeanz=Asybs
  
  # trick for dim 1
  postcovy<-array(postcovy, dim = dim(postcovy))
  postcovz<-array(postcovz, dim = dim(postcovz))
  
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
  
  mixy<-list("Mu"=postmeany, "S"=postcovy, "Pi"=Piy)
  mixz<-list("Mu"=postmeanz, "S"=postcovz, "Pi"=Piz)
  
  
  #### step2 : Compute MW2 and L2 between postvgiven y and post given z
  MW2dist<-was2mixL1(mixy,mixz)
  L2dist<- L2normalL1(mixy,mixz)
  
  c(MW2dist, L2dist)
  #list("MW2"=MW2dist, "L2"=L2dist)
}
