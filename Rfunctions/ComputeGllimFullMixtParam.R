ComputeGllimFullMixtParam<-function(modg,z,wthr=0){
  #################################################################################
  # Case L> 1 
  # This is version is for arbitrary Sigmak's, see ComputeGllimMixtParam for the
  # isotropic case (faster)
  # Return the parameters of the Gllim posterior for z 
  # For faster computation when K is large, the mixture components with very small 
  # weights can be removed as they are unlikely to impact the distance values
  #################################################################################
  ### Input:
  #  - modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
  #  - z is  D x 1 observation (vector in one-column matrix format)
  #  - wthr= weight threshold eg 0.001 recommended: we keep the mixture weights 
  #    higher than wthr. Default is 0 meaning no threshold
  ### Output: a list with the Pi, Mu and S parameters
  #################################################################################
  
  #### step1 : compute posteriors parameters for  z (mixture)
  # notational shortcut: gllim direct parameters
  modeleg<-modg$mod
  #Gammaa<-modeleg$Gamma
  Aa<-modeleg$A
  #Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  covstara<-modg$covstar
  invGammaa<-modg$invGamma
  invSigmaa<-modg$invSigma
  
  L<-dim(invGammaa)[1]
  
  # xllim inverse model
  gllimpredz<-gllim_inverse_map(z,modeleg)
  # post mixture weights: dim is 1 x K 
  Piz=t(gllimpredz$alpha)
  leftkz<-seq(1,K)[Piz>wthr]
  dimaz<-sum(Piz>wthr)
  Piz<-Piz[Piz>wthr]
  
  # quantites that depend on z
  #Aks z + bks: 
  postcovz<-array(0,c(L,L,dimaz))
  postmeanz<-matrix(0,L,dimaz)
  
  # covariance matrices are independent on y or z
  postcovz[,,1:dimaz]=covstara[,,leftkz]
  # dim L x K : the mixtures means
  # for z 
  Asybs<-NULL
  for (k in leftkz){
    # ISOTROPIC case
    Asybs <- cbind(Asybs,t(t(covstara[,,k]))%*%(t(Aa[,,k])%*%invSigmaa[,,k]%*%(z-ba[,k]) +invGammaa[,,k]%*%t(t(ca[,k])) ))  
  }
  # L x dimaz  ( z)
  postmeanz=Asybs
  
  list("Mu"=postmeanz, "S"=postcovz, "Pi"=Piz)
  
}
