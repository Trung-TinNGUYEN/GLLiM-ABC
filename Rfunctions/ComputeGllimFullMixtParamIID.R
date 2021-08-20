ComputeGllimFullMixtParamIID<-function(modg,R,z,wthr=0){
  #################################################################################
  # Case L> 1  and FULL covariance matrices Sigmak
  # Return the parameters of the Gllim IID posterior for z 
  # For faster computation when K is large, the mixture components with very small 
  # weights can be removed as they are unlikely to impact the distance values
  #################################################################################
  ### Input:
  #  - modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
  #  - z is  DR x 1 observation (vector in one-column matrix format)
  #  - wthr= weight threshold eg 0.001 recommended: we keep the mixture weights 
  #    higher than wthr. Default is 0 meaning no threshold
  ### Output: a list with the Pi, Mu and S parameters
  #################################################################################
  
  #### step1 : compute posteriors parameters for  z (mixture)
  # notational shortcut: gllim direct parameters
  modeleg<-modg$mod
  #Gammaa<-modeleg$Gamma
  Aa<-modeleg$A
  Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  pia<-modeleg$pi
  
  covstarR<-modg$covstarR
  logdetVR<-modg$logdetVR
  invGammaa<-modg$invGamma
  invSigmaa<-modg$invSigma
  
  
  L<-dim(invGammaa)[1]
  D= dim(ba)[1]
  #L=dim(ca)[1]
  K<-dim(ba)[2]
  
  # z is a DR-dim vector turned into a D x R matrix
  zDR<-matrix(z, ncol=R)
  
  # Compute post mixture weights for  z: Piz dim is 1 x K 
  # and Pimatz M x K
  logPiz<-NULL
  
  #tmp<-matrix(0,D,K)
  #tmpM<-array(0,c(D,D,K))
  #t(t()) trick to handle L=1 case 
  for(k in 1:K){
    tmpk<-Aa[,,k]%*%t(t(ca[,k])) + ba[,k]
    tmpzDR<-t(rowSums(zDR-tmpk[,1]))
    # tmpM pourrait etre calculer avant
    #tmpMk<-Aa[,,k]%*%t(t(covstarR[,,k]))%*%t(Aa[,,k])/(Sigmaa[1,1,k]^2)
    #logpizk<--0.5*sum((zDR-tmpk[,1])^2/Sigmaa[1,1,k])+0.5*quad.form(tmpMk,t(tmpzDR))-0.5*logdetVR[k]
    tmpMk<-invSigmaa[,,k]%*%Aa[,,k]%*%t(t(covstarR[,,k]))%*%t(Aa[,,k])%*%invSigmaa[,,k]
    logpizk<--0.5*sum(quad.diag(invSigmaa[,,k],(zDR-tmpk[,1]))) + 0.5*quad.form(tmpMk,t(tmpzDR))-0.5*logdetVR[k]
    #
    logPiz<-c(logPiz, log(pia[k])+logpizk)
    #tmp[,k]<-tmpk
    #tmpM[,,k]<-tmpMk
  }
  
  
  # added for very small values? 
  #logPiy<-logPiy-min(logPiy)
  
  den=mylogsumexp(logPiz);
  logPiz= logPiz-den 
  
  Piz=exp(logPiz); 
  
  # in case of very low weights
  seuil<-min(wthr,sort(Piz, decreasing=TRUE)[3])
  leftkz<-seq(1,K)[Piz>seuil]
  dimaz<-sum(Piz>seuil)
  Piz<-Piz[Piz>seuil]
  #Piy<-Piy/sum(Piy)
  
  
  postcovz<-array(0,c(L,L,dimaz))
  postmeanz<-matrix(0,L,dimaz)
  
  postcovz[,,1:dimaz]=covstarR[,,leftkz]
  
  Aszbs<-NULL
  for (k in leftkz){
    # Full Sigma case
    Aszbs<-cbind(Aszbs,t(t(covstarR[,,k]))%*%(t(Aa[,,k])%*%invSigmaa[,,k]%*%t(t(rowSums(zDR-ba[,k]))) + invGammaa[,,k]%*%t(t(ca[,k]))) )
  }
  # L x dimay  (y)
  postmeanz=Aszbs
  
  list("Mu"=postmeanz, "S"=postcovz, "Pi"=Piz)
  
}
