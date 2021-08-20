LogPostVarGllimFull<-function(modg,z){
  ###############################################################################
  # For ABC with posterior log variances as summary statistics: 
  # compute the vector of posterior logvariances in the GLLiM full Sigmak case for an 
  # observation z or a matrix D x M of M observations
  # Remark: similar computation in Compute1Dist but independent
  ###############################################################################
  ###### Input: 
  # - z is D x M
  # - modg is the output of GllimFit, ie a gllim model and precomputed matrices
  ###### output:
  #  - postlog is L x M containing the L posterior log variances for each input vector
  ###############################################################################
  modeleg<-modg$mod
  Aa<-modeleg$A
  Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  covstara<-modg$covstar
  invGammaa<-modg$invGamma
  invSigmaa<-modg$invSigma
  
  # posterior mixture weights from xllim inverse model
  postweightz<-gllim_inverse_map(z,modeleg)$alpha
  # posterior  mixture weights: dim is M x K 
  dimK<-dim(postweightz)[2]
  dimM<-dim(postweightz)[1]
  L<-dim(invGammaa)[1]
  
  # quantites that depend on z
  #Aks y + bks: for all y in z
  postlog<-matrix(0,L,dimM)
  for (ny in 1:dimM){
    # rem: the mixture covariances do not depend on z
    # dim L x K : the mixture means
    Asybs<-ca
    for (k in 1:K){
      # ISOTROPIC case
      #Asybs[,k]= covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(z[,ny]-ba[,k]) +invGammaa[,,k]%*%ca[,k]) 
      Asybs[,k]<-t(t(covstara[,,k]))%*%(t(Aa[,,k])%*%invSigmaa[,,k]%*%(z[,ny]-ba[,k]) + invGammaa[,,k]%*%t(t(ca[,k]))) 
      }
    
    # L x K matrix with posterior logvariances in column
    Firstdiagterm<-(apply(covstara,3,diag) + Asybs^2)%*%postweightz[ny,]
    Seconddiagterm<-(Asybs%*%postweightz[ny,])^2
    
    postlog[,ny]<-log(Firstdiagterm-Seconddiagterm)
  }
  return(postlog)
}
