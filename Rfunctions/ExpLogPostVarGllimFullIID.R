ExpLogPostVarGllimFullIID<-function(modg,matz,R){
  ###############################################################################
  # For ABC with posterior expectations and/or log variances as summary statistics: 
  # compute the vector of posterior expectations and
  # posterior logvariances in the GLLiM IID FULL covariance case for an 
  # observation matz or a matrix DR x M of M observations of size DR where R is the
  # number of iid repetitions
  # Remark: similar computation in PreCompDistxxxx but independent
  ###############################################################################
  ###### Input: 
  # - R is the number of replications
  # - z is DR x M
  # - modg is the output of GllimIsoFitIID, ie a gllim IID model and precomputed matrices
  ###### output:
  #  - postexp is L x M containing the L posterior expectations for each input vector
  #  - postlog is L x M containing the L posterior log variances for each input vector
  ###############################################################################
  modeleg<-modg$mod
  Aa<-modeleg$A
  Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  pia<-modeleg$pi
  
  covstarR<-modg$covstarR
  logdetVR<-modg$logdetVR
  invGammaa<-modg$invGamma
  invSigmaa<-modg$invSigma
  
  # posterior mixture weights from GLLiM IID inverse model
  D= dim(ba)[1]
  L=dim(ca)[1]
  K<-dim(ba)[2]
  M<-dim(matz)[2]
  
  # matz is a DRxM matrix turned into an array (D,R,M)
  matzDR<-array(matz,c(D,R,M))
  
  # Compute post exp and log var for all the z:  M x 2L 
  ###logPimatz<-matrix(0,M,K)
  
  ExpLogVarmatz<-matrix(0,M, 2*L)
  
  # compute common quantites before
  tmp<-matrix(0,D,K)
  tmpM<-array(0,c(D,D,K))
  # t(t()) trick for L=1 case
  for(k in 1:K){
    tmp[,k]<-Aa[,,k]%*%t(t(ca[,k]))+ba[,k]
    tmpM[,,k]<-invSigmaa[,,k]%*%Aa[,,k]%*%t(t(covstarR[,,k]))%*%t(Aa[,,k])%*%invSigmaa[,,k]
  }
  
  #for(i in 1:M){
  # ExpLogVarmatz is  M x 2L  
  ExpLogVarmatz<- foreach(i=1:M, .combine=rbind) %dopar% {
    logPiz<-NULL
    Asybs<-NULL
    # trick just for R=1! 
    matzDRi<-t(t(matzDR[,,i]))
    for (k in 1:K){
      # unormalized log weights ( K)
      tmpzDR<-t(rowSums(matzDRi-tmp[,k]))
      logpizk<--0.5*sum(quad.diag(invSigmaa[,,k],(matzDRi-tmp[,k]))) +0.5*quad.form(tmpM[,,k],t(tmpzDR))-0.5*logdetVR[k]
      #logpizk<--0.5*sum((matzDRi-tmp[,k])^2/Sigmaa[1,1,k])+0.5*quad.form(tmpM[,,k],t(tmpzDR))-0.5*logdetVR[k]
      logPiz<-c(logPiz, log(pia[k])+logpizk)
      #logPimatz[i,]<-logPiz
      # for the means (L x K)
      #Asybs<-cbind(Asybs,t(t(covstarR[,,k]))%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%t(t(rowSums(matzDRi-ba[,k]))) + invGammaa[,,k]%*%t(t(ca[,k])) ) )
      Asybs<-cbind(Asybs,t(t(covstarR[,,k]))%*%(t(Aa[,,k])%*%invSigmaa[,,k]%*%t(t(rowSums(matzDRi-ba[,k]))) + invGammaa[,,k]%*%t(t(ca[,k]))) )
      }
    # logPiz 
    den=mylogsumexp(logPiz);
    logPiz= logPiz-den
    # post weights K 
    Piz<-exp(logPiz)
    postE<-Asybs%*%Piz # L x K x K = L x 1
    # for post var L x K with posterior logvariances in column
    Firstdiagterm<-(apply(covstarR,3,diag) + Asybs^2)%*%Piz
    Seconddiagterm<-postE^2
    postV<-log(Firstdiagterm-Seconddiagterm)
    # combine mean and var
    c(postE, postV)
  }
  # End  ##compute posterior weights in Pimatz M x K
  
  return(ExpLogVarmatz)
}
