ComputeCovstarRdetVIso<-function(modg,R){
##################################################################################
# Auxiliary function that goes with Compute1DistMarginIID, to compute the Sigmak*
# and log det Vk in the iid case. These quantites are independent  of y
# Remark: could be included in a new GllimFitIID ...
##################################################################################
### Input:
# - modg is the result of GllimFit
# - R is  the number or iid y
### Output:
# the Sigmastark matrices for the GLLiM-IID and the log of det of Vk matrices
#################################################################################
 
  modeleg<-modg$mod
  Aa<-modeleg$A
  L<-dim(Aa)[2]
  D<-dim(Aa)[1]
  Sigmaa<-modeleg$Sigma
  Gammaa<-modeleg$Gamma
  invGammaa<-modg$invGamma
  # just for initial allocation
  # dim L x L x K
  covstarRa<-invGammaa
  logdetVa<-NULL
  for (k in 1:K){
    covstarRa[,,k]=chol2inv(chol(invGammaa[,,k]+ R*t(Aa[,,k])%*%Aa[,,k]/Sigmaa[1,1,k]))
    
    logdetVa<-c(logdetVa, R*D*log(Sigmaa[1,1,k])+log(det(diag(L)+R*t(Aa[,,k])%*%Aa[,,k]%*%Gammaa[,,k]/Sigmaa[1,1,k])))
  }
  list("covstarR"=covstarRa, "logdetVR"=logdetVa)
}


