CompDistFullParIIDK1<-function(modg,R,y,matz){
  #################################################################################
  # For K=1 very particular case, where gllim produces single Gaussians and
  # distances are closed form
  # A version of PreCompDistPar when y and z are iid samples to be compared
  # Case L >1 and FULL covariances 
  # 
  # y and z are now of length DxR
  # modg: result of gllimIIDK1 for K=1, direct gllim parameter theta and Sigmastar
  # precomputed
  # same structure but the mixtures parameters are modified
  # to account for the R observations
  #
  # Compute both the W2 and L2 distances between the GLLiM posterior for y and 
  # the Gllim posterior forall z in matz,  both are Gaussians
  # The specificity is to exploit the specific structure of the GLLiM mixtures
  # to avoid repeating the same computations and accelerate computation using 
  # cores and foreach package. Requires to run before the following cmd: 
  # numCores <- detectCores() #numCores = 8
  # registerDoParallel(numCores)
  #
  # REMARK: W2 reduces to the squared euclidian distance between the means because
  # the covariances matrices are all equal to covstarR. The resulting W2 is 
  # proportional to the distance between empirical means!
  # || Sigma* A^T Sigma^-1 (sum y^r - sum z^r) ||^2
  #
  #################################################################################
  ### Input:
  #  - modg: gllim  IID parameter obtained in the K=1 + Sigmastar (precomputed)
  #  - y and matz are DR x 1 and DR x M observations 
  #  - covstarR is Sigma*
  #  - R is the number of iid replications
  #
  ### Output: a list with the W2 and L2 distance values
  #################################################################################
  
  #### step1 : compute posteriors parameters for y and z (mixtures)
  # notational shortcut: gllim direct parameters
  Aa<-modg$A
  Sigmaa<-modg$Sigma
  ca<-modg$c
  ba<-modg$b
  
  covstarR<-modg$covstarR
  invGammaa<-modg$invGamma
  invSigmaa<-modg$invSigma
  
  
  D= dim(ba)[1]
  L=dim(ca)[1]
 # K<-1
  
  M<-dim(matz)[2]
  
  # y is a DR-dim vector turned into a D x R matrix
  yDR<-matrix(y, ncol=R)
  # matz is a DRxM matrix turned into an array (D,R,M)
  matzDR<-array(matz,c(D,R,M))
  
  postcovy<-covstarR
  
  # for K=1, the trace part is 0 in W2
  # and the b part in the means cancels also because they are the same
  
    # FULL case without the constant part...
  # L x 1  (y)
  postmeany<-t(t(covstarR))%*%t(Aa)%*%invSigmaa%*%t(t(rowSums(yDR)))
  
  #first constant part in L2 : can be made simpler...
  postmean0<-rep(0,L)
  L2cst<-2*L2scal2normal(postmean0,postmean0,postcovy,postcovy);
  
  ### for all other z 
  listdist<-NULL
  # foreach is used to use the cores of the computer)
  
  listdist<- foreach (i=1:M, .combine=cbind) %dopar% {
    postcovz<-covstarR
    # L x 1
    postmeanz<-t(t(covstarR))%*%t(Aa)%*%invSigmaa%*%t(t(rowSums(t(t(matzDR[,,i])))))
    
    MW2dist <- sum((postmeany - postmeanz)^2)
    L2dist <- L2cst-2*L2scal2normal(postmeanz[,1],postmeany[,1],postcovz,postcovy)
    
    c(MW2dist,L2dist) 
  }
}
