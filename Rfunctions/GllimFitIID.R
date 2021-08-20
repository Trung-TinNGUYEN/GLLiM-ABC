GllimFitIID<-function(thetadata, ydata,iR, K, constr=list(Sigma="")){
  # %%%%%%%% Gaussian Locally Linear Mapping and pre-computed quantities %%%%%%%%%
  # % Description: Fit a GLLiM + IID model using the xLLiM package with constraints
  # % cstr= constr on covariance matrices Sigmak, can be
  # constr=list(Sigma="") for FULL covoriances Sigmak (DEFAULT)
  # Remark: this constraint is almost equivalent to bloc diagonal covariances with blocs
  # of size D is y is of size DxR. 
  # constr = list(Sigma="d") for diagonal covariances
  # constr = list(Sigma="i") for iso cov (but fot this one used the ISO version faster)
  #  etc... see xlllim: 'd'=diag., 'i'=iso., '*'=equal for all k, 'v'=equal det. for all k
  # using N associated observations
  # % ydata and parameters thetadata.
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%%% Input %%%%
  # %- thetadata (LxN)        % N training parameters values (L)
  # %- y (DxN)                % N training observations of dimension D
  #### DR x N , R iid replications (R=iR but R is already used)
  #### iR = number of replications
  # %- K (integer)             % number of components in the mixtures
  # %%%% Output %%%%
  # % -mod : the output of the gllim function in xllim that contains
  # %   the parameter Psi: mod$c, etc... 
  # %   - c (LxK)       % Gaussian means of theta
  # %   - Gamma (LxLxK) % Gaussian covariances of theta
  # %   - pi (1xK)      % Gaussian weights 
  # %   - A (DxLxK)     % Affine transformation matrices
  # %   - b (DxK)       % Affine transformation vectors
  # %   - Sigma (DxDxK) % Error covariances assumed FULL here 
  # % -invGamma (LxLxK) % inverses of matrices Gammak , do not depend on y or z
  # % -invSigma DxDxK inverses of Sigmak, do not depend on y or z
  ##### % -covstar (LxLxK)  % matrices Sigma_k^* covariance matrices of the posterior
  ##### the log determinants can be precomputed too
  ##### the Sigmak* and log det Vk in the iid case,quantites independent  of y
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # modified fromm xllim to provide a GLLiM iid estimation of the direct
  # parameters, constraint on covariances can be chosen here
  
  # The case D=1 has to be treated separately due to the 1xM matric turned into
  # a vector by R.. To be merged later
  D=nrow(ydata)/iR
  
  if (D==1) {resgllim <- gllimIID1D(thetadata, ydata,iR, in_K=K,cstr=constr,verb=0)
    # pre-computation (to save time) indep on y (or z), 
  # same covstar and log det and  inverse Gamma matrices for all y
  Aa<-resgllim$A # D x L x K
  L<-dim(Aa)[2]
  #D<-dim(Aa)[1]
  # in case some cluster disappears...
  K<-dim(Aa)[3]
  #K<-dim(Aa)[3]
  Sigmaa<-resgllim$Sigma
  Gammaa<-resgllim$Gamma
  invSigmaa<-Sigmaa
  invGammaa<-Gammaa
  covstarRa<-Gammaa
  # just for initial allocation
  # dim L x L x K
  logdetVa<-NULL
  # inverse of Gammak and Sigmakstar and its log det 
  for (k in 1:K){
    invGammaa[,,k]<-chol2inv(chol(Gammaa[,,k]))
    invSigmaa[,,k]<-chol2inv(chol(Sigmaa[,,k]))
    # attention numerical issue when mat non symmetric before inversion
    # no transpose because D=1 make it a vector anyway
    tempMat<-iR*Aa[,,k]%*%t(Aa[,,k])/invSigmaa[,,k]
    covstarRa[,,k]=chol2inv(chol(invGammaa[,,k]+tempMat))
    #tempMat<-diag(L)+iR*t(Aa[,,k])%*%invSigmaa[,,k]%*%Aa[,,k]%*%Gammaa[,,k]
    ##covstarRa[,,k]=chol2inv(chol(invGammaa[,,k]+ R*t(Aa[,,k])%*%Aa[,,k]/Sigmaa[1,1,k]))
    #covstarRa[,,k]=Gammaa[,,k]%*%chol2inv(chol(tempMat))
    dettemp<-det(diag(L)+tempMat%*%Gammaa[,,k])
    
    # since Sigma is 1x1
    logdetVa<-c(logdetVa, iR*log(Sigmaa[,,k])+log(dettemp))
    #logdetVa<-c(logdetVa, iR*log(det(Sigmaa[,,k]))+log(dettemp))
  } # end for
  } # end if
  else{
  resgllim <- gllimIID(thetadata, ydata,iR, in_K=K,cstr=constr,verb=0)
   
  # pre-computation (to save time) indep on y (or z), 
  # same covstar and log det and  inverse Gamma matrices for all y
  Aa<-resgllim$A
  L<-dim(Aa)[2]
  D<-dim(Aa)[1]
  # in case some cluster disappears...
  K<-dim(Aa)[3]
  #K<-dim(Aa)[3]
  Sigmaa<-resgllim$Sigma
  Gammaa<-resgllim$Gamma
  invSigmaa<-Sigmaa
  invGammaa<-Gammaa
  covstarRa<-Gammaa
  # just for initial allocation
  # dim L x L x K
  logdetVa<-NULL
  # inverse of Gammak and Sigmakstar and its log det 
  for (k in 1:K){
    invGammaa[,,k]<-chol2inv(chol(Gammaa[,,k]))
    invSigmaa[,,k]<-chol2inv(chol(Sigmaa[,,k]))
    # attention numerical issue when mat non symmetric before inversion
    tempMat<-iR*t(Aa[,,k])%*%invSigmaa[,,k]%*%Aa[,,k]
    covstarRa[,,k]=chol2inv(chol(invGammaa[,,k]+tempMat))
    #tempMat<-diag(L)+iR*t(Aa[,,k])%*%invSigmaa[,,k]%*%Aa[,,k]%*%Gammaa[,,k]
    ##covstarRa[,,k]=chol2inv(chol(invGammaa[,,k]+ R*t(Aa[,,k])%*%Aa[,,k]/Sigmaa[1,1,k]))
    #covstarRa[,,k]=Gammaa[,,k]%*%chol2inv(chol(tempMat))
    dettemp<-det(diag(L)+tempMat%*%Gammaa[,,k])
    
    logdetVa<-c(logdetVa, iR*log(det(Sigmaa[,,k]))+log(dettemp))
    
  }
  } # end else
  list("mod"=resgllim, "invGamma"=invGammaa ,"invSigma"=invSigmaa , "covstarR"=covstarRa,"logdetVR"=logdetVa)
  
}

