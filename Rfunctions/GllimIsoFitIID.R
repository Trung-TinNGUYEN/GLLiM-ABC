GllimIsoFitIID<-function(thetadata, ydata,iR, K){
  # %%%%%%%% Gaussian Locally Linear Mapping and pre-computed quantities %%%%%%%%%
  # % Description: Fit a GLLiM + IID model using the xLLiM package with constraints
  # % cstr= "i" (isotropic covariance Sigmak), using N associated observations
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
  # %   - Sigma (DxDxK) % Error covariances assumed isotropic here ('i'=iso option)
  # % -invGamma (LxLxK) % inverses of matrices Gammak , do not depend on y or z
  ##### % -covstar (LxLxK)  % matrices Sigma_k^* covariance matrices of the posterior
  ##### the log determinants can be precomputed too
  ##### the Sigmak* and log det Vk in the iid case,quantites independent  of y
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # modified fromm xllim to provide a GLLiM iid estimation of the direct
  # parameters, isotropic covariances by default
  resgllim <- gllimIID(thetadata, ydata,iR, in_K=K, verb=0); 
  
  # pre-computation (to save time) indep on y (or z), 
  # same covstar and log det and  inverse Gamma matrices for all y
  Aa<-resgllim$A
  L<-dim(Aa)[2]
  D<-dim(Aa)[1]
  #K<-dim(Aa)[3]
  Sigmaa<-resgllim$Sigma
  Gammaa<-resgllim$Gamma
  invGammaa<-Gammaa
  covstarRa<-Gammaa
  # just for initial allocation
  # dim L x L x K
  logdetVa<-NULL
  # inverse of Gammak and Sigmakstar and its log det 
  for (k in 1:K){
    invGammaa[,,k]<-chol2inv(chol(Gammaa[,,k]))
    #covstarRa[,,k]=chol2inv(chol(invGammaa[,,k]+ R*t(Aa[,,k])%*%Aa[,,k]/Sigmaa[1,1,k]))
    covstarRa[,,k]=Sigmaa[1,1,k]*chol2inv(chol(Sigmaa[1,1,k]*invGammaa[,,k]+ iR*t(Aa[,,k])%*%Aa[,,k]))
    
    logdetVa<-c(logdetVa, iR*D*log(Sigmaa[1,1,k])+log(det(diag(L)+iR*t(Aa[,,k])%*%Aa[,,k]%*%Gammaa[,,k]/Sigmaa[1,1,k])))
  }
  list("mod"=resgllim, "invGamma"=invGammaa , "covstarR"=covstarRa,"logdetVR"=logdetVa)

}

