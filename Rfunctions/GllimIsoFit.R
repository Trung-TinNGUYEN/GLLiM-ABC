GllimIsoFit<-function(thetadata, ydata, K){
  # %%%%%%%% Gaussian Locally Linear Mapping and pre-computed quantities %%%%%%%%%
  # % Description: Fit a GLLiM model using the xLLiM package with constraints
  # % cstr= "i" (isotropic covariance Sigmak), using N associated observations
  # % ydata and parameters thetadata.
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%%% Input %%%%
  # %- thetadata (LxN)        % N training parameters values (L)
  # %- y (DxN)                % N training observations of dimension D
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
  # % -covstar (LxLxK)  % matrices Sigma_k^* covariance matrices of the posterior
  # %                   % mixtures, do not depend on y or z 
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  # xllim result, isotropic covariances by default
  resgllim <- gllim(thetadata, ydata, in_K=K); 
  
  # pre-computation (to save time) indep on y (or z), 
  # same covstar and inverse Gamma matrices for all y
  
  # notation shortcut
  Gammaa<-resgllim$Gamma
  Aa<-resgllim$A
  Sigmaa<-resgllim$Sigma
  
  # just for initial allocation
  # dim L x L x K
  covstara<-Gammaa
  invGammaa<- Gammaa
  
  # inverse of Gammak and Sigmakstar
  for (k in 1:K){
    # invGammaa[,,k]<-solve(Gammaa[,,k]) # "solve" is slower
    # covstara[,,k]=solve(invGammaa[,,k]+t(Aa[,,k])%*%Aa[,,k]/Sigmaa[1,1,k])
    invGammaa[,,k]<-chol2inv(chol(Gammaa[,,k]))
    # Sigmaa[1,1,k] is ok because Sigmak is isotropic in this model
    covstara[,,k]=chol2inv(chol(invGammaa[,,k]+t(Aa[,,k])%*%Aa[,,k]/Sigmaa[1,1,k]))
  }
  
  list("mod"=resgllim, "invGamma"=invGammaa , "covstar"=covstara)
}

# %%%%%%%% Not really used for now: Same as GllimFit with a dopar on k%%%%%%%%%
# could be useful to save time when K is large but does not seem to speed up ??
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GllimIsoFitdoP<-function(thetadata, ydata, K){
  
  # xllim result
  resgllim <- gllim(thetadata, ydata, in_K=K); 
  
  # pre-computation (to save time) indep on y (or z), 
  # same covstar and inverse Gamma matices for all y
  Gammaa<-resgllim$Gamma
  Aa<-resgllim$A
  Sigmaa<-resgllim$Sigma
  
  # just for initial allocation
  # dim L x L x K
  covstara<-Gammaa
  invGammaa<- Gammaa
  
  # inverse of Gammak and Sigmakstar
  # This variant does not seem to speed up??
  foreach (k=1:K) %dopar% {
    invGammaa[,,k]<-chol2inv(chol(Gammaa[,,k]))
    # Sigmaa[1,1,k] is ok because Sigmak is isotropic in this model
    covstara[,,k]=chol2inv(chol(invGammaa[,,k]+t(Aa[,,k])%*%Aa[,,k]/Sigmaa[1,1,k]))
  }
  
  list("mod"=resgllim, "invGamma"=invGammaa , "covstar"=covstara)
}

