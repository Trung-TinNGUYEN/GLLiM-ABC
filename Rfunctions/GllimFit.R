GllimFit<-function(thetadata, ydata, K, constr=list(Sigma="")){
  # %%%%%%%% Gaussian Locally Linear Mapping and pre-computed quantities %%%%%%%%%
  # % Description: Fit a GLLiM model using the xLLiM package with the choice
  # % for the constraints on covariance matrices  Sigmak, 
  # % constr = list(Sigma="d") for diagonal covariances
  # % constr = list(Sigma="i") for iso cov (but fot this one used the ISO version faster)
  # % etc... see xlllim: 'd'=diag., 'i'=iso., '*'=equal for all k, 'v'=equal det. for all k
  # % using N associated observations ydata and parameters thetadata.
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
  # % - invGamma (LxLxK) % inverses of matrices Gammak , do not depend on y or z
  # % - invSigma (D x D x K) % inverses of matrices Sigmak , do not depend on y or z
  # % - covstar (LxLxK)  % matrices Sigma_k^* covariance matrices of the posterior
  # %                    % mixtures, do not depend on y or z 
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  # xllim result 
  resgllim <- gllim(thetadata, ydata, in_K=K,cstr=constr); 
  
  # pre-computation (to save time) of quantities independent on y (or z), 
  # same: covstar and inverse Gamma matrices for all y
  
  # notation shortcut
  Gammaa<-resgllim$Gamma
  Aa<-resgllim$A
  Sigmaa<-resgllim$Sigma
  
  # just for initial allocation
  # dim L x L x K
  covstara<-Gammaa
  invGammaa<- Gammaa
  
  # dim D x D x K
  invSigmaa<- Sigmaa
  
  # inverse of Gammak and Sigmakstar
  for (k in 1:K){
    # invGammaa[,,k]<-solve(Gammaa[,,k]) # "solve" is slower
    # covstara[,,k]=solve(invGammaa[,,k]+t(Aa[,,k])%*%Aa[,,k]/Sigmaa[1,1,k])
    invSigmaa[,,k]<-chol2inv(chol(Sigmaa[,,k]))
    invGammaa[,,k]<-chol2inv(chol(Gammaa[,,k]))
    # Sigmaa[1,1,k] is ok because Sigmak is isotropic in this model
    covstara[,,k]=chol2inv(chol(invGammaa[,,k]+t(Aa[,,k])%*%invSigmaa[,,k]%*%Aa[,,k]))
  }
  
  list("mod"=resgllim, "invGamma"=invGammaa ,"invSigma"=invSigmaa ,  "covstar"=covstara)
}

