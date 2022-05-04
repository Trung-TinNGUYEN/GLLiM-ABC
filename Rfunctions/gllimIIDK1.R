gllimIIDK1 = function(tapp,yapp,iR,in_r=NULL,maxiter=150,Lw=0,cstr=NULL,verb=0,in_theta=NULL){
  # Adapted by F. Forbes from standard GLLiM
  # Remark: hybrid case not done 
  # MaximizationIID Faster than first (old) version
  # %%%%%%%% General EM Algorithm for Gaussian Locally Linear Mapping %%%%%%%%%
  # %%% Author: Antoine Deleforge (April 2013) - antoine.deleforge@inria.fr %%%
  # % Description: Compute maximum likelihood parameters theta and posterior
  # % probabilities r=p(z_n=k|x_n,y_n;theta) of a gllim model with constraints
  # % cstr using N associated observations t and y.
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%%% Input %%%%
  # %- t (LtxN)               % Training latent variables
  # %- y (DRxN)                % Training observed variables: the R iid D-dim are stacked
  # DR x N , R iid replications (R=iR but R is already used)
  # iR = number of replications
  # %- in_K (int)             % Initial number of components
  # % <Optional>
  # %- Lw (int)               % Dimensionality of hidden components (default 0)
  # %- maxiter (int)          % Maximum number of iterations (default 100)
  # %- in_theta (struct)      % Initial parameters (default NULL)
  # %                         | same structure as output theta
  # %- in_r (NxK)             % Initial assignments (default NULL)
  # %- cstr (struct)          % Constraints on parameters theta (default NULL,'')
  # %   - cstr$ct             % fixed value (LtxK) or ''=uncons.
  # %   - cstr$cw             % fixed value (LwxK) or ''=fixed to 0
  # %   - cstr$Gammat         % fixed value (LtxLtxK) or ''=uncons.
  # %                         | or {'','d','i'}{'','*','v'} [1]
  # %   - cstr$Gammaw         % fixed value (LwxLwxK) or ''=fixed to I
  # %   - cstr$pi             % fixed value (1xK) or ''=uncons. or '*'=equal
  # %   - cstr$A             % fixed value (DxL) or ''=uncons.
  # %   - cstr$b             % fixed value (DxK) or ''=uncons.
  # %   - cstr$Sigma         % fixed value (DxDxK) or ''=uncons.
  # %                         | or {'','d','i'}{'','*'} [1]
  # %- verb {0,1,2}           % Verbosity (default 1)
  # %%%% Output %%%%
  # %- theta  (struct)        % Estimated parameters (L=Lt+Lw)
  # %   - theta.c (LxK)       % Gaussian means of X
  # %   - theta.Gamma (LxLxK) % Gaussian covariances of X
  # %   - theta.pi (1xK)      % Gaussian weights of X
  # %   - theta.A (DxLxK)     % Affine transformation matrices
  # %   - theta.b (DxK)       % Affine transformation vectors
  # %   - theta.Sigma (DxDxK) % Error covariances
  # %- r (NxK)                % Posterior probabilities p(z_n=k|x_n,y_n;theta) 
  # %%% [1] 'd'=diag., 'i'=iso., '*'=equal for all k, 'v'=equal det. for all k
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # % ======================Input Parameters Retrieval=========================
  # A faire plus tard mais pas forcement indispensable maintenant? 
  # [Lw, maxiter, in_theta, in_r, cstr, verb] = ...
  #     process_options(varargin,'Lw',0,'maxiter',100,'in_theta',NULL,...
  #                              'in_r',NULL,'cstr',struct(),'verb',1);
  
  
  #### Example: toto<-gllimIIDK1(nlmparams, nlmdata,100, cstr=list(Sigma=""),verb=0)
  ####  provides c, Gamma, A, b, Sigma
  
  # % ==========================Default Constraints============================
  if(! "ct" %in% names(cstr)) cstr$ct=NULL;
  if(! "cw" %in% names(cstr)) cstr$cw=NULL;
  if(! "Gammat" %in% names(cstr)) cstr$Gammat=NULL;
  if(! "Gammaw" %in% names(cstr)) cstr$Gammaw=NULL;
  if(! "pi" %in% names(cstr)) cstr$pi=NULL;
  if(! "A" %in% names(cstr)) cstr$A=NULL;
  if(! "b" %in% names(cstr)) cstr$b=NULL;
  if(! "Sigma" %in% names(cstr)) cstr$Sigma="i";
  
  if (ncol(tapp) != ncol(yapp)) {stop("Observations must be in columns and variables in rows")}
  
  # iid adptation
  D= nrow(yapp)/iR
  N = ncol(yapp)  
  
  # yapp is a DR x N matrix turned into an array (D,N,R)
  yappR<-array(0,c(D,N,iR))
  for (i in 1:iR){
    yappR[,,i]<-yapp[(1+D*(i-1)):(D*i),]
  }
  yapp1<-yappR[,,1] # for emgmIID and init
  
  ######## K=1 ##########
  K=1
  
  rk<-matrix(1,1,N)
  
  
  MaximizationIIDK1 = function(tapp,yapp,yappR,iR,muw,Sw,cstr,verb){
    # yapp is a DR x N matrix turned into an array yappR (D,N,R)
    if(verb>=1) print('  M'); 
    if(verb>=3) print(' k='); 
    # rem: when K=1, then r=1 and ncol(r) =NULL
    K = 1;
    D = nrow(yapp)/iR;
    N=ncol(yapp)
    Lt = nrow(tapp)
    Lw = ifelse(is.null(muw),0,nrow(muw))
    L=Lt+Lw;
    
    th = list()
    th$c=matrix(NaN,nrow=L,ncol=K)
    th$Gamma=array(0,dim=c(L,L,K));
    
    #if(Lw>0)
    #   {th$c[(Lt+1):L,]=cstr$cw; #% LwxK
    #   th$Gamma[(Lt+1):L,(Lt+1):L,]=cstr$Gammaw;} #% LwxLwxK}
    
    #th$pi=rep(NaN,K);    
    th$A=array(NaN,dim=c(D,L,K));
    th$b=matrix(NaN,nrow=D,ncol=K);
    th$Sigma= array(NaN,dim=c(D,D,K));  
    
    
    
    # IID case special here
    #sumth = list()
    sumthA=array(0,dim=c(D,L,K));
    sumthb=matrix(0,nrow=D,ncol=K);
    sumthSigma= array(0,dim=c(D,D,K));  
    
    # for each replication we can apply the usual gllim and then take the
    # mean of the estimations for bk, Sigmak. the ck, Gammak parameters do not change
    # Attention Ak requires some care, also a mean in the end...
    
    # xbar_k and ybar_k can be computed before, they do not vary with the replication
   # rk_bar=rep(0,K);
    rk_bar=N
    # rk = matrix(0,N,K)
    yk_bar<-matrix(0,D,K)
    xk_bar<-matrix(0,Lt,K)
    
    for (k in 1:K){
      #rk=r[,k]; #% 1xN        
      #rk_bar[k]=sum(rk); #% 1x1
      yk_bartemp<-0
      # trick to be improved to deal with case D=1 remove and put in gllimIID1D for now
      #if(D==1){for (iirep in 1:iR){yk_bartemp<-yk_bartemp + rowSums(sweep(t(yappR[,,iirep]),2,rk,"*"))/rk_bar[k]}
      #        }
      #else{
      for (iirep in 1:iR){yk_bartemp<-yk_bartemp + rowSums(sweep(yappR[,,iirep],2,rk,"*"))/N}
      #  }
      yk_bar[,k]<-yk_bartemp/iR
      #% Dx1
      xk_bar[,k]<-rowSums(sweep(tapp,2,rk,"*"))/N # L x 1
      
      ### M steps that do not depend on iR, ie pi_k, ck, Gammak
      if(verb>=3) print(k);    
      #  % Posteriors' sums
      ##rk=r[,k]; #% 1xN         
      ##rk_bar[k]=sum(rk); #% 1x1
      
      #debut if
      if(Lt>0)
      {
        if(verb>=3) {print('c');}
        #% Compute optimal mean ctk  
        if(is.null(cstr$ct))
        {th$c[1:Lt,k]=rowSums(sweep(tapp,2,rk,"*"))/N;}# % Ltx1
        else {th$c[1:Lt,k]=cstr$ct[,k];}
        #% Compute optimal covariance matrix Gammatk
        if(verb>=3) {print('Gt');}
        diffGamma= sweep(sweep(tapp,1,th$c[1:Lt,k],"-"),2,sqrt(rk),"*");    #% LtxN
        if( is.null(cstr$Gammat) || (length(cstr$Gammat)==1 & cstr$Gammat=='*')) # | ou ||?
          # %%%% Full Gammat
        {th$Gamma[1:Lt,1:Lt,k]=tcrossprod(diffGamma)/N; #% DxD
        }#th$Gamma[1:Lt,1:Lt,k]=th$Gamma[1:Lt,1:Lt,k];                    
        else
        {
          if( !is.character(cstr$Gammat))
            #%%%% Fixed Gammat
          {th$Gamma[1:Lt,1:Lt,k]=cstr$Gammat[,,k];  }          
          else
          {
            if(cstr$Gammat[1]=='d' | cstr$Gammat[1]=='i')
              #% Diagonal terms   
            {gamma2=rowSums(diffGamma^2)/N; #%Ltx1
            if(cstr$Gammat[1]=='d')
              #%%% Diagonal Gamma
            {th$Gamma[1:Lt,1:Lt,k]=diag(gamma2);} #% LtxLt  
            else
              #%%% Isotropic Gamma
            {th$Gamma[1:Lt,1:Lt,k]=mean(gamma2)*diag(Lt);} #% LtxLt
            }
            else
            {if(cstr$Gammat[1]=='v')
              #%%%% Full Gamma
            {th$Gamma[1:Lt,1:Lt,k]=tcrossprod(diffGamma)/N;} #% LtxLt
              else {# cstr$Gammat,
                stop('  ERROR: invalid constraint on Gamma.'); }
            }
          }
        }				
      } # fin if       
      
      # % Compute optimal weight pik
      #th$pi[k]=rk_bar[k]/N; #% 1x1
    }
    
    ### M steps that depend on the yappR
    for (irep in 1:iR){
      #trick for D=1 removed for now
      #if(D==1){yappi<-t(yappR[,,irep])} 
      #else{
      yappi<-yappR[,,irep]
      #}
      # not used but kept for the else # here x = tapp
      if(Lw>0)
      {x=rbind(tapp,muw[,,k]); #% LxN
      Skx=rbind(cbind(matrix(0,Lt,Lt),matrix(0,Lt,Lw)),cbind(matrix(0,Lw,Lt),Sw[,,k])); }#% LxL    
      else
      {x=tapp; #% LtxN
      Skx=matrix(0,Lt,Lt);} #%LtxLt
      # end if else
      for (k in 1:K){
       # rk=r[,k]
        if(verb>=3) {print('A');}
        # if else
        if(is.null(cstr$b))
        {# % Compute weighted means of y and x
          #### These means are changing in the iid case
          #### yk_bar is now the mean over replications (iR) 
          #### This part has been factored out for speed...
          if(L==0)
          {xk_bar[,k]=NULL;}
        }
        else
        {yk_bar[,k]=cstr$b[,k];
        xk_bar[,k]=rep(0,L);
        th$b[,k]=cstr$b[,k]; 
        } # end if else
        
        #% Compute weighted, mean centered y and x
        weights=matrix(1,1,N); #% 1xN        
        y_stark=sweep(yappi,1,yk_bar[,k],"-"); #% DxN #col or row? 
        y_stark= sweep(y_stark,2,weights,"*"); #% DxN  #col or row?     
        if(L>0)
        { x_stark=sweep(tapp,1,xk_bar[,k],"-"); #% LxN  
        x_stark= sweep(x_stark,2,weights,"*"); #% LxN
        }            
        else
        {x_stark=NULL;}
        
        #% Robustly compute optimal transformation matrix Ak
        #warning off MATLAB:nearlySingularMatrix;
        if(!all(Skx==0))
        {if(N>=L & det(Skx+tcrossprod(x_stark)) >10^(-8))
        {th$A[,,k]=tcrossprod(y_stark,x_stark) %*% qr.solve(Skx+tcrossprod(x_stark));} #% DxL
          else
          {th$A[,,k]=tcrossprod(y_stark,x_stark) %*% ginv(Skx+tcrossprod(x_stark));} #%DxL
        }
        else
        {if(!all(x_stark==0))
        {if(N>=L & det(tcrossprod(x_stark))>10^(-8))
        {th$A[,,k]=tcrossprod(y_stark,x_stark) %*% qr.solve(tcrossprod(x_stark));} #% DxL
          else
          {if(N<L && det(crossprod(x_stark))>10^(-8)) 
          {th$A[,,k]=y_stark %*% solve(crossprod(x_stark)) %*% t(x_stark);} #% DxL
            else
            {if(verb>=3) print('p') 
              th$A[,,k]=y_stark %*% ginv(x_stark);}  #% DxL
          }}
          else
          {#% Correspond to null variance in cluster k or L=0:
            if(verb>=1 & L>0) print('null var\n');
            th$A[,,k]=0; # % DxL
          }
        } 
        # end else
        
        if(verb>=3)print('b'); 
        # % Intermediate variable wk=y-Ak*x
        if(L>0)
        {wk=yappi-th$A[,,k]%*%x;} #% DxN  
        else
        {wk=yappi;}
        
        #% Compute optimal transformation vector bk
        if(is.null(cstr$b))
          th$b[,k]=rowSums(sweep(wk,2,rk,"*"))/N; #% Dx1 
        
        if(verb>=3) print('S');
        #% Compute optimal covariance matrix Sigmak
        
        # if(Lw>0)
        #   { Awk=th$A[,(Lt+1):L,k];
        #   Swk=Sw[,,k];                
        #  ASAwk=Awk%*%tcrossprod(Swk,Awk);}
        #else
        ASAwk=0;
        
        diffSigma=sweep(sweep(wk,1,th$b[,k],"-"),2,sqrt(rk),"*"); #%DxN
        
        if (cstr$Sigma %in% c("","*")) 
        {#%%%% Full Sigma  
          th$Sigma[,,k]=tcrossprod(diffSigma)/N; #% DxD
          th$Sigma[,,k]=th$Sigma[,,k]+ASAwk;  }                  
        else 
        {
          if(!is.character(cstr$Sigma))
            #%%%% Fixed Sigma
          {th$Sigma=cstr$Sigma;}
          else {
            if(cstr$Sigma[1]=='d' || cstr$Sigma[1]=='i')
              #% Diagonal terms   
            {sigma2=rowSums(diffSigma^2)/N; #%Dx1
            if(cstr$Sigma[1]=='d')
            {#%%% Diagonal Sigma
              th$Sigma[,,k]=diag(sigma2,ncol=D,nrow=D); #% DxD
              if (is.null(dim(ASAwk))) {th$Sigma[,,k]=th$Sigma[,,k] + diag(ASAwk,ncol=D,nrow=D)}
              else {th$Sigma[,,k]=th$Sigma[,,k]+diag(diag(ASAwk));}    
            }            
            else
            {#%%% Isotropic Sigma
              th$Sigma[,,k]=mean(sigma2)*diag(D); #% DxD
              if (is.null(dim(ASAwk))) {th$Sigma[,,k]=th$Sigma[,,k]+sum(diag(ASAwk,ncol=D,nrow=D))/D*diag(D);}
              else {th$Sigma[,,k]=th$Sigma[,,k]+sum(diag(ASAwk))/D*diag(D);}
            }  
            }                       
            else {	cstr$Sigma ;
              stop('  ERROR: invalid constraint on Sigma.');}
          }
        }
        
        #% Avoid numerical problems on covariances:
        if(verb>=3) print('n');
        if(! is.finite(sum(th$Gamma[1:Lt,1:Lt,k]))) {th$Gamma[1:Lt,1:Lt,k]=0;}
        th$Gamma[1:Lt,1:Lt,k]=th$Gamma[1:Lt,1:Lt,k]+1e-8*diag(Lt);
        if(! is.finite(sum(th$Sigma[,,k]))) {th$Sigma[,,k]=0;}
        th$Sigma[,,k]=th$Sigma[,,k]+1e-8*diag(D);
        if(verb>=3) print(',');
      } # end for on k
      
      if(verb>=3) print('end');
      
      if (cstr$Sigma=="*")
      {#%%% Equality constraint on Sigma
        th$Sigma=sweep(th$Sigma ,3,N,"*"); 
        th$Sigma=array(apply(th$Sigma,c(1,2),mean),dim=c(D,D,K)) 
      }
      
      if( !is.null(cstr$Gammat) && cstr$Gammat=='v')
      {#%%% Equal volume constraint on Gamma
        detG=rep(0,K);
        for (k in 1:K){
          if (D==1) {detG[k]=th$Gamma[1:Lt,1:Lt,k]}
          else {detG[k]=det(th$Gamma[1:Lt,1:Lt,k]);} #% 1x1
          th$Gamma[1:Lt,1:Lt,k] = th$Gamma[1:Lt,1:Lt,k] / detG[k]
        }
        th$Gamma[1:Lt,1:Lt,]=sum(detG^(1/Lt)*th$pi)*th$Gamma[1:Lt,1:Lt,];
      }
      
      if(is.character(cstr$Gammat) && !is.null(cstr$Gammat) && cstr$Gammat[length(cstr$Gammat)]=='*')
      {#%%% Equality constraint on Gammat
        for (k in 1:K){
          th$Gamma[1:Lt,1:Lt,k]=th$Gamma[1:Lt,1:Lt,k]%*%matrix(N,1,1);    
          th$Gamma[1:Lt,1:Lt,k]=matrix(1,Lt,Lt) * sum(th$Gamma[1:Lt,1:Lt,k])/N;  
        }  
      }
      
     # if( ! is.character(cstr$pi) || is.null(cstr$pi))
     # {if(! is.null(cstr$pi)) {th$pi=cstr$pi;}} else {
      #  if (!is.null(cstr$pi) && cstr$pi[1]=='*') 
       # {th$pi=1/K*rep(1,K);} else {stop('  ERROR: invalid constraint on pi.');} 
      #} 
      
      sumthSigma <- sumthSigma + th$Sigma
      sumthA <- sumthA + th$A
      sumthb <- sumthb + th$b
    } # end for on irep 
    
    th$A<-sumthA/iR
    th$b<-sumthb/iR
    th$Sigma<-sumthSigma/iR
    #print(th$A[,,1])
    return(th)
  }  # end function Maximization
  
  
  
  
  # % ==========================EM initialization==============================
  Lt=nrow(tapp)
  L=Lt+Lw;
  D = nrow(yapp)/iR ; N = ncol(yapp);
  
  if(Lw==0) {Sw=NULL; muw=NULL;} 
  
  # % =====================MAXIMIZATION STEP===========================
  theta = MaximizationIIDK1(tapp,yapp,yappR,iR,muw,Sw,cstr,verb);    
  
  theta$Gamma = theta$Gamma[,,1]
  theta$Sigma = theta$Sigma[,,1]
  theta$A = theta$A[,,1]
  
  # Compute Sigma* and all useful to compute the posteriors later
  invGamma<-chol2inv(chol(theta$Gamma))
  invSigma<-chol2inv(chol(theta$Sigma))
  tempMat<-iR*t(theta$A)%*%invSigma%*%(theta$A)
  covstar<-chol2inv(chol(invGamma+ tempMat))
  
  theta$invGamma<-invGamma
  theta$invSigma<-invSigma
  theta$covstarR = covstar

  
  # compute determinant: not necessary when K=1
 # dettemp<-det(diag(L)+tempMat%*%(theta$Gamma))
  #logdetVa<- iR*log(det(theta$Sigma))+log(dettemp)
 # theta$logdetVR = logdetVa
  
  #theta$r = r
  #theta$LLf=LLf
  #theta$LL = LL[1:iter]
  return(theta)
}


