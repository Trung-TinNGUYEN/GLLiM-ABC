GllimISMA2<-function(thetaMat, pival, muval,sval,yobs){
  # compute means, covariance matrix, correlation matrix  of the true posterior of a MA(2)
  # model when the observation is yobs
  # the computation is made using the unnormalized unormp expression of the 
  # posterior and by importance sampling using a sample from a K-Gaussian
  # mixture proposal distribution, on the same space as theta
  #
  # INPUT: 
  # the K-Gaussian mixture used as IS proposal is defined by its parameters: 
  # pival (proportions): vector of length K
  # muval means: matrix dim of theta x K
  # sval covariances: array dim theta x dim theta x K
  # thetaMat: matrix Mdim x dim of theta, theta values used for the monte-carlo sum
  #
  # OUTPUT
  # ISpds: normalised importance weights
  # ISmean, ISvar, IScor: estimations for marginals means, variances and correlations
  # in the MA2 case
  #
     Kdim<-length(pival)
      Mdim<-dim(thetaMat)[1]
      
      # compute unnormalized IS weights in listpds
      listpds<- foreach (i=1:Mdim, .combine=c) %dopar% {
        # LikeMA2 is the Likelihood  x the prior uniform on the triangle
        unormp<-LikeMA2(thetaMat[i,],yobs)
        if (unormp==0)
          pds<-0 else{
            dmixt<-0
            for (k in 1:Kdim){
              dmixt<-dmixt + pival[k]*dmvnorm(thetaMat[i,],muval[,k], sval[,,k])
            }
            pds<-unormp/dmixt
          }
        pds
      }
      # normalized importance weights
      listpds<-listpds/sum(listpds)
      # compute IS means, std et correlation
      # method= "unbiased" (default) or "ML" 
      compest<-cov.wt(thetaMat,listpds, cor=TRUE, center=TRUE)
      
      # other ways to compute without cov.wt but not tested:
      #ISmean<-listpds%*%thetaMat
      #thetaMattemp<-thetaMat
      #thetaMattemp[,1]<-thetaMat[,1]-ISmean[1]
      #thetaMattemp[,2]<-thetaMat[,2]-ISmean[2]
      #ISstd<-sqrt(listpds%*%(thetaMattemp)^2)
      #thetaMattemp<-thetaMattemp[,1]*thetaMattemp[,2]
      #IScor<-(listpds%*%thetaMattemp)/(ISstd[1]*ISstd[2])
      
      list("ISpds"=listpds, "ISmean"= compest$center, "ISvar"=compest$cov, "IScor"=compest$cor)
}
