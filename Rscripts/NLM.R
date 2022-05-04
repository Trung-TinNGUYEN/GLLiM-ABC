##############################################################################################
####  Example Normal Location model, D=2, L=2, R=100 (from WABC Berton et al paper)
##############################################################################################
#### Functions specific to this example, to create data sets, computation of the 
#### true posterior, all saved in a dataNLM.Rdata file in the end
##############################################################################################

# iR=100  D=2 and L=2 nsimu=10^4 (=nb_iters)
# iR is n the number of iid 2-dim simulation for the same theta
# sigmalik<-matrix(c(1,0.5,0.5,1), ncol=2)

NLM_dataset_simu <- function(iR, nsimu, sigmalik) {
  simu_data <- matrix(NA, nsimu, 2*iR)  # here D=2
  stored_Theta <- matrix(NA, nsimu, 2) # here L=2
  
 # matrix nsimu x 2 of thetas from the prior
    rand_Theta <- rmvnorm(nsimu, mean=rep(0,2), sigma= diag(25,2))
    
    # Gllim training set nsimu x (iR x 2)
    for (k in 1:nsimu) {
    yk<-rmvnorm(iR,mean=rand_Theta[k,], sigmalik)
    simu_data[k,] <- as.vector(t(yk))
  }
  
  list(Data = simu_data, Theta = rand_Theta)
}

# True param (-0.71, 0.09) (cf WABC) and simulated observartions
truemean<-c(-0.71, 0.09)
sigmalik<-matrix(c(1,0.5,0.5,1), ncol=2)
iR=100
nlmobs<-as.vector(t(rmvnorm(iR,mean=truemean, sigmalik)))  # 1x 200

###################################### 
# simulation of 10^4 NLM learning set
Nsimu = 10^4  # (N)
learnset<-NLM_dataset_simu(iR=100,Nsimu, sigmalik)

nlmdata<-learnset$Data  # N x 2*iR (Nx 200)
nlmparams<-learnset$Theta  # N x L (N x2)
rm(learnset)
# transpose for xllim
nlmdata<-t(nlmdata) # 200 x 10^4
nlmparams<-t(nlmparams) # 2 x 10^4
nlmobs<-as.matrix(nlmobs)  # 200 x 1

###################################### 
# simulation of 10^4 NLM learning set for ABC
Msimu = 10^6  # (M)
learnset<-NLM_dataset_simu(iR=100,Msimu, sigmalik)

nlmdata2<-learnset$Data  # M x 2*iR (Mx 200)
nlmparams2<-learnset$Theta  # M x L (M x2)
rm(learnset)
# transpose for xllim
nlmdata2<-t(nlmdata2) # 200 x 10^6
nlmparams2<-t(nlmparams2) # 2 x 10^6


######################################
#True Posterior
# yobs is 2iR x 1 needs to be turned into 2x iR
# thetaseq is 2 x nbsample
PostPdfNLM<-function(thetaseq,yobs,sigmalik,iR){
  # post moments
  ysum<-rowSums(matrix(yobs, ncol=iR))
  invsigmalik<-solve(sigmalik)
    postsigma<-solve(diag(0.04,2)+iR*invsigmalik) 
      postmean<-postsigma%*%invsigmalik%*%ysum
      
      Ns<-dim(thetaseq)[2]
      ## Marginal posteriors
      margpdf1<-rep(0,Ns)
      margpdf2<-rep(0,Ns)
      jointpdf<-rep(0,Ns)
      for(i in 1:Ns){
        jointpdf[i]<-dmvnorm(thetaseq[,i],mean=postmean, sigma=postsigma)
        margpdf1[i]<-dnorm(thetaseq[1,i], mean=postmean[1],sd=sqrt(postsigma[1,1]) )
        margpdf2[i]<-dnorm(thetaseq[2,i], mean=postmean[2],sd=sqrt(postsigma[2,2]) )
      }
      list("marg1"=margpdf1,"marg2"=margpdf2,"postpdf"= jointpdf)
}

##################################################################################
# To save the data for later use: create a file dataBBM.Rdata
save(list = c("nlmobs","nlmdata","nlmparams","nlmdata2","nlmparams2"), file = "dataNLM.RData")
##################################################################################




###################################### 