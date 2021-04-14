
##############################################################################################
####  Example Sum of 2 MA(1) with opposite theta, D=10, L=1
##############################################################################################
#### Functions specific to this example, to create data sets, computation of the 
#### true posterior, all saved in a dataSMA1.Rdata file in the end


#################################################################################################
###
### Step 1: Setting: target y (observations to be inverted, dimension is D) ;
###         Data sets simulation, one for GLLiM, one for ABC if wanted
###
#################################################################################################
# input:
# output: ytarget; (thetasimu1, ysimu1) size N; (thetasimu2, ysimu2) size M ;
#  true posterior if possible

# To simulate from a MA(1), ny is the length (dimension) of the simulation
# default is ny =10 (D), requires package mvtnorm
# NOT REALLY USED here
MA1<-function(theta,ny=10){
  firstrow<-c(theta^2+1, theta, rep(0,ny-2))
  Sigmat<-toeplitz(firstrow)
  ysim<-rmvnorm(1,rep(0,ny),Sigmat)
  ysim
}

# To simulate from a sum of MA(1) with opposite theta
# Remark: This sum of MA1 is actually a diagonal Gaussian!! 
SMA1<-function(theta,ny=10){
  firstrow1<-c(theta^2+1, theta, rep(0,ny-2))
  Sigmat1<-toeplitz(firstrow1)
  firstrow2<-c(theta^2+1, -theta, rep(0,ny-2))
  Sigmat2<-toeplitz(firstrow2)
  Sigmat<-Sigmat1+Sigmat2
  ysim<-rmvnorm(1,rep(0,ny),Sigmat)
  ysim
}

########################################################
# Numerical computation of the true posterior for a SMA1 model
# when the observation is y.
# thetaseq are the theta values at which the pdf is computed and that
# depends on the chosen prior, eg thetaseq<-seq(-2,2,length=200)
PostPdfSMA1<-function(thetaseq, y){
# numerical normalisation for plot of the true posterior
likevalseq=NULL
thetaN<-length(thetaseq)
ny=length(y)

for(i in 1:thetaN){
  
  # likevalseq: unnormalized likelihood (prior Unif(-2,2) x Likelihood)
  theta<-thetaseq[i]
  firstrow1<-c(theta^2+1, theta, rep(0,ny-2))
  Sigmat1<-toeplitz(firstrow1)
  firstrow2<-c(theta^2+1, -theta, rep(0,ny-2))
  Sigmat2<-toeplitz(firstrow2)
  Sigmat<-Sigmat1+Sigmat2
  # to account for the prior here unif(-2,2)
  if ((theta< -2) | (theta>2) )
    # attention delta = mode not mean
    likeval<-0 else likeval<-dmvnorm(y,rep(0,ny),Sigmat)
  
  likevalseq<-c(likevalseq, likeval)
}
normc<-integrate.xy(thetaseq, likevalseq)
return(likevalseq/normc)
}
### Eg Application for the plot : plot(thetaseq,PostPdfSMA1(thetaseq,ytarget))

##################################################################################
## target observation simulated: D=10, thetatrue=1
ytargetSMA1<-SMA1(1,10)

##################################################################################
# training set size N, for GLLIM, to learn a parametric representation of the posteriors
# simulate N theta  uniformly in  [-2,2]
N =10^5
thetasimu1SMA1<-matrix(runif(N,-2,2), 1, N)
# dim L=1 x 10^5
# simulate as many corresponding  observations in dimension 10 (ny)
ysimu1SMA1<-apply(thetasimu1SMA1,2,SMA1)

# a second training set, size M for ABC 
M=10^5
thetasimu2SMA1<-matrix(runif(M,-2,2), 1, M)
# dim L=1 x 10^5
# simulate as many corresponding  observations in dimension 10 (ny)
ysimu2SMA1<-apply(thetasimu2SMA1,2,SMA1)
##################################################################################
# True posterior 
thetaseqSMA1<-seq(-2,2,length=200)
postpdfSMA1<-PostPdfSMA1(thetaseqSMA1,ytargetSMA1)


##################################################################################
# To save the data for later use: create a file dataSMA1.Rdata
save(list = c("ytargetSMA1","ysimu1SMA1","thetasimu1SMA1","ysimu2SMA1","thetasimu2SMA1", "thetaseqSMA1","postpdfSMA1" ), file = "dataSMA1.RData")
##################################################################################
 



