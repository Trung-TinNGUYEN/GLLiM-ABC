L2scal2normalL1<-function(mu1, mu2, sigma1,sigma2){
#########################################################################
# L=1 case
# L2scal2normalL1 computes the L2 scalar product between two univariate
# Gaussians 
#########################################################################
  
  dnorm(mu1, mu2, sqrt(sigma1+sigma2))
}