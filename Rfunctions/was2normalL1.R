was2normalL1 <- function(mu1,mu2,sigma1,sigma2) {
##############################################################################
## Case L = 1   FASTER  
## Compute the square Wasserstein distance for a L2 cost, between 2 univariate 
## Gaussian distributions.
##############################################################################
##### Input :
## - mu_i : means of the Gaussians
## - sigma1: variances of the Gaussians
##### Output:
## - sq : squared Wasserstein distance
##############################################################################

  sq. <- (mu1-mu2)^2 + sigma1+sigma2 - 2*sqrt(sigma1*sigma2)
  return(max(sq.,0))
}
