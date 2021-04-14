was2normal <- function(mu1,mu2,sigma1,sigma2) {
##############################################################################
## Compute the square Wasserstein distance for a L2 cost, between 2 Gaussian
## distributions.
## Rem: with a faster way to compute square roots for sdp matrices
##############################################################################
##### Input :
## - mu_i : mean vectors of the Gaussians
## - sigma1: covariance matrices of the Gaussians
##### Output:
## - sq : squared Wasserstein distance
##############################################################################
 
  E2 <- eigen(sigma2)
  V2 <- E2$vectors
  U2 <- solve(V2)
  D2 <- diag(E2$values) 
  sqrt2 <- V2 %*% D2^(1/2) %*% U2
  E <- eigen(sqrt2 %*%sigma1 %*% sqrt2)
  V <- E$vectors
  U <- solve(V)
  D <- diag(E$values) 
  sqrtout <- V %*% D^(1/2) %*% U
  
  sq. <- sum((mu1-mu2)^2) + sum(diag(sigma1+sigma2 - 2*sqrtout))
  
  return(max(sq.,0))
}
