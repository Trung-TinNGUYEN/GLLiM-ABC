L2normal<-function(mix1,mix2){
##################################################################################
##### L2normal computes the squared L2 distance between two Gaussian mixtures
##################################################################################
### Input:
# - mix1 or mix2 must contain: Mu mean of L x K, S covariance matrix of L x L x K
# and Pi weight vector of size K
#### Ouput:
# - squared L2 distance 
##################################################################################
  
  # Compute the Gram matrices (scalar products)
  K1=length(mix1$Pi)
  K2=length(mix2$Pi)
  gram1 <- matrix(NA,nrow=K1,ncol=K1) # symmetic
  gram2 <- matrix(NA,nrow=K2,ncol=K2) # symmetric
  gram12 <- matrix(NA,nrow=K1,ncol=K2) # not symmetric
  
  for (ii in 1:K1) {
    for (jj in ii:K1) {
      gram1[ii,jj] <- L2scal2normal(mix1$Mu[,ii],
                                    mix1$Mu[,jj],
                                    mix1$S[,,ii],
                                    mix1$S[,,jj]);
      gram1[jj,ii] <- gram1[ii,jj]
    }
  }
  for (ii in 1:K2) {
    for (jj in ii:K2) {
      gram2[ii,jj] <- L2scal2normal(mix2$Mu[,ii],
                                    mix2$Mu[,jj],
                                    mix2$S[,,ii],
                                    mix2$S[,,jj]);
      gram2[jj,ii] <- gram2[ii,jj]
    }
  }
  for (ii in 1:K1) {
    for (jj in 1:K2) {
      gram12[ii,jj] <- L2scal2normal(mix1$Mu[,ii],
                                     mix2$Mu[,jj],
                                     mix1$S[,,ii],
                                     mix2$S[,,jj])
    }
  }
  
  # Compute the L2 between the 2 mixtures
  # if not using emulator package replace quad.form by
  # t(x)%*%M%*%x etc...
  L2dist <- quad.form(gram1,mix1$Pi) + quad.form(gram2,mix2$Pi) - 2 * quad.3form(gram12,mix1$Pi,mix2$Pi) 
  
  L2dist
}

