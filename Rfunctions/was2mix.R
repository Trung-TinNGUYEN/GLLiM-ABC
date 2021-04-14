was2mix<-function(mix1,mix2){
##################################################################################
## Squared Mixture Wasserstein distance (MW2 in Delon and Desolneux 2020) between two
## Gaussian mixtures
##################################################################################
## ATTENTION: there is a problem when the mixture has only one component because 
# the mean needs to be 1 x L matrix and not a vector
# this is solved by a trick in ComputeDist functions that consists in writting 
# a one component mixture into a two component mixture with equal components
##################################################################################
#### Input:
# - mix1 or mix2 must contain: Mu mean of L x K, S covariance matrix of L x L x K
# and Pi weight vector of size K
#### Ouput:
# - squared MW2 distance ie MW2^2
##################################################################################
  
  # Compute the cost matrix
  cost. <- matrix(NA,nrow=length(mix1$Pi),ncol=length(mix2$Pi))
  for (ii in 1:length(mix1$Pi)) {
    for (jj in 1:length(mix2$Pi)) {
      cost.[ii,jj] <- was2normal(mix1$Mu[,ii],
                                 mix2$Mu[,jj],
                                 mix1$S[,,ii],
                                 mix2$S[,,jj])
    }
  }
  
  # Compute the optimal transportation cost
  plan. <- transport(mix1$Pi,mix2$Pi,costm = cost.,
                     fullreturn = TRUE,
                     #method = "revsimplex")
                     method = 'networkflow')
  plan.$cost
}
