L2scal2normal<-function(mu1, mu2, sigma1,sigma2){
##### L2scal2normal computes the L2 scalar product between two gaussians
##### needed for L2 computation
  dmvnorm(mu1, mu2, sigma1+sigma2)
}