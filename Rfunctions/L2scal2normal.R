L2scal2normal<-function(mu1, mu2, sigma1,sigma2){
##### L2scal2normal computes the L2 scalar product between two gaussians
##### needed for L2 computation
  # mu1 is a row 1 x dim and mu2 is a vector dim x 1 or a row
  # March 2022: add of transpose due to a change in dmvnorm from mclust
  dmvnorm(t(mu1), t(mu2), sigma1+sigma2)
}