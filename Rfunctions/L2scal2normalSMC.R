L2scal2normalSMC<-function(mu1, mu2, sigma1,sigma2){
  ##### L2scal2normal computes the L2 scalar product between two gaussians
  ##### needed for L2 computation
  # mu1 is a row 1 x dim and mu2 is a vector dim x 1 or a row
  dmvnorm(t(mu1), t(mu2), sigma1+sigma2)
}