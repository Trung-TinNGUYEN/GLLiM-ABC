logsumexp = function(x){
  ######################################################################
  # A simplified version of the one in xLLiM
  # % Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
  # %   By default dim = 2 (columns).
  # if (is.null(dimension)) dimension=1 ;
  #  % subtract the largest in each column
  xmax = max(x)
  x = x-xmax #sweep(x,1,y,"-")
  s = xmax+log(sum(exp(x)))
  i = is.infinite(xmax)
  if (i > 0) s=xmax
  return(s)
}