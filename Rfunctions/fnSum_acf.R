fnSum_acf = function(x, Ka){
  # To compute autocovariance eg for a MA(2)  time series
   # x is time series data
   # Ka is number of acf lags
  T = length(x)
  S = rep(0,Ka)
  for (i in 1:Ka){
    S[i] = sum(x[(1+i):T]*x[1:(T-i)])/T
  }
  return (S)
}