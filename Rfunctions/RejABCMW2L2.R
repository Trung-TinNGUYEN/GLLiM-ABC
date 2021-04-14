RejABCMW2L2<-function(distq=0.01,distMW2,distL2,thetasimu2){
################################################################################
# To compute the selected samples with rejection ABC using MW2 and L2 for:
#  1) a posterior margin say i, using output of ComputeDistMargin for instance
#  or 2) the joint posterior (all L dimensions)
#  or 3) a subset of dimensions ...
### Input:
#    - distq: chosen quantile, default=0.01 (1% quantile)
#    - distMW2 : vector of MW2 distance values,length M
#    - distL2 : vector of L2 distance values, length M
#    - thetasimu2: theta simulated values, dim L x M (L=1, L=2  is possible)
#     !!REM!! if L=1 need to be a 1 x M matrix and not a vector
### Output: list with 2 matrices, each of dim L x S (S = distq * M)
#   - vector of selected theta sample for MW2
#   - vector of selected theta sample for L2
################################################################################
  
  cutoffMW2<-quantile(distMW2,distq)
  cutoffL2=quantile(distL2,distq)
  # Selected sample of size 100, 1000 or 500...depending on M and quantile choice
  MW2postval=thetasimu2[ ,distMW2<cutoffMW2]
  L2postval=thetasimu2[ ,distL2<cutoffL2]
  return(list("MW2postval"= MW2postval,"L2postval"= L2postval ))
}

