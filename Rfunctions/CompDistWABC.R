CompDistWABC<-function(iR,yobs,zdata){
  # compute W2 distance between target sample and simulated samples for 
  # WABC. All samples are of size iR (=100) and dim D(=2)
  # yobs is eg nlmobs is DiR x 1
  # zdata is eg nlmdata2 is DiR x M
  # requires package transport
  # ex : CompDistABC(100, nlmbos,nlmdata2)
  
  # samples need to b iR x D
  Msim<-dim(zdata)[2]
  sobs<-t(matrix(yobs[,1],ncol=iR))
  aobs<-wpp(sobs, rep(1/iR, iR))
  
  ### for all other z 
  listdist<-NULL
  # foreach is used to use the cores of the computer)
  listdist<- foreach (i=1:Msim, .combine=cbind) %dopar% {
    sz<-t(matrix(zdata[,i],ncol=iR))
    az<-wpp(sz, rep(1/iR, iR))
    W2dist <- wasserstein(aobs,az,p=2)
  
    c(W2dist) 
  }
}

## all samples are of size 1000
#N0<-1000
#
#ComputeDist<-function(truesample){
#  a00<-wpp(truesample, rep(1/N0, N0))
 # a0<-wpp(MHvalITDTdf1s001, rep(1/N0, N0))
#  a1<-wpp(rej3ITDTdf1s001105$unadj.values, rep(1/N0, N0))
#  a2<-wpp(rej3ITDTlogvdf1s001$unadj.values,rep(1/N0, N0) )
#  a3<-wpp(postvaluesred3ITDTdf1s001, rep(1/N0, N0))
 # a3L2<-wpp(L2postvaluesred3ITDTdf1s001, rep(1/N0, N0))
 # a3MMDs1<-wpp(MMDpostvaluesred3ITDTdf1s001 , rep(1/N0, N0))
 # a4<-wpp(samplemixgllim, rep(1/N0, N0))
  
#  wMH<-wasserstein(a00,a0,p=2)
#  wexp<-wasserstein(a00,a1,p=2)
#  wlogv<-wasserstein(a00,a2,p=2)
#  wdistW<-wasserstein(a00,a3,p=2)
#  wdistL2<-wasserstein(a00,a3L2,p=2)
#  wdistMMDs1<-wasserstein(a00,a3MMDs1,p=2)
 # wgllim<-wasserstein(a00,a4,p=2)
  