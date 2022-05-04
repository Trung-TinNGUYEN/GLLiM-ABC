# March 2022: script to compute MSE for estimators in the MA(2) case using SMC-ABC

### Required packages: TO BE UPDATED need winference and BEWARE of the change
### with mvtnorm and the dmvnorm function...

# methods and algorithms
library(xLLiM)
library(abc)
library(mcmc)
library(transport)
library(abctools)

# simulations and pdf and misc tools 
library(uniformly)
library(mvtnorm)
library(MixSim)
library(expm)
library(emulator)  # for quadratic forms computation (optional)
library(sfsmisc) # to normalize true posteriors

library("corpcor")  # for is.positive.definite in xLLiM and GLLiM 

# to use multiple cores
library(foreach)
library(doParallel)
library(parallel)
library(MASS)

# For plots
library(ggplot2)

# Auxiliary function to load all functions in a directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
#Then for example to load all ABC-GLLiM functions
sourceDir("../Rfunctions")

# multiple core for foreach use
numCores <- detectCores() #numCores = 8
registerDoParallel(numCores)

## reload data: for the learning set  but SMC works differently so not
## useful
#load("../data-et-al4examples/dataMA2.RData")

##### use same data ymatobsMA2150 and ymatobsMA2AC2 as rejection ABC: see below
##### 
# save the 100 targets for re-use later
#save(list = c("ymatobsMA2150"), file = "ymatobsMA2150K20.RData")
#save(list = c("ymatobsMA2AC2"), file = "ymatobsMA2AC2K20.RData")
# and same gllim model:
modgllimused<-modgllimMA2R5IIDFullK20bis


# not re-done: 
#####################################
# Comparison with 2 Autocovariances as summary statistics
#####################################
# requires function fnSum_acf.R

# Here we use Ka= 2 autocovariances, ie =2 different lags
# target 
ytargetACF2MA2150<-fnSum_acf(ytargetMA2150,2)

# Learning and ABC set, from the same simulations as before
ysimu2ACF2MA2150<-apply(ysimu2MA2150, 2, fnSum_acf, Ka=2)
# dim is 2 x 10^5

# Relaod the desired gllim model
load("modgllimMA2R5IIDFullK30.RData")
# redo with K=20
load("modgllimMA2R5IIDFullK20.RData")



#modgllimused<-modgllimMA2R5IIDFullK20
modgllimused<-modgllimMA2R5IIDFullK20bis

## For MW2 and L2: Gllim model used ################
modg<-modgllimused
# notational shortcut: gllim direct parameters
modeleg<-modg$mod
Aa<-modeleg$A
Sigmaa<-modeleg$Sigma
ca<-modeleg$c
ba<-modeleg$b
pia<-modeleg$pi
covstarR<-modg$covstarR
logdetVR<-modg$logdetVR
invGammaa<-modg$invGamma
invSigmaa<-modg$invSigma

D= dim(ba)[1]
L=dim(ca)[1]
K<-dim(ba)[2]

wthr<-0.0001
#################################


# simulation of more observations (100) from theta=(0.6,0.2)
# requires function MA2 which is in MA2.R file

#MA2MSE<-function(modg,R,ymatobs, ymatobsAC2, ysimuset, ysimusetAC2, thetasimuset,wt=10^(-4)){
# modg: gllim model as computed by 
# ysimuset D x N simulation set for ABC eg N=10^5
# ysimusetAC2
# thetasimuset L x N simulated params
# ymatobs: D x M matix of M observation vector generated from the same theta
# ymatobsAC2 autocov 2 x M 
# eg D=150 and M=100
# ex of use:
ymatobsMA2150<-matrix(0,150,100)
for (ni in 1:100){ymatobsMA2150[,ni]<- MA2(c(0.6,0.2),ny=150)}
ymatobsMA2AC2<-apply(ymatobsMA2150, 2, fnSum_acf, Ka=2)
# res100<-MA2MSE(modgllimMA2R5IIDFullK30, R=5,ymatobsMA2150, ymatobsMA2AC2, ysimu2MA2150,ysimu2ACF2MA2150, thetasimu2MA2)

M=dim(ymatobsMA2150)[2]

# save the 100 targets for re-use later
save(list = c("ymatobsMA2150"), file = "ymatobsMA2150K20.RData")
save(list = c("ymatobsMA2AC2"), file = "ymatobsMA2AC2K20.RData")

#########################################################
### set SMC
library(winference)
registerDoParallel(cores = detectCores())
#rm(list = ls())
setmytheme()
gllim_colors <- get_gllim_colors()

set.seed(11)

doRun <- FALSE
max_time <- 30*60
d <- 30  # D 
len<-150

target <- get_ma2(len)
#target$parameters$tau <- 5
nobservations <- 5  # R or 
nparticles <- 2048
p <- 1
prefix <- ""

# function to simulate data
target$simulate <- function(theta){
  return(target$robservation(nobservations, theta))
}

# common algorithmic parameters
param_algo <- list(nthetas = nparticles, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)


#####################################################


#########################################################
## MSE loop starts here 
#########################################################
#
# plantw with foreach, essaye avec for
#system.time(listestim<- foreach (i=1:M, .combine=cbind) %dopar% { 
# 

# run each method separately already 16 to 17 h
#listestimSMC<-array(0, c(3,5,M))

D= dim(ba)[1]
L=dim(ca)[1]
K<-dim(ba)[2]
wthr<-0.0001

listestimSMC<-array(0, c(3,5,M))
listestimSMC100<-array(0, c(3,5,M))

for(i in 1:M){
 # for(i in 95:100){
  ytarget<-t(ymatobsMA2150[,i]) # 1 x  150
  
  # SMC WABC and GLLiM
  # obs the observations D x R
  obs<-matrix(t(ytarget), ncol=5)
  
  # wasserstein distance to obs
 ##1 wdistance <- get_transport_to_y(obs, p = p)
  #
  
  # L2 and MW2 distance to obs: depend on obs !!!!
  # Compute post mixture weights for y (=obs): Piy dim is 1 x K 
  logPiy<-NULL
  tmp<-matrix(0,D,K)
  tmpM<-array(0,c(D,D,K))
  # t(t()) trick for L=1 case
  for(k in 1:K){
    tmpk<-Aa[,,k]%*%t(t(ca[,k]))+ba[,k]
    tmpyDR<-t(rowSums(obs-tmpk[,1]))
    # tmpM pourrait etre calculer avant but is reused via tmpM anyway
    tmpMk<-invSigmaa[,,k]%*%Aa[,,k]%*%t(t(covstarR[,,k]))%*%t(Aa[,,k])%*%invSigmaa[,,k]
    logpiyk<--0.5*sum(quad.diag(invSigmaa[,,k],(obs-tmpk[,1]))) +0.5*quad.form(tmpMk,t(tmpyDR))-0.5*logdetVR[k]
    #logpiyk<--0.5*sum((yDR-tmpk[,1])^2/Sigmaa[1,1,k])+0.5*quad.form(tmpMk,t(tmpyDR))-0.5*logdetVR[k]
    logPiy<-c(logPiy, log(pia[k])+logpiyk)
    tmp[,k]<-tmpk
    tmpM[,,k]<-tmpMk
  }
  # added for very small values? 
  #logPiy<-logPiy-min(logPiy)
  #postweightyz<-gllimpredyz$alpha
  
  den=mylogsumexp(logPiy);
  logPiy= logPiy-den 
  Piy=exp(logPiy); 
  
  # in case of very low weights
  seuil<-min(wthr,sort(Piy, decreasing=TRUE)[3])
  leftky<-seq(1,K)[Piy>seuil]
  dimay<-sum(Piy>seuil)
  Piy<-Piy[Piy>seuil]
  #Piy<-Piy/sum(Piy)
  ### end wieghts computation for obs
  
  postcovy<-array(0,c(L,L,dimay))
  postmeany<-matrix(0,L,dimay)
  postcovy[,,1:dimay]=covstarR[,,leftky]
  
  # Pre-computation of the trace part in the Wassertein distance, does not
  # depend on z or y and could even be pre-computed out of the function like
  # covstara (todo)
  # if computed here, it depends on y only via leftky 
  # tracecost is leftky x K 
  tracecost<-matrix(0,dimay,K)
  for (ii in 1:dimay) {
    sigma1<-postcovy[,,ii]
    for (jj in 1:K){
      sigma2<-covstarR[,,jj]
      E2 <- eigen(sigma2)
      V2 <- E2$vectors
      U2 <- solve(V2)
      D2 <- diag(E2$values) 
      sqrt2 <- V2 %*% D2^(1/2) %*% U2
      E <- eigen(sqrt2 %*%sigma1 %*% sqrt2)
      V <- E$vectors
      U <- solve(V)
      DD <- diag(E$values) 
      sqrtout <- V %*% DD^(1/2) %*% U
      tracecost[ii,jj] <-  sum(diag(sigma1+sigma2 - 2*sqrtout))
    }
  }
  #
  
  Asybs<-NULL
  for (k in leftky){
    # FULL case
    Asybs<-cbind(Asybs,t(t(covstarR[,,k]))%*%(t(Aa[,,k])%*%invSigmaa[,,k]%*%t(t(rowSums(obs-ba[,k]))) + invGammaa[,,k]%*%t(t(ca[,k])) ) )
  }
  # L x dimay  (y)
  postmeany=Asybs
  # problem with networkflow when K=1
  # to avoid we duplicate the mixture in this case
  if (dimay==1) {
    postmeany<-cbind(postmeany,postmeany)
    temp=array(0,c(L,L,2))
    temp[,,1]<-postcovy
    temp[,,2]<-postcovy
    postcovy<-temp
    Piy=c(0.5,0.5)
    dimay<-2
  } 
  
  #
  # Pre computation of the first quadratic form (dep on y) in the L2 distance
  # could be improved probably (todo)
##1  gramy <- matrix(NA,nrow=dimay,ncol=dimay) # symmetic
##1  for (ii in 1:dimay) {
##1    for (jj in ii:dimay) {
##1      gramy[ii,jj] <- L2scal2normalSMC(postmeany[,ii],
##1                                       postmeany[,jj],
##1                                       postcovy[,,ii],
##1                                       postcovy[,,jj]);
##1      gramy[jj,ii] <- gramy[ii,jj]
##1    }
##1  }
##1  qgramy<-quad.form(gramy,Piy)
  
  ####
  #### Function for computing the MW2 distance between obs and some z
  # MW2 for K>1
  dgllimMW2<-function(z){
    logPiz<-NULL
    for (k in 1:K){
      tmpzDR<-t(rowSums(z-tmp[,k]))
      logpizk<--0.5*sum(quad.diag(invSigmaa[,,k],(z-tmp[,k]))) +0.5*quad.form(tmpM[,,k],t(tmpzDR))-0.5*logdetVR[k]
      logPiz<-c(logPiz, log(pia[k])+logpizk)
    }
    # logPiz 
    den=mylogsumexp(logPiz);
    logPiz= logPiz-den
    #Pimatz was M x K---Piz
    Piz=exp(logPiz)
    seuil<-min(wthr,sort(Piz, decreasing=TRUE)[3])
    leftkz<-seq(1,K)[Piz>seuil]
    dimaz<-sum(Piz>seuil)
    Piz<-Piz[Piz>seuil]
    
    postcovz<-array(0,c(L,L,dimaz))
    postmeanz<-matrix(0,L,dimaz)
    postcovz[,,1:dimaz]=covstarR[,,leftkz]
    
    Asybs<-NULL
    for (k in leftkz){
      # ISOTROPIC case
      Asybs<-cbind(Asybs,t(t(covstarR[,,k]))%*%(t(Aa[,,k])%*%invSigmaa[,,k]%*%t(t(rowSums(t(t(z))-ba[,k]))) + invGammaa[,,k]%*%t(t(ca[,k])) ) )
    }
    # L x dimaz  (y)
    postmeanz=Asybs
    
    if (dimaz==1) {
      postmeanz<-cbind(postmeanz,postmeanz)
      temp=array(0,c(L,L,2))
      temp[,,1]<-postcovz
      temp[,,2]<-postcovz
      postcovz<-temp
      Piz=c(0.5,0.5)
      dimaz<-2
      leftkz<-c(leftkz,leftkz)
    } 
    
    mixy<-list("Mu"=postmeany, "S"=postcovy, "Pi"=Piy)
    mixz<-list("Mu"=postmeanz, "S"=postcovz, "Pi"=Piz)
    
    # tracecost has to be reduced to leftky x leftkz before the call
    #tracecostz<-matrix(0,dimay,dimaz)
    tracecostz<-tracecost[,leftkz]
    
    #MW2dist <- was2mixPreComp(mixy,mixz,tracecostz)
    return(was2mixPreComp(mixy,mixz,tracecostz) )
  } # end function distance
  

  ####
  #### Function for computing the L2 distance between obs and some z
  # L2 for K>1
##1  dgllimL2<-function(z){
##1    logPiz<-NULL
##1    for (k in 1:K){
##1      tmpzDR<-t(rowSums(z-tmp[,k]))
##1      logpizk<--0.5*sum(quad.diag(invSigmaa[,,k],(z-tmp[,k]))) +0.5*quad.form(tmpM[,,k],t(tmpzDR))-0.5*logdetVR[k]
##1      logPiz<-c(logPiz, log(pia[k])+logpizk)
##1    }
##1    # logPiz 
##1    den=mylogsumexp(logPiz);
##1    logPiz= logPiz-den
##1    #Pimatz was M x K---Piz
##1    Piz=exp(logPiz)
##1    seuil<-min(wthr,sort(Piz, decreasing=TRUE)[3])
##1    leftkz<-seq(1,K)[Piz>seuil]
##1    dimaz<-sum(Piz>seuil)
##1    Piz<-Piz[Piz>seuil]
    
##1    postcovz<-array(0,c(L,L,dimaz))
##1    postmeanz<-matrix(0,L,dimaz)
##1    postcovz[,,1:dimaz]=covstarR[,,leftkz]
    
##1    Asybs<-NULL
##1    for (k in leftkz){
      # ISOTROPIC case
##1      Asybs<-cbind(Asybs,t(t(covstarR[,,k]))%*%(t(Aa[,,k])%*%invSigmaa[,,k]%*%t(t(rowSums(t(t(z))-ba[,k]))) + invGammaa[,,k]%*%t(t(ca[,k])) ) )
##1    }
    # L x dimaz  (y)
##1    postmeanz=Asybs
    
##1    if (dimaz==1) {
##1      postmeanz<-cbind(postmeanz,postmeanz)
##1      temp=array(0,c(L,L,2))
##1      temp[,,1]<-postcovz
##1      temp[,,2]<-postcovz
##1      postcovz<-temp
##1      Piz=c(0.5,0.5)
##1      dimaz<-2
##1      leftkz<-c(leftkz,leftkz)
##1    } 
    
##1    mixy<-list("Mu"=postmeany, "S"=postcovy, "Pi"=Piy)
##1    mixz<-list("Mu"=postmeanz, "S"=postcovz, "Pi"=Piz)
    
    
    #L2dist <- L2normalPreComp(mixy,mixz,qgramy)
##1    return(L2normalPreCompSMC(mixy,mixz,qgramy))
##1  } # end function distance
  
  ###### end prep MW2 and L2 
  
  
  ###### SMC WABC, MW2 and L2
##1  resultsma2.wabc <- wsmc(wdistance, target, param_algo, maxsimulation = 10^5)
  resultsma2.MW2 <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^5)
##1  resultsma2.L2 <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^5)
  
  # list of 2048 selected thetas, reduced to the best 100 distances
##1  step<-length(resultsma2.wabc$thetas_history)
  # 2048 distances
##1  dist<-resultsma2.wabc$distances_history[[step]]
##1  ordist<-order(dist)
  # 2048 x 2
##1  thetas<-resultsma2.wabc$thetas_history[[step]]
##1  thetas100.wabc<-thetas[ordist[1:100],]
  #
  step<-length(resultsma2.MW2$thetas_history)
  # 2048 distances
  dist<-resultsma2.MW2$distances_history[[step]]
  ordist<-order(dist)
  # 2048 x 2
  thetas.MW2<-resultsma2.MW2$thetas_history[[step]]
  thetas100.MW2<-thetas.MW2[ordist[1:100],]
  #
##1  step<-length(resultsma2.L2$thetas_history)
  # 2048 distances
##1  dist<-resultsma2.L2$distances_history[[step]]
##1  ordist<-order(dist)
  # 2048 x 2
##1  thetas<-resultsma2.L2$thetas_history[[step]]
##1  thetas100.L2<-thetas[ordist[1:100],]
  #####
  
  
  # computing various estimators : thetas are 100 x2
##1  estimSMC<-matrix(0,3,5) # methods x estim 2 means, 2 std, 1 cor
  estimSMC100<-matrix(0,3,5) # methods x estim 2 means, 2 std, 1 cor
  estimSMC<-matrix(0,3,5)
  #methods order: SA150, SAAC2, AC2,GM, E, EV, L2, MW2
  # means
##1 estimSMC[1,1:2]<-colMeans(thetas100.wabc)
  estimSMC100[2,1:2]<-colMeans(thetas100.MW2)
  estimSMC[2,1:2]<-colMeans(thetas.MW2)
 ##1 estimSMC[3,1:2]<-colMeans(thetas100.L2)
  # std
 ##1 estimSMC[1,3:4]<-sqrt(diag(cov(thetas100.wabc) ) )
  estimSMC100[2,3:4]<-sqrt(diag(cov(thetas100.MW2) )) 
  estimSMC[2,3:4]<-sqrt(diag(cov(thetas.MW2) )) 
 ##1 estimSMC[3,3:4]<-sqrt(diag(cov(thetas100.L2) )) 
  # cor cor(t(rejsampleMA2R5$L2postval)) but cor(rejsampleEMA2R5$unadj.values)
 ##1 estimSMC[1,5]<-cor(thetas100.wabc)[1,2]
  estimSMC100[2,5]<- cor(thetas100.MW2)[1,2]
  estimSMC[2,5]<- cor(thetas.MW2)[1,2]
 ##`` estimSMC[3,5]<-cor(thetas100.L2)[1,2] 
  
  listestimSMC100[,,i]<-estimSMC100
  listestimSMC[,,i]<-estimSMC
} # end for
#) # end system.time
# listestim
#}  
#
#  pour M=1 environ 10 min pour chaque methode! MW2: 694 s 

save(list = c("listestimSMC100"), file = "listestimSMC100MA2K20bis.RData")
save(list = c("listestimSMC"), file = "listestimSMC2048MA2K20bis.RData")
#

###### use listtrue values already computed for other methods
##### repeat  with Gllim K30 for IS that for K20bis
# Compute numerically true values for the 100 observations ymatobsMA2150
#save(list = c("listtrue"), file = "listtrue100MA2150K20bis.RData")
#



truematSMC<-array(0,c(3,5,M))
for(i in 1:M){
  for(j in 1:3){
    truematSMC[j,,i]<-listtrue[,i]
  }
}

# compute MSE with 100 best for SMC
msematSMC100<-(listestimSMC100 -truematSMC)^2
msemSMC100<-matrix(0,3,5)
for(i in 1:3){
  msemSMC100[i,]<-rowMeans(msematSMC100[i,,])
}

# compute MSE 2048
msematSMC<-(listestimSMC -truematSMC)^2
msemSMC<-matrix(0,3,5)
for(i in 1:3){
  msemSMC[i,]<-rowMeans(msematSMC[i,,])
}

#save(list = c("msem"), file = "MSE100MA2K20.RData")
# REM: only the second line in the matrix is ok, others have not been computed
# MW2-SMC
save(list = c("msemSMC"), file = "MSESMC2048MA2K20bis.RData")
save(list = c("msemSMC100"), file = "MSESMC100MA2K20bis.RData")


