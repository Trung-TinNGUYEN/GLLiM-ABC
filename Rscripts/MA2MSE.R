####################################################################
##
## July 2021: script to compute MSE over 100 repetitions of the same experiment
## for estimators in the MA(2) case 
##
####################################################################

### Required packages

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

## reload data
load("../data-et-al4examples/dataMA2.RData")


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
  
  #############################################################
  #
  # Do not run with foreach, so switch to a for loop
  #####system.time(listestim<- foreach (i=1:M, .combine=cbind) %dopar% { 
  
  listestim<-array(0, c(8,5,M))
  #### For loop
  ####
  for(i in 1:M){
  ytarget<-t(ymatobsMA2150[,i])
  # L2 and MW2 distance computation
  distMW2L2<-PreCompDistFullParIID(modgllimMA2R5IIDFullK30,R=5,t(ytarget),ysimu2MA2150,wthr = 0.0001)
  
  # for 100 selected points 
  rejsampledist<-RejABCMW2L2(0.001,distMW2L2[1,],distMW2L2[2,],thetasimu2MA2)
  
  #### comppute gllim posterior expectation and log variances in the IID case
  ExpVarMatz<-ExpLogPostVarGllimFullIID(modgllimMA2R5IIDFullK30,ysimu2MA2150,R=5)
  # for the current target 
  ExpVarTarget<-ExpLogPostVarGllimFullIID(modgllimMA2R5IIDFullK30,t(ytarget),R=5)
  #L=2 for MA2
  espysimu<-ExpVarMatz[,1:2]
  # espytarget is 1 x L 
  espytarget<-ExpVarTarget[1:2]
  
  #ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
  # rejection is the one used for wasserstein ...
  # tol =0.01 select 1000 samples when Nsimu =10^5]
  # thetasimu2MA2 is L x M
  # selected L=1sample of size L x 100 is in rejsampleESMA1$unadj.values
  rejsampleE <-abc(target=espytarget, param=t(thetasimu2MA2), sumstat=espysimu, tol=.001, method ="rejection")
  ### exp + log var 
  rejsampleEV <-abc(target=ExpVarTarget, param=t(thetasimu2MA2), sumstat=ExpVarMatz, tol=.001, method ="rejection")
  
  mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
  #mytf<-list(function(x){cbind(x,x^2,x^3,x^4,x^5,x^6,x^7,x^8 )})
  # with x to x^8  is the same ?
  # DONT Forget to change: for MA2 threshold is 0.001
  saabc.MA2<-semiauto.abc(obs=ytarget,param=t(thetasimu2MA2), sumstats=t(ysimu2MA2150),satr=mytf, tol=0.001, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)
  
  # gllim mixture for yobs
  ## 2D plot GLLiM mixture for target 100 points
  mixp<-ComputeGllimFullMixtParamIID(modgllimMA2R5IIDFullK30,R=5,t(ytarget),0)
  samplemixg<-simdataset(100,mixp$Pi, t(mixp$Mu),mixp$S)$X
  
  #####################################
  # Comparison with 2 Autocovariances as summary statistics
  #####################################
 # Here we use Ka= 2 autocovariances, ie =2 different lags
  # target 
  ytargetAC2<-t(ymatobsMA2AC2[,i] ) 
 
  # ABC et SA 
  # ABC rejection
  rejsampleAC2 <-abc(target=ytargetAC2, param=t(thetasimu2MA2), sumstat=t(ysimu2ACF2MA2150), tol=.001, method ="rejection")
  
  #### SA on data of size D=2 AC 
  ## SEMI AUTO ABC-- 
  saabc.AC2<-semiauto.abc(obs=ytargetAC2,param=t(thetasimu2MA2), sumstats=t(ysimu2ACF2MA2150),satr=mytf, tol=0.001, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)
  
  # computing various estimators
  estim<-matrix(0,8,5) # methods x estim 2 means, 2 std, 1 cor
  #methods order: SA150, SAAC2, AC2,GM, E, EV, L2, MW2
  # means
  estim[1,1:2]<-colMeans(saabc.MA2$post.sample)
  estim[2,1:2]<-colMeans(saabc.AC2$post.sample)
  estim[3,1:2]<-colMeans(rejsampleAC2$unadj.values)
  estim[4,1:2]<-colMeans(samplemixg)
  estim[5,1:2]<-colMeans(rejsampleE$unadj.values)
  estim[6,1:2]<-colMeans(rejsampleEV$unadj.values)
  estim[7,1:2]<-rowMeans(rejsampledist$L2postval)
  estim[8,1:2]<-rowMeans(rejsampledist$MW2postval)
  # std
  estim[1,3:4]<-diag(sqrt(cov(saabc.MA2$post.sample) ) )
  estim[2,3:4]<-diag(sqrt(cov(saabc.AC2$post.sample) )) 
  estim[3,3:4]<-diag(sqrt(cov(rejsampleAC2$unadj.values) )) 
  estim[4,3:4]<-diag(sqrt(cov(samplemixg) ))
  estim[5,3:4]<-diag(sqrt(cov(rejsampleE$unadj.values)))
  estim[6,3:4]<-diag(sqrt(cov(rejsampleEV$unadj.values)))
  estim[7,3:4]<-diag(sqrt(cov(t(rejsampledist$L2postval)) ))
  estim[8,3:4]<-diag(sqrt(cov(t(rejsampledist$MW2postval)) ))
  # cor cor(t(rejsampleMA2R5$L2postval)) but cor(rejsampleEMA2R5$unadj.values)
  estim[1,5]<-cor(saabc.MA2$post.sample)[1,2]
  estim[2,5]<- cor(saabc.AC2$post.sample)[1,2]
  estim[3,5]<-cor(rejsampleAC2$unadj.values)[1,2] 
  estim[4,5]<-cor(samplemixg)[1,2]
  estim[5,5]<-cor(rejsampleE$unadj.values)[1,2]
  estim[6,5]<-cor(rejsampleEV$unadj.values)[1,2]
  estim[7,5]<-cor(t(rejsampledist$L2postval))[1,2]
  estim[8,5]<-cor(t(rejsampledist$MW2postval))[1,2]
  
  listestim[,,i]<-estim
  } ##### end for loop
   
  #) # end system.time
 # listestim
#}  
  #     user   system  elapsed 
  # 2053.367   45.293  680.174 --> 11min pour M=2
  # For M=100
   save(list = c("listestim"), file = "listestim100MA2K30.RData")
  #
   ########################################################################
   ## 
   ## Compute numerically true values for the 100 observations ymatobsMA2150
   ##
   listtrue<-matrix(0,5,M)
   for(i in 1:M){
     # for(i in 35:M){
     ytarget<-t(ymatobsMA2150[,i])
   ## sample from GLLiM mixture for target: 10000 or 10^5 points
   mixp<-ComputeGllimFullMixtParamIID(modgllimMA2R5IIDFullK30,R=5,t(ytarget),0)
   samplemixg<-simdataset(100000,mixp$Pi, t(mixp$Mu),mixp$S)$X
   compes<-GllimISMA2(samplemixg,mixp$Pi, mixp$Mu,mixp$S,ytarget)
   
   #TRUE values:
   listtrue[1:2,i]<-compes$ISmean
   #[1] 0.6348713 0.2027734 # close to numerical integration ones
   listtrue[3:4,i]<-diag(sqrt(compes$ISvar))
   #0.08024054 0.07577829 close to numerical integration ones
   listtrue[5,i]<-compes$IScor[1,2]
   }
   #####
   save(list = c("listtrue"), file = "listtrue100MA2150.RData")
   #####
   
   ################################################################
   # Check that the 4 first values are the same with numerical integration
   # because IS is made with the gllim mixture...
   # compute TRUE posterior means  by numerical integration
   ###############################################################
   # check on weird values (M=100) are the same with or without IS: 
   # not run for all observations then
   # listtrue[,100]
   #[1] 0.49646654 0.04714440 0.08862423 0.08902668 0.32607564
   #> listtruenum[,100]
   #[1] 0.49629898 0.04729529 0.08872330 0.08896931
   
   ###############################################################
   # True post for a MA2 model D=150: very concentrated 
   thetaseq1<-seq(-2,2,length=400)
   thetaseq2<-seq(-1,1,length=200)
   
   listtruenum<-matrix(0,4,M)
   
   for(i in 1:M){
     ytarget<-t(ymatobsMA2150[,i])
   postpdf<-PostPdfMA2(thetaseq1,thetaseq2, ytarget)
   #MA2_df<-postpdfMA2150$postdf
   
   mpdf1<-postpdf$marg1
   mpdf2<-postpdf$marg2
   meant1<-integrate.xy(thetaseq1, thetaseq1*mpdf1)
   meant2<-integrate.xy(thetaseq2, thetaseq2*mpdf2)
   # 0.6348032 and  meantheta2= 0.2027278
   stdt1<-sqrt(integrate.xy(thetaseq1, (thetaseq1-meant1)^2*mpdf1))
   stdt2<-sqrt(integrate.xy(thetaseq2, (thetaseq2-meant2)^2*mpdf2))
   # stdtheta1= 0.0801874 and stdtheta2= 0.07540907
  
   listtruenum[1,i]<-meant1
   listtruenum[2,i]<-meant2
   listtruenum[3,i]<-stdt1
   listtruenum[4,i]<-stdt2
   }
   
   ################################################################
   # compute MSE
   truemat<-array(0,c(8,5,M))
   for(i in 1:M){
   for(j in 1:8){
     truemat[j,,i]<-listtrue[,i]
   }
   }
   msemat<-(listestim -truemat)^2
   msem<-matrix(0,8,5)
   for(i in 1:8){
     msem[i,]<-rowMeans(msemat[i,,])
   }
   
  