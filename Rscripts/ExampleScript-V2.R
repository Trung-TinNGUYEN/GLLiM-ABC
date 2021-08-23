# FF: From May to July 2021, Cleaned up in August 2021


#####################################################################################
##
## 2 Additional examples Bivariate Beta and MA(2) models
## using GLLiM-IID, added to the previous examples in ExampleScript.R 
## Read ExampleScript.R first, to see what libraries to upload and so on
##
######################################################################################


#####################################################################################
## 
## Moving Average MA(2) model:  
## target and simulations are of dimension  D=150
## Gllim learn with  Gllim-IID considering that the series are chunks of length 30, ie
## D=30 x R=5 iid observations
##
######################################################################################
#
# DATA: generated with MA2.R
# recover with load(dataMA2.Rdata): ysimu2MA2150 and ytargetMA2150 and thetasimu2MA2
# M=10^5 D =150

#
# True posterior for a MA2 model D=150: the posterior is very concentrated
# 
thetaseq1<-seq(-2,2,length=400)
thetaseq2<-seq(-1,1,length=200)
postpdfMA2150<-PostPdfMA2(thetaseq1,thetaseq2, ytargetMA2150)
MA2_df<-postpdfMA2150$postdf

### PLOT of contours of true posterior  for y=ytargetMA2150  (0.6, 0.2)
v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()


#######################################################################################
# Gllim IID with K =30 and Full covariances, ie the blocs have no specific structure (full sub matrices) 
# Priv-FF: Attention, GllimFitIID replaces GllimFullFitIID
#######################################################################################
### D=30 R=5
#############
system.time(modgllimMA2R5IIDFullK30 <- GllimFitIID(thetasimu2MA2, ysimu2MA2150,5, K, constr=list(Sigma="")))
#   user   system  elapsed 
# 921.249  146.975 1072.647 -->18 minutes
# Rem: Gllim mixture looks better than with R=2 or R=15

save(list = c("modgllimMA2R5IIDFullK30"), file = "modgllimMA2R5IIDFullK30.RData")
#

##############################
# Compute L2 and MW2 distances
##############################
system.time(distMW2L2MA2R5FullK30<-
              PreCompDistFullParIID(modgllimMA2R5IIDFullK30,R=5,t(ytargetMA2150),ysimu2MA2150,wthr = 0.0001)
)
#user    system   elapsed 
# 1038.259   36.907  184.252  thr=0.0001

save(list = c("distMW2L2MA2R5FullK30"), file = "distMW2L2MA2R5FullK30.RData")

##############################
# Select samples from distances 
# 0.001 for 100 points, 0.01 for 1000
rejsampleMA2R5<-RejABCMW2L2(0.001,distMW2L2MA2R5FullK30[1,],distMW2L2MA2R5FullK30[2,],thetasimu2MA2)

##############################
# compute gllim posterior expectation and log variances in the IID case
system.time(ExpVarMatzMA2R5FullK30<-ExpLogPostVarGllimFullIID(modgllimMA2R5IIDFullK30,ysimu2MA2150,5)
)
# user  system elapsed 
# 856.717  22.825 133.730 

# for the target (M=1)
# Remark: if target generated as IID, the t(t()) is needed, target has to be
# DR x 1 BUT HERE generated in one shot so should be t(ytargetMA2150)

ExpVarTargetMA2R5FullK30<-ExpLogPostVarGllimFullIID(modgllimMA2R5IIDFullK30,t(ytargetMA2150),5)

# espysimu2 is of size M x L 
#
L=2
espysimu2MA2R5IID<-ExpVarMatzMA2R5FullK30[,1:L]
# espytarget is 1 x L 
espytargetMA2R5IID<-ExpVarTargetMA2R5FullK30[1:L]

# ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for Wasserstein ...
# tol =0.001 select 100 samples when Nsimu =10^5]

####### GLLiM-E-ABC
# thetasimu2MA2 is L x M
# selected  L x 100 is in rejsampleEMA2R5$unadj.values
rejsampleEMA2R5 <-abc(target=espytargetMA2R5IID, param=t(thetasimu2MA2), sumstat=espysimu2MA2R5IID, tol=.001, method ="rejection")

#
####### GLLiM-EV-ABC: exp + log var 
#
rejsampleEVMA2R5 <-abc(target=ExpVarTargetMA2R5FullK30, param=t(thetasimu2MA2), sumstat=ExpVarMatzMA2R5FullK30, tol=.001, method ="rejection")


########################################
##
## SEMI AUTO ABC- SA on data of size 150 
## 
########################################
mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
#mytf<-list(function(x){cbind(x,x^2,x^3,x^4,x^5,x^6,x^7,x^8 )})
# with x to x^8  is the same ...
# DONT Forget to change: for MA2 threshold is 0.001
saabc.MA2150<-semiauto.abc(obs=ytargetMA2150,param=t(thetasimu2MA2), sumstats=t(ysimu2MA2150),satr=mytf, tol=0.001, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)



########################################
##
## PLOTS
##
########################################

dfL2MA2<-data.frame(t(rejsampleMA2R5$L2postval))
dfMW2MA2<-data.frame(t(rejsampleMA2R5$MW2postval))
dfEVMA2<-data.frame(rejsampleEVMA2R5$unadj.values)
dfEMA2 <-data.frame(rejsampleEMA2R5$unadj.values)
dfSAMA2150 <-data.frame(saabc.MA2150$post.sample)

########################################
### 2D Plots with superposed true post 
########################################

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=dfMW2MA2, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=dfL2MA2, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=dfEMA2, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=dfEVMA2, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=dfSAMA2150, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()

#########################
# gllim mixture for yobs
#########################
## 2D plot GLLiM mixture for target 100 points
mixparamPlot<-ComputeGllimFullMixtParamIID(modgllimMA2R5IIDFullK30,R=5,t(ytargetMA2150),0)
samplemixgllimPlot<-simdataset(100,mixparam$Pi, t(mixparam$Mu),mixparam$S)$X

dfgmixt <-data.frame(samplemixgllimPlot)

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=dfgmixt, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()


############################################################
##
## Comparison with 2 Autocovariances as summary statistics
##
############################################################
## requires function fnSum_acf.R

# Here we use Ka= 2 autocovariances, ie =2 different lags
# target 
ytargetACF2MA2150<-fnSum_acf(ytargetMA2150,2)

#############
# Learning and ABC set, from the same simulations as before
#############

ysimu2ACF2MA2150<-apply(ysimu2MA2150, 2, fnSum_acf, Ka=2)
# dim is 2 x 10^5

#############
# ABC and SA 
#############

# ABC rejection
rejsampleACF2MA2 <-abc(target=ytargetACF2MA2150, param=t(thetasimu2MA2), sumstat=t(ysimu2ACF2MA2150), tol=.001, method ="rejection")

#### SA on data of size D=2 AC 
## SEMI AUTO ABC-- 
mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
#mytf<-list(function(x){cbind(x,x^2,x^3,x^4,x^5,x^6,x^7,x^8 )})
# with x to x^8  is the same ?
# DONT Forget to change: for MA2 threshold is 0.001
saabc.ACF2MA2150<-semiauto.abc(obs=ytargetACF2MA2150,param=t(thetasimu2MA2), sumstats=t(ysimu2ACF2MA2150),satr=mytf, tol=0.001, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)

dfACF2MA2 <-data.frame(rejsampleACF2MA2$unadj.values)
dfSAACF2MA2150 <-data.frame(saabc.ACF2MA2150$post.sample)

#############
#Plots
#############

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=dfACF2MA2, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=dfSAACF2MA2150, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()


###########################################################
##
## EMPIRICAL Estimators for samples (ABC, gllim mixture...)
##
###########################################################

# parameters empirical means from samples, all equal weights...
rowMeans(rejsampleMA2R5$L2postval)
#[1] 0.7403013 0.2128260
 rowMeans(rejsampleMA2R5$MW2postval)
#[1] 0.7415236 0.2054516
 colMeans(rejsampleEVMA2R5$unadj.values)
 # 0.6892256 0.1625062
 colMeans(rejsampleEMA2R5$unadj.values)
 # 0.7371831 0.2079816
 colMeans(samplemixgllimPlot) 
 # 0.5205398 0.1817777
 colMeans(saabc.MA2150$post.sample)
# 0.02562663 0.02688614
 colMeans(rejsampleACF2MA2$unadj.values)
 # 0.6352177 0.2458847
 colMeans(saabc.ACF2MA2150$post.sample)
 # 0.6368611 0.2502731
 
 # to compute the std on the diagonal of
 sqrt(cov(t(rejsampleMA2R5$MW2postval))) 
 #
 #[,1]       [,2]
 #[1,] 0.10435316 0.06801877
 #[2,] 0.06801877 0.08413486
 
 sqrt(cov(t(rejsampleMA2R5$L2postval)))
 #[,1]       [,2]
 #[1,] 0.1028442 0.06292410
 #[2,] 0.0629241 0.08827913
 
 diag(sqrt(cov(rejsampleEMA2R5$unadj.values)))
 #0.10350075 0.08444854
 diag(sqrt(cov(rejsampleEVMA2R5$unadj.values)))
# 0.1099820 0.1204021
 diag(sqrt(cov(samplemixgllimPlot)))
 #0.4987139 0.2938641
 diag(sqrt(cov(saabc.MA2150$post.sample)))
 # 0.4061703 0.4757548
 diag(sqrt(cov(rejsampleACF2MA2$unadj.values)))
 #0.1069371 0.1644097
 diag(sqrt(cov(saabc.ACF2MA2150$post.sample)))
 # 0.1088209 0.1592887
 
 # for correlations
  cor(t(rejsampleMA2R5$L2postval))
 #[,1]      [,2]
 #[1,] 1.0000000 0.4361101
 #[2,] 0.4361101 1.0000000
  cor(t(rejsampleMA2R5$MW2postval))
  #[,1]     [,2]
  #[1,] 1.000000 0.526958
 # [2,] 0.526958 1.000000
# 
  cor(rejsampleEMA2R5$unadj.values)
#  1,] 1.0000000 0.4543587
# [2,] 0.4543587 1.0000000
  cor(rejsampleEVMA2R5$unadj.values)
#  1,] 1.0000000 0.4742813
#[2,] 0.4742813 1.0000000
  cor(samplemixgllimPlot)
 # [,1]         [,2]
 # [1,]  1.000000000 -0.006872794
 # [2,] -0.006872794  1.000000000
cor(saabc.MA2150$post.sample)
#[1,]  1.0000000 -0.1123789
#[2,] -0.1123789  1.0000000

cor(rejsampleACF2MA2$unadj.values)
#[1,]  1.0000000 -0.0448014
#[2,] -0.0448014  1.0000000
cor(saabc.ACF2MA2150$post.sample)
#[1,] 1.00000000 0.01833507
#[2,] 0.01833507 1.00000000


#######################################################################
##
##  Compute TRUE posterior means  by numerical integration
##
#######################################################################

 margpdf1<-postpdfMA2150$marg1
 margpdf2<-postpdfMA2150$marg2
 meantheta1<-integrate.xy(thetaseq1, thetaseq1*margpdf1)
 meantheta2<-integrate.xy(thetaseq2, thetaseq2*margpdf2)
 # 0.6348032 and  meantheta2= 0.2027278
 stdtheta1<-sqrt(integrate.xy(thetaseq1, (thetaseq1-meantheta1)^2*margpdf1))
 stdtheta2<-sqrt(integrate.xy(thetaseq2, (thetaseq2-meantheta2)^2*margpdf2))
 # stdtheta1= 0.0801874 and stdtheta2= 0.07540907
 
 ######################################################################
 ##
 ## Compute TRUE posterior means, std, cor by importance sampling with
 ## gllim mixture as proposal 
 ##
 ######################################################################
 
 ## sample from GLLiM mixture for target: 10000 or 10^5 points
 mixparam<-ComputeGllimFullMixtParamIID(modgllimMA2R5IIDFullK30,R=5,t(ytargetMA2150),0)
 samplemixgllim<-simdataset(100000,mixparam$Pi, t(mixparam$Mu),mixparam$S)$X
 compes<-GllimISMA2(samplemixgllim,mixparam$Pi, mixparam$Mu,mixparam$S,ytargetMA2150)

 #TRUE values:
 compes$ISmean
 #[1] 0.6348713 0.2027734 # close to numerical integration ones
  diag(sqrt(compes$ISvar))
  #0.08024054 0.07577829 close to numerical integration ones
  compes$IScor
 #[1,] 1.000000 0.472139
 #[2,] 0.472139 1.000000
 
 
  


  ######################################################################################
  ##
  ## Bivariate Beta model like in ES paper
  ## L=5, D=2 R=500 N=10^5
  ## simulations initially with R=500 then reduced to R=100 for tests in GLLiM-ABC paper
  ##
  ######################################################################################
  ## Data et al generated via BBM.R  with R=500
  ## from load(dataBBM.Rdata)
  ######################################################################################
  

  ######################################################################################
  # K=100 et Covariance with full blocs and reduced iR to 100, ie iRD=200
  #
  # max iter set to 150 in gllimIID.R  iR=100
  #
  system.time(modgllimBBMIIDFullK100n100 <- GllimFitIID(bbmparams, bbmdata[1:200,],100, K=100, constr=list(Sigma="")))
  #   user   system  elapsed 
  # 30019.14 11256.22 41419.40  11.5 h  cv en 80+ iterations
  #
  ## Second run, first one lost
  # 40927.69 14752.16 55853.74 15.51 h cv en 120 iterations! 
  #
  save(list = c("modgllimBBMIIDFullK100n100"), file = "modgllimBBMIIDFullK100n100.RData")
  #
  ## third run with  K=100
  system.time(modgllimBBMIIDFullK100n100Run3 <- GllimFitIID(bbmparams, bbmdata[1:200,],100, K=100, constr=list(Sigma="")))
  # 30123.05 10190.69 40414.79 --> 11.22 h -- 85 iter
  #
  save(list = c("modgllimBBMIIDFullK100n100Run3"), file = "modgllimBBMIIDFullK100n100Run3.RData")
  

  ########### 
  # Check gllim mixture for bbmobs all theta at 1 or bbmobs2 with all true theta to 2
  ###########

  mixparamBBM<-ComputeGllimFullMixtParamIID(modgllimBBMIIDFullK100n100,R=100,bbmobs[1:200,],0)
  samplemixgllimBBM<-simdataset(1000,mixparamBBM$Pi, t(mixparamBBM$Mu),mixparamBBM$S)$X

  
  # check for all theta=2
  mixparamBBM<-ComputeGllimFullMixtParamIID(modgllimBBMIIDFullK100n100,R=100,bbmobs2[1:200,],0)
  samplemixgllimBBM<-simdataset(1000,mixparamBBM$Pi, t(mixparamBBM$Mu),mixparamBBM$S)$X
  
  
    
  ################################
  ### Compute MW2 and L2 distances
  ################################
  # K=100 n=100 for bbmobs[1:200,] (all theta=1)
  #
  system.time(distMW2L2BBMFullK100n100obs1<-
                PreCompDistFullParIID(modgllimBBMIIDFullK100n100,R=100,bbmobs[1:200,],bbmdata[1:200,],wthr = 10^(-10))
  )
  # 94038.429   454.763 12511.779 3.4 h
  save(list = c("distMW2L2BBMFullK100n100obs1"), file = "distMW2L2BBMFullK100n100obs1.RData")
  
  ### Same setting but for second Gllim run adn with a trheshold at 10^-4:
  system.time(distMW2L2BBMFullK100n100obs1Run2<-
                PreCompDistFullParIID(modgllimBBMIIDFullK100n100,R=100,bbmobs[1:200,],bbmdata[1:200,],wthr = 10^(-4))
  )
  # 16838.158    84.982  2284.880 --38 min
  save(list = c("distMW2L2BBMFullK100n100obs1Run2"), file = "distMW2L2BBMFullK100n100obs1Run2.RData")
  
  ################################
  ### Select samples from distances 
  # 0.001 for 100 points, 0.01 for 1000
  # 0.0005 for 50 pts
 
  rejsampleBBM<-RejABCMW2L2(0.0005,distMW2L2BBMFullK100n100obs1[1,],distMW2L2BBMFullK100n100obs1[2,],bbmparams)
  # K100 n100 Run2
  # rejsampleBBM<-RejABCMW2L2(0.0005,distMW2L2BBMFullK100n100obs1Run2[1,],distMW2L2BBMFullK100n100obs1Run2[2,],bbmparams)
  

  ################################
  ### compute gllim posterior expectation and log variances in the IID case 
  ################################ 
  # K100 n 100 
  # 
  system.time(ExpVarMatzBBMFullK100n100<-ExpLogPostVarGllimFullIID(modgllimBBMIIDFullK100n100,bbmdata[1:200,],100) )
  # 
  # t(t()) to keep the matrix format...
  ExpVarTargetBBMFullK100n100<-ExpLogPostVarGllimFullIID(modgllimBBMIIDFullK100n100,t(t(bbmobs[1:200,])),100) 
  #
   L=5
  espysimu2BBMIID<-ExpVarMatzBBMFullK100n100[,1:L]
  # espytarget is 1 x L 
  espytargetBBMIID<-ExpVarTargetBBMFullK100n100[1:L]
  #
  ################################
  ### GLLiM-E-ABC and 
  rejsampleEBBM <-abc(target=espytargetBBMIID, param=t(bbmparams), sumstat=espysimu2BBMIID, tol=.0005, method ="rejection")
  ### GLLiM-EV-ABC, exp + log var 
  rejsampleEVBBM <-abc(target=ExpVarTargetBBMFullK100n100, param=t(bbmparams), sumstat=ExpVarMatzBBMFullK100n100, tol=.0005, method ="rejection")
  
  
    
  ##################################################################################
  ## SEMI AUTO ABC (SA) done with 7 quantiles for first dimension, concatenated with 7 quantiles
  ## for second dimension
  ##################################################################################
  
  ## 
  ## Compute quantiles with function below
  ##
   qprob=seq(0,1,by=1/6)
   # for the reduced target: bbmobs[1:200,]

   QuantTransf<-function(dataRDN,qprob){
     Nsim<-dim(dataRDN)[[2]] # N ou 1
     ql<-length(qprob)
     quantmat<-matrix(0,2*ql,Nsim)
     for (i in 1:Nsim){
       temp<-matrix(dataRDN[,i],nrow=2)  #length 2R turned into a 2 x R matrix
       quantmat[,i]<-c(quantile(temp[1,], qprob), quantile(temp[2,],qprob))
     }
     quantmat
   }
   # ex 
   quantbbmobs<-QuantTransf(bbmobs,qprob) #14 x1
   # ex 
   quantbbmdata<-QuantTransf(bbmdata,qprob) # 14x N
   
  ##
  ## SA on data of size 150 reduced to 14 quantiles + their powers 
  ## 
  ##################################################################################
  mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
  saabc.BBM<-semiauto.abc(obs=t(quantbbmobs),param=t(bbmparams), sumstats=t(quantbbmdata),satr=mytf, tol=0.0005, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)
  


  ####################################################################
  #### GLLiM-ABC with quantiles, K=40 and standard GLLiM (non iid)
  #### Full covariances
  ####################################################################

  K=40
  system.time(modgllimQuantFullK40 <- GllimFit(bbmparams, quantbbmdata, K))
  #user   system  elapsed 
  #450.613 168.960 620.412 
  save(list = c("modgllimQuantFullK40"), file = "modgllimQuantGFullK40.RData")
  #

  # check gllim mixture quality: 
  mixparamQuant<-ComputeGllimFullMixtParam(modgllimQuantFullK40,quantbbmobs,0)
  samplemixgllimQuant<-simdataset(50,mixparamQuant$Pi, t(mixparamQuant$Mu),mixparamQuant$S)$X
  
  # Compute MW2 et L2 dist
  system.time(distMW2L2QuantFullK40<-
                PreCompDistFullPar(modgllimQuantFullK40,quantbbmobs,quantbbmdata,wthr = 0.0001)
  )
  #4531.449 1287.609  927.275
  rejsampleQuant<-RejABCMW2L2(0.0005,distMW2L2QuantFullK40[1,],distMW2L2QuantFullK40[2,],bbmparams)
  # 

  ################################
  # GLLiM-E-ABC and GLLiM-EV-ABC
  ################################

  # espyQuant is of size N x L 
  espyQuant<-t(gllim_inverse_map(quantbbmdata,modgllimQuantFullK40$mod)$x_exp)
  
  # sum stat = posterior expectation for the target (obs or ytest)
  # ytarget is 1 x D
  # espytarget is 1 x L 
  #
  # 1) inversion of the simulated observation
  espytargetQuant=t(gllim_inverse_map(quantbbmobs,modgllimQuantFullK40$mod)$x_exp)
  # 
  rejsampleEQuant <-abc(target=espytargetQuant, param=t(bbmparams), sumstat=espyQuant, tol=.0005, method ="rejection")
  
  ### EV  needs a specific function not done yet
  ##
  # simulated observation
  logvarytargetQuant<-LogPostVarGllimFull(modgllimQuantFullK40,quantbbmobs)
  
  system.time(logvaryQuant<-LogPostVarGllimFull(modgllimQuantFullK40, quantbbmdata))
  #user  system elapsed 
  #184.613   4.454 190.320
  
  ####### Gllim-EV-ABC
  #ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
  # rejection is the one used for wasserstein ...
  # tol =0.0005 select 50 samples when Nsimu =10^5]
  
  # thetasimu2 is L x M
  # selected sample of size L x 100 is in rejsampleEV$unadj.values
  # esplogytarget  1 x 2L 
  # esplogysimu2   M x 2L
  esplogytargetQuant<-c(espytargetQuant, t(logvarytargetQuant))
  esplogyQuant<-cbind(espyQuant, t(logvaryQuant))
  
  # threshold 0.01---> 1000 samples
  rejsampleEVQuant <-abc(target=esplogytargetQuant, param=t(bbmparams), sumstat=esplogyQuant, tol=.0005, method ="rejection")
  
  
  
  #############################################
  ## PLOTS posterior margins
  #############################################

  # margin 1 theta1
  plot(density(rejsampleBBM$MW2postval[1,], from=0,to=5), xlim=c(0,5), ylim=c(0,0.8), col="black", main="Theta1")
  lines(density(rejsampleEBBM$unadj.values[,1], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVBBM$unadj.values[,1], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density( saabc.BBM$post.sample[,1] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleBBM$L2postval[1,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimBBM[,1], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  # margin 2 theta2
  plot(density(rejsampleBBM$MW2postval[2,], from=0,to=5), xlim=c(0,5), ylim=c(0,0.83), col="black", main="Theta2")
  lines(density(rejsampleEBBM$unadj.values[,2], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVBBM$unadj.values[,2], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density( saabc.BBM$post.sample[,2] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleBBM$L2postval[2,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimBBM[,2], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  # margin 3 theta3
  plot(density(rejsampleBBM$MW2postval[3,], from=0,to=5), xlim=c(0,5), ylim=c(0,1.2), col="black", main="Theta3")
  lines(density(rejsampleEBBM$unadj.values[,3], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVBBM$unadj.values[,3], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density( saabc.BBM$post.sample[,3] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleBBM$L2postval[3,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimBBM[,3], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  # margin 4 theta4
  plot(density(rejsampleBBM$MW2postval[4,], from=0,to=5), xlim=c(0,5), ylim=c(0,1), col="black", main="Theta4")
  lines(density(rejsampleEBBM$unadj.values[,4], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVBBM$unadj.values[,4], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density( saabc.BBM$post.sample[,4] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleBBM$L2postval[4,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimBBM[,4], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  # margin 5 theta5
  plot(density(rejsampleBBM$MW2postval[5,], from=0,to=5), xlim=c(0,5), ylim=c(0,1.1), col="black", main="Theta5")
  lines(density(rejsampleEBBM$unadj.values[,5], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVBBM$unadj.values[,5], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density(saabc.BBM$post.sample[,5] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleBBM$L2postval[5,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimBBM[,5], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  


   #######################################################
   # Empirical means and rmse for each ABC samples etc.
   #######################################################

  ## Empirical means for each ABC samples etc....
  rowMeans(rejsampleBBM$MW2postval)
  #[1] 1.211681 1.319846 0.871997 1.004025 1.235376
   rowMeans(rejsampleBBM$L2postval)
  #[1] 1.597159 1.700024 1.534706 1.627089 1.827051
   colMeans(rejsampleEBBM$unadj.values)
  #[1] 1.1420153 1.8989545 0.8717550 0.7867029 1.4726572
   colMeans(rejsampleEVBBM$unadj.values)
  #[1] 0.9900871 1.8675065 0.7462293 0.5943075 1.3852693
   colMeans(saabc.BBM$post.sample)
  #[1] 0.7708875 0.8253933 0.9472275 0.7567767 0.9177751
   colMeans(samplemixgllimBBM)
  #[1] 1.5049117 1.7360966 0.8904948 0.9896936 1.6168874
   
   ## # compute RMSE
   Mm<-50 #pts
   truemat<-matrix(1,5,Mm)
   msemat<-(rejsampleBBM$MW2postval -truemat)^2
   rmsem<-sqrt(rowMeans(msemat))
   #W2: 0.7903403 0.8200935 0.4386980 0.5236739 0.4262541
   #L2: 1.1444669 1.2243447 0.9681786 1.1166781 1.2953570
   msemat<-(t(rejsampleEBBM$unadj.values) -truemat)^2
   rmsem<-sqrt(rowMeans(msemat))
   # E: 0.6773737 1.1806217 0.4977470 0.5680962 0.6259876
   # EV: 0.5053208 1.0765529 0.4682663 0.5618459 0.5168076
   # SA: 0.4925638 0.4725109 0.5238772 0.4684273 0.5228873
   # Mix: 0.9257485 1.2764569 0.8247102 0.8478881 1.0208375
   

   #######################################################
   # Empirical means and rmse - 14 quantiles case
   #######################################################

   rowMeans(rejsampleQuant$MW2postval)
   #[1] 1.2574198 0.9094898 0.9219186 1.0429116 0.9499574
   rowMeans(rejsampleQuant$L2postval)
   #[1] 3.390369 3.009987 3.467197 3.361721 2.653772
    colMeans(rejsampleEQuant$unadj.values)
   #[1] 1.2660946 0.9056010 0.8727178 1.1055710 1.0822442
    colMeans(rejsampleEVQuant$unadj.values)
   #[1] 1.5301143 0.8524538 1.0950464 0.7271341 0.9045728
   colMeans(saabc.BBM$post.sample)
   #[1] 0.7708875 0.8253933 0.9472275 0.7567767 0.9177751
   colMeans(samplemixgllimQuant)
   #[1] 0.4489356 0.8585833 0.7395903 0.5519630 0.5777553
   
   ### compute RMSE
   Mm<-50 #pts
   truemat<-matrix(1,5,Mm)
   msemat<-(rejsampleQuant$MW2postval -truemat)^2
   rmsem<-sqrt(rowMeans(msemat))
   #W2: 0.5726653 0.5000068 0.5240591 0.5280579 0.4638565
   # L2: 2.747297 2.492395 2.692861 2.731984 1.994603
   msemat<-(t(rejsampleEQuant$unadj.values) -truemat)^2
   rmsem<-sqrt(rowMeans(msemat))
   # E: 0.6286144 0.5011295 0.5496683 0.5411596 0.5135395
   # EV: 0.8083779 0.5040375 0.5774495 0.5649924 0.4497382
   # SA: 0.4925638 0.4725109 0.5238772 0.4684273 0.5228873
   # mix: 1.685064 1.464130 1.367335 1.299079 1.213878
   
  #########################################################
  ##.           Plots posterior margins - 14 Quantiles case
  #########################################################
  ## margins
  # margin 1 theta1
  plot(density(rejsampleQuant$MW2postval[1,], from=0,to=5), xlim=c(0,5), ylim=c(0,0.7), col="black", main="Theta1")
  lines(density(rejsampleEQuant$unadj.values[,1], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVQuant$unadj.values[,1], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density( saabc.BBM$post.sample[,1] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleQuant$L2postval[1,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimQuant[,1], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  # margin 2 theta2
  plot(density(rejsampleQuant$MW2postval[2,], from=0,to=5), xlim=c(0,5), ylim=c(0,0.83), col="black", main="Theta2")
  lines(density(rejsampleEQuant$unadj.values[,2], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVQuant$unadj.values[,2], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density( saabc.BBM$post.sample[,2] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleQuant$L2postval[2,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimQuant[,2], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  # margin 3 theta3
  plot(density(rejsampleQuant$MW2postval[3,], from=0,to=5), xlim=c(0,5), ylim=c(0,0.75), col="black", main="Theta3")
  lines(density(rejsampleEQuant$unadj.values[,3], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVQuant$unadj.values[,3], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density( saabc.BBM$post.sample[,3] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleQuant$L2postval[3,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimQuant[,3], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  # margin 4 theta4
  plot(density(rejsampleQuant$MW2postval[4,], from=0,to=5), xlim=c(0,5), ylim=c(0,0.85), col="black", main="Theta4")
  lines(density(rejsampleEQuant$unadj.values[,4], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVQuant$unadj.values[,4], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density( saabc.BBM$post.sample[,4] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleQuant$L2postval[4,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimQuant[,4], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  # margin 5 theta5
  plot(density(rejsampleQuant$MW2postval[5,], from=0,to=5), xlim=c(0,5), ylim=c(0,0.8), col="black", main="Theta5")
  lines(density(rejsampleEQuant$unadj.values[,5], from=0,to=5), xlim=c(0,5), col="red")
  lines(density(rejsampleEVQuant$unadj.values[,5], from=0,to=5), xlim=c(0,5), col="red", lty=2)
  lines(density(saabc.BBM$post.sample[,5] , from=0,to=5), xlim=c(0,5), col="green")
  lines(density(rejsampleQuant$L2postval[5,], from=0,to=5), xlim=c(0,5), col="blue")
  #lines(thetaseq1,postpdfMA2150$marg1, col="purple")
  lines(density(samplemixgllimQuant[,5], from=0,to=5), xlim=c(0,5), col="green", lty=2)
  #
  abline(v=1, lty=2)
  
  ####################################################################################
  ### Compute true posterior means by Importance Sampling is not possible because the likelihood
  # is not available
  # TODO instead:  compare to true values used for simulation, like Table 3 in the ES paper
  ####################################################################################

  
    
  
  
  


