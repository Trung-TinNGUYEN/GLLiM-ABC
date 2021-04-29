##############################################################################################
####
#### GLLiM+ABC 
#### 
#### F. Forbes 
#### Update April 29, 2021: Faster L2 & MW2 computation with PreCompDistIsoPar
#### (illustrated on the MA(2) example)
##############################################################################################


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

##############################################################################################
####  Example Sum of 2 MA(1) with opposite theta, D=10, L=1
##############################################################################################


### Step 1: Setting: target y (observations to be inverted, dimension is D) ;
###         Data sets simulation, one for GLLiM, one for ABC if wanted

# output: ytarget; (thetasimu1, ysimu1) size N; (thetasimu2, ysimu2) size M ;
#  true posterior if possible

# After having run the SMA1.R code to produce the data etc...
load("dataSMA1.RData")
# recover c("ytargetSMA1","ysimu1SMA1","thetasimu1SMA1","ysimu2SMA1","thetasimu2SMA1")


### Step 2: Fitting GLLiM on (thetasim1, ysim1), model selection (K, constraints on Sigmas) 
### and computing required quantities not directly provided by xLLiM
###########################################################################################################
# xllim  K (chosen manually here) gaussian components  and  (default) isotropic constraints on the Sigma_k
# data and param should be D x N and L x N
# output: a list with the output of xllim + inverse Gammak and Sigmakstar

##### application SMA1
K=20
system.time(modgllimSMA1 <- GllimIsoFit(thetasimu1SMA1, ysimu1SMA1, K))
#user  system elapsed 
#264.904  68.117 334.399 --> 5.6 minutes


### Step 3.1: GLLiM-E and EV using package abc for the ABC scheme 
###
#################################################################################################
# ABC avec sum stat from gllim: 1) espectation
##################################
# summary stat = posterior expectations computed from gllim for the training examples
# ysimu2 is D x M 
# modg is the output of GllimIsoFit
# espysimu2 is of size M x L 
espysimu2SMA1<-t(gllim_inverse_map(ysimu2SMA1,modgllimSMA1$mod)$x_exp)


# sum stat = posterior espectation for the target (obs or y test)
# ytarget is 1 x D
# espytarget is 1 x L 
espytargetSMA1=t(gllim_inverse_map(t(ytargetSMA1),modgllimSMA1$mod)$x_exp)

#ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for wasserstein ...
# tol =0.001 select 100 samples when Nsimu =10^5]

# thetasimu2 is L x M
# selected sample of size L x 100 is in rejsampleESMA1$unadj.values
rejsampleESMA1 <-abc(target=espytargetSMA1, param=t(thetasimu2SMA1), sumstat=espysimu2SMA1, tol=.001, method ="rejection")

# example: rejsampleSMA2 <-abc(target=espytargetSMA2, param=t(thetasimu2SMA2), sumstat=espysimu2SMA2, tol=.001, method ="rejection")
# plot(rejsampleSMA2$unadj.values, xlim=c(-2,2),ylim=c(-1,1), pch=20, cex=.5)
# 
##################################
# ABC avec sum stat from gllim: 2) espectation + log variances
##################################

logvarytargetSMA1<-LogPostVarGllimIso(modgllimSMA1,t(ytargetSMA1))

system.time(logvarysimu2SMA1<-LogPostVarGllimIso(modgllimSMA1, ysimu2SMA1))
#user  system elapsed 
#29.764   0.549  30.479  in seconds

####### Gllim-EV
#ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for wasserstein ...
# tol =0.001 select 100 samples when Nsimu =10^5]

# thetasimu2 is L x M
# selected sample of size L x 100 is in rejsampleEV$unadj.values
# esplogytarget  1 x 2L 
# esplogysimu2   M x 2L
esplogytargetSMA1<-c(espytargetSMA1, t(logvarytargetSMA1))
esplogysimu2SMA1<-cbind(espysimu2SMA1, t(logvarysimu2SMA1))

rejsampleEVSMA1 <-abc(target=esplogytargetSMA1, param=t(thetasimu2SMA1), sumstat=esplogysimu2SMA1, tol=.001, method ="rejection")

# Example
#esplogytarget<-c(espytargetSMA2, t(logvarytarget))
#esplogysimu2<-cbind(espysimu2SMA2, t(logvarysimu2))
#rejsampleEV <-abc(target=esplogytarget, param=t(thetasimu2SMA2), sumstat=esplogysimu2, tol=.001, method ="rejection")

#plot(rejsampleEV$unadj.values, xlim=c(-2,2),ylim=c(-1,1), pch=20, cex=.5)

## SEMI AUTO ABC
mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
saabc.SMA1<-semiauto.abc(obs=ytargetSMA1,param=t(thetasimu2SMA1), sumstats=t(ysimu2SMA1),satr=mytf, tol=0.001, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)


#################################################################################################
###
### Step 3.2: computing distances, L2, MW2 
###
#################################################################################################
# Cases : L=1 or L>1 (to be unified later?) ; L2, MW2 

####################
##   SMA1 is a Case L=1     ##
####################
## for L=1 distances are computed one at a time with dopar
## TODO : could be faster to compute all distances at once like in ComputeDistMargin
##        rather than usind dopar and the cores

numCores <- detectCores() #numCores = 8
registerDoParallel(numCores)

### Possibility to reduce the mixtures before computing distances to save time
Wthresh<-0.001
#M=10^5
M<-dim(ysimu2SMA1)[2]

system.time(distMW2L2SMA1<- foreach (i=1:M, .combine=cbind) %dopar% {
  Compute1DistL1Iso(modgllimSMA1,t(ytargetSMA1),matrix(ysimu2SMA1[,i]),Wthresh)
})
#user   system  elapsed 
#3027.860   61.114  437.567 --> 7.3 minutes

# Comparison without dopar no core used
system.time(distMW2L2SMA1test<-
  ComputeDistL1Iso(modgllimSMA1,t(ytargetSMA1),ysimu2SMA1,Wthresh)
)
#user   system  elapsed 
#1184.802  237.887 1424.786 ---> 23 min

# Comparison with accelerated version and dopar
system.time(distMW2L2SMA1test<-
              PreCompDistL1IsoPar(modgllimSMA1,t(ytargetSMA1),ysimu2SMA1,Wthresh)
)
#user   system  elapsed 
#1986.699  554.566  367.248 ---> 6 min (vs 7 a bit faster)

### Step 4: ABC scheme, default rejection ABC using distance quantiles
#################################################################################################

rejsampleSMA1<-RejABCMW2L2(0.01,distMW2L2SMA1[1,],distMW2L2SMA1[2,],thetasimu2SMA1)

# return(list("MW2postval"= MW2postval,"L2postval"= L2postval ))

#### Plot 
plot(density(rejsampleSMA1$L2postval, from=-2,to=2), xlim=c(-2,2), ylim=c(0,0.43), col="blue", main="rho")
lines(density(rejsampleESMA1$unadj.values, from=-2,to=2), xlim=c(-2,2), col="red")
lines(density(rejsampleEVSMA1$unadj.values, from=-2,to=2), xlim=c(-2,2), col="red", lty=2)
lines(density( saabc.SMA1$post.sample , from=-2,to=2), xlim=c(-2,2), col="green")
lines(density(rejsampleSMA1$MW2postval, from=-2,to=2), xlim=c(-2,2), col="black")
lines(thetaseqSMA1,postpdfSMA1, col="purple")
# true param is 1 here 
abline(v=1, lty=2)
abline(v=-1, lty=2)

##############################################################################################
##
## Example Sum of 2 MA(2) with same theta2 and opposite theta1, D=10, L=2
##
##############################################################################################
### Step 1: Setting: target y (observations to be inverted, dimension is D) ;
###         Data sets simulation, one for GLLiM, one for ABC if wanted
# output: ytarget; (thetasimu1, ysimu1) size N; (thetasimu2, ysimu2) size M ;
#  true posterior if possible

load("dataSMA2.RData")
# recover c("ytargetSMA2","ysimu1SMA2","thetasimu1SMA2","ysimu2SMA2","thetasimu2SMA2","ysimu3SMA2","thetasimu3SMA2" )

### Step 2: Fitting GLLiM on (thetasim1, ysim1), model selection (K, constraints on Sigmas) 
### and computing required quantities not directly provided by xLLiM
###########################################################################################################
# xllim  K (chosen manually here) gaussian components  and  (default) isotropic constraints on the Sigma_k
# data and param should be D x N and L x N
# output: a list with the output of xllim + inverse Gammak and Sigmakstar

##### application SMA2
K=80
system.time(modgllimSMA2 <- GllimIsoFit(thetasimu1SMA2, ysimu1SMA2, K))
#  user  system elapsed 
#915.988  327.920 1248.644 --> 21 min

save(list = c("modgllimSMA2"), file = "modgllimSMA2K80.RData")
##################################################################################

### Step 3.1: GLLiM-E and EV using package abc for the ABC scheme 
###
#################################################################################################
# ABC avec sum stat from gllim: 1) espectation
##################################
# summary stat = posterior expectations computed from gllim for the training examples
# ysimu2 is D x M 
# modg is the output of GllimIsoFit
# espysimu2 is of size M x L 
espysimu2SMA2<-t(gllim_inverse_map(ysimu2SMA2,modgllimSMA2$mod)$x_exp)

# M=10^6 : espysimu3 is of size M1 x L 
espysimu3SMA2<-t(gllim_inverse_map(ysimu3SMA2,modgllimSMA2$mod)$x_exp)

# sum stat = posterior espectation for the target (obs or y test)
# ytarget is 1 x D
# espytarget is 1 x L 
espytargetSMA2=t(gllim_inverse_map(t(ytargetSMA2),modgllimSMA2$mod)$x_exp)

#ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for wasserstein ...
# tol =0.001 select 100 samples when Nsimu =10^5]

# thetasimu2 is L x M
# selected sample of size L x 100 is in rejsampleESMA1$unadj.values
# Thrshold 0.01 for 1000 samples
rejsampleESMA2 <-abc(target=espytargetSMA2, param=t(thetasimu2SMA2), sumstat=espysimu2SMA2, tol=.01, method ="rejection")

# thetasimu3 is L x M1
# selected sample of size L x 1000 is in rejsampleESMA1$unadj.values
# Thrshold 0.001 for 1000 samples
rejsampleESMA2M106 <-abc(target=espytargetSMA2, param=t(thetasimu3SMA2), sumstat=espysimu3SMA2, tol=.001, method ="rejection")


# example: rejsampleSMA2 <-abc(target=espytargetSMA2, param=t(thetasimu2SMA2), sumstat=espysimu2SMA2, tol=.001, method ="rejection")
# plot(rejsampleSMA2$unadj.values, xlim=c(-2,2),ylim=c(-1,1), pch=20, cex=.5)
# 
##################################
# ABC avec sum stat from gllim: 2) espectation + log variances
##################################

logvarytargetSMA2<-LogPostVarGllimIso(modgllimSMA2,t(ytargetSMA2))

system.time(logvarysimu2SMA2<-LogPostVarGllimIso(modgllimSMA2, ysimu2SMA2))
#user  system elapsed 
#130.618   2.404 133.693 -->2.2 min

# M1=10^6
system.time(logvarysimu3SMA2<-LogPostVarGllimIso(modgllimSMA2, ysimu3SMA2))
#user   system  elapsed 
#1616.417   19.441 1638.812 ---> 27 min

####### Gllim-EV
#ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for wasserstein ...
# tol =0.001 select 100 samples when Nsimu =10^5]

# thetasimu2 is L x M
# selected sample of size L x 100 is in rejsampleEV$unadj.values
# esplogytarget  1 x 2L 
# esplogysimu2   M x 2L
esplogytargetSMA2<-c(espytargetSMA2, t(logvarytargetSMA2))
esplogysimu2SMA2<-cbind(espysimu2SMA2, t(logvarysimu2SMA2))

# threshold 0.01---> 1000 samples
rejsampleEVSMA2 <-abc(target=esplogytargetSMA2, param=t(thetasimu2SMA2), sumstat=esplogysimu2SMA2, tol=.01, method ="rejection")

# with M1=10^6 and epsilon =0.001
esplogysimu3SMA2<-cbind(espysimu3SMA2, t(logvarysimu3SMA2))
rejsampleEVSMA2M106 <-abc(target=esplogytargetSMA2, param=t(thetasimu3SMA2), sumstat=esplogysimu3SMA2, tol=.001, method ="rejection")


#plot(rejsampleEV$unadj.values, xlim=c(-2,2),ylim=c(-1,1), pch=20, cex=.5)

## SEMI AUTO ABC-- TODO a fair comparison in the usage of data 1 and 2 
mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
# DONT Forget to change: for SMA2 threshold is 0.01
saabc.SMA2<-semiauto.abc(obs=ytargetSMA2,param=t(thetasimu2SMA2), sumstats=t(ysimu2SMA2),satr=mytf, tol=0.01, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)
# with M1=10^6 and epsilon =0.001
saabc.SMA2M106<-semiauto.abc(obs=ytargetSMA2,param=t(thetasimu3SMA2), sumstats=t(ysimu3SMA2),satr=mytf, tol=0.001, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)

#################################################################################################
###
### Step 3.2: computing distances, L2, MW2 
###
#################################################################################################
# Cases : L=1 or L>1 (to be unified later?) ; L2, MW2 

###########################################################################
##   SMA2 is a Case L=2 : 1) compute post margins , 2) compute full post     
##########################################################################

######### 1) Compute posterior margins is much FASTER

### Possibility to reduce the mixtures before computing distances to save time
Wthresh<-0.001


## ComputeDistMargin: Compute all z without foreach faster? similar in speed
## Not faster
#system.time(distMW2L2SMA2marg12test<-
#                ComputeDistMarginIso(modgllimSMA2,t(ytargetSMA2),ysimu2SMA2,wthr = 0.001))
##user  system elapsed 
##Timing stopped at: 2786 355.4 3150---> 52 min ???

# To compare in speed with 
numCores <- detectCores() #numCores = 8
registerDoParallel(numCores)

M<-10^5
system.time(distMW2L2SMA2marg12<- foreach (i=1:M, .combine=cbind) %dopar% {
  Compute1DistMarginIso(modgllimSMA2,t(ytargetSMA2),matrix(ysimu2SMA2[,i]),wthr = 0.001)
})
#user  system elapsed 
#15663.188   266.815  2158.846 ---> 36 min



save(list = c("distMW2L2SMA2marg12"), file = "distSMA2marg12K80.RData")

############################
# M<-10^6 plus change to thetasimu3 and ysimu3
M1<-10^6
system.time(distMW2L2SMA2marg12M106<- foreach (i=1:M1, .combine=cbind) %dopar% {
  Compute1DistMarginIso(modgllimSMA2,t(ytargetSMA2),matrix(ysimu3SMA2[,i]),wthr = 0.001)
})
#user  system elapsed 
#159222.829   2545.986  21829.198 --> 6 heures

save(list = c("distMW2L2SMA2marg12M106"), file = "distSMA2marg12K80M106.RData")


### True posterior if not computed and restored from before
thetaseq1<-seq(-2,2,length=400)
#thetaseq2<-seq(-1,1,length=400) # provides plots with bumps??
thetaseq2<-seq(-1,1,length=200) # no bumps 
# function PostPdfSMA2 is in SMA2.R
postpdfSMA2<-PostPdfSMA2(thetaseq1,thetaseq2, ytargetSMA2)
# output list("marg1"=margpdf1,"marg2"=margpdf2,"postdf"= full_df)


### Step 4: ABC scheme, default rejection ABC using distance quantiles
#################################################################################################

MW2SMA2m1<-distMW2L2SMA2marg12[1,]
MW2SMA2m2<-distMW2L2SMA2marg12[2,]
L2SMA2m1<-distMW2L2SMA2marg12[3,]
L2SMA2m2<-distMW2L2SMA2marg12[4,]

# 1st margin
rejsampleSMA2m1<-RejABCMW2L2(0.01,MW2SMA2m1,L2SMA2m1,t(thetasimu2SMA2[1,]))
# 2nd margin
rejsampleSMA2m2<-RejABCMW2L2(0.01,MW2SMA2m2,L2SMA2m2,t(thetasimu2SMA2[2,]))

# return(list("MW2postval"= MW2postval,"L2postval"= L2postval ))


#### Plots of posterior margins from distances on surrogate posterior margins
#################################################################################################

# theta1 1st margin
plot(density(rejsampleSMA2m1$L2postval, from=-2,to=2), xlim=c(-2,2), ylim=c(0,0.53), col="blue", main="theta1")
lines(density(rejsampleESMA2$unadj.values[,1], from=-2,to=2), xlim=c(-2,2), col="red")
lines(density(rejsampleEVSMA2$unadj.values[,1], from=-2,to=2), xlim=c(-2,2), col="red", lty=2)
lines(density( saabc.SMA2$post.sample[,1] , from=-2,to=2), xlim=c(-2,2), col="green")
lines(density(rejsampleSMA2m1$MW2postval, from=-2,to=2), xlim=c(-2,2), col="black")
lines(thetaseq1,postpdfSMA2$marg1, col="purple")
# true param is 0.7 here 
abline(v=0.7, lty=2)
abline(v=-0.7, lty=2)

# theta2 2nd margin
plot(density(rejsampleSMA2m2$L2postval, from=-1,to=1), xlim=c(-1,1), ylim=c(0,1.1), col="blue", main="theta2")
lines(density(rejsampleESMA2$unadj.values[,2], from=-1,to=1), xlim=c(-1,1), col="red")
lines(density(rejsampleEVSMA2$unadj.values[,2], from=-1,to=1), xlim=c(-1,1), col="red", lty=2)
lines(density( saabc.SMA2$post.sample[,2] , from=-1,to=1), xlim=c(-1,1), col="green")
lines(density(rejsampleSMA2m2$MW2postval, from=-1,to=1), xlim=c(-1,1), col="black")
lines(thetaseq2,postpdfSMA2$marg2, col="purple")
# true param is 0.5 here 
abline(v=0.5, lty=2)


##### M1=10^6 and thetasimu3 etc...
MW2SMA2m1M106<-distMW2L2SMA2marg12M106[1,]
MW2SMA2m2M106<-distMW2L2SMA2marg12M106[2,]
L2SMA2m1M106<-distMW2L2SMA2marg12M106[3,]
L2SMA2m2M106<-distMW2L2SMA2marg12M106[4,]

# the threshold can be changed to 0.001 (0.1% quantile)
# 1st margin
rejsampleSMA2m1M106<-RejABCMW2L2(0.001,MW2SMA2m1M106,L2SMA2m1M106,t(thetasimu3SMA2[1,]))
# 2nd margin
rejsampleSMA2m2M106<-RejABCMW2L2(0.001,MW2SMA2m2M106,L2SMA2m2M106,t(thetasimu3SMA2[2,]))

# theta1 1st margin
plot(density(rejsampleSMA2m1M106$L2postval, from=-2,to=2), xlim=c(-2,2), ylim=c(0,0.53), col="blue", main="theta1")
lines(density(rejsampleESMA2M106$unadj.values[,1], from=-2,to=2), xlim=c(-2,2), col="red")
lines(density(rejsampleEVSMA2M106$unadj.values[,1], from=-2,to=2), xlim=c(-2,2), col="red", lty=2)
lines(density( saabc.SMA2M106$post.sample[,1] , from=-2,to=2), xlim=c(-2,2), col="green")
lines(density(rejsampleSMA2m1M106$MW2postval, from=-2,to=2), xlim=c(-2,2), col="black")
lines(thetaseq1,postpdfSMA2$marg1, col="purple")
# true param is 0.7 here 
abline(v=0.7, lty=2)
abline(v=-0.7, lty=2)

# theta2 2nd margin
plot(density(rejsampleSMA2m2M106$L2postval, from=-1,to=1), xlim=c(-1,1), ylim=c(0,1.1), col="blue", main="theta2")
lines(density(rejsampleESMA2M106$unadj.values[,2], from=-1,to=1), xlim=c(-1,1), col="red")
lines(density(rejsampleEVSMA2M106$unadj.values[,2], from=-1,to=1), xlim=c(-1,1), col="red", lty=2)
lines(density( saabc.SMA2M106$post.sample[,2] , from=-1,to=1), xlim=c(-1,1), col="green")
lines(density(rejsampleSMA2m2M106$MW2postval, from=-1,to=1), xlim=c(-1,1), col="black")
lines(thetaseq2,postpdfSMA2$marg2, col="purple")
# true param is 0.5 here 
abline(v=0.5, lty=2)


#################################################################################################
### 2D plot /contours of the true posterior
#################################################################################################
# SMA2_df<-PostPdfSMA2(thetaseq1,thetaseq2, ytargetSMA2)$postdf
 SMA2_df<-postpdfSMA2$postdf
#
###  contours of true posterior  for y=ytargetSMA2 
v <- ggplot(SMA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.7, linetype="dashed", size=.5) + geom_hline(yintercept=0.5, linetype="dashed", size=.5)
v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()
v
#v + geom_contour(breaks=c(  0.9, 1.1, 1.2, seq(1.23,  max(z), length.out = 10)))
#v + geom_contour(bins=20)
#v + geom_density_2d(size = .5, aes(color = ..level..)) 
#

# To save a  pdf
#pdf(file = "name.pdf", width = 5.5, height = 4)
#[...]
#dev.off()

############################################################################
#
######## 2) compute distances between full post: 
# With Compute1DistISo and foreach : about 10 hours  for M=10^5!!!!
# With PreCompDistIsoPar and foreach: about 2,5 hours (with 8 cores)
#
### Possibility to reduce the mixtures before computing distances to save time
#Wthresh<-0.001

# First Speed test with M=10^4
system.time(distMW2L2SMA2M104v1<- foreach (i=1:M, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimSMA2,t(ytargetSMA2),matrix(ysimu2SMA2[,i]),wthr = 0.001)
})
# user  system elapsed 
#26863.739   151.663  3512.463 -->58 min

# first 10^4 distances
save(list = c("distMW2L2SMA2M104v1"), file = "distSMA2M104v1K80.RData")

#M<-10^5
system.time(distMW2L2SMA2<- foreach (i=1:M, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimSMA2,t(ytargetSMA2),matrix(ysimu2SMA2[,i]),wthr = 0.001)
})
#user     system    elapsed 
#268552.615    861.171  34733.420 --> 9.64 hours

# all 10^5 distances SAVED
save(list = c("distMW2L2SMA2"), file = "distSMA2K80.RData")

#Acceleration TEST : first attempt 
# test avec M<-10^5
# Comparison without dopar and accelerated  computation with Pre computation
system.time(distMW2L2SMA2test<-
              ComputeDistIsoG(modgllimSMA2,t(ytargetSMA2),ysimu2SMA2[,1:M],Wthresh)
)
#user   system  elapsed 
#33787.17   801.56 34712.92 --> 9.64 SAME TIME but user time smaller 
# but the cores are not used

# Second test with other function using cores: FASTER
system.time(distMW2L2SMA2testbis<-PreCompDistIsoPar(modgllimSMA2,t(ytargetSMA2),ysimu2SMA2[,1:1000],Wthresh)
             )
#user  system elapsed 
#613.461   9.949  90.898 instead of 300s 
# M=10000
# 5995.733   93.222  932.139 ---> 15 min instead of 58--> should be 2.5 h for 10^5
# 

###########################################################################
# Plot of Marginal posteriors  but with samples selected from 2D distances 
##########################################################################
MW2SMA2<-distMW2L2SMA2[1,]
L2SMA2<-distMW2L2SMA2[2,]

rejsampleSMA2<-RejABCMW2L2(0.01,MW2SMA2,L2SMA2,thetasimu2SMA2)

# theta1 1st margin
plot(density(rejsampleSMA2$L2postval[1,], from=-2,to=2), xlim=c(-2,2), ylim=c(0,0.53), col="blue", main="theta1")
lines(density(rejsampleESMA2$unadj.values[,1], from=-2,to=2), xlim=c(-2,2), col="red")
lines(density(rejsampleEVSMA2$unadj.values[,1], from=-2,to=2), xlim=c(-2,2), col="red", lty=2)
lines(density( saabc.SMA2$post.sample[,1] , from=-2,to=2), xlim=c(-2,2), col="green")
lines(density(rejsampleSMA2$MW2postval[1,], from=-2,to=2), xlim=c(-2,2), col="black")
lines(thetaseq1,postpdfSMA2$marg1, col="purple")
# true param is 0.7 here 
abline(v=0.7, lty=2)
abline(v=-0.7, lty=2)

# theta2 2nd margin
plot(density(rejsampleSMA2$L2postval[2,], from=-1,to=1), xlim=c(-1,1), ylim=c(0,1.1), col="blue", main="theta2")
lines(density(rejsampleESMA2$unadj.values[,2], from=-1,to=1), xlim=c(-1,1), col="red")
lines(density(rejsampleEVSMA2$unadj.values[,2], from=-1,to=1), xlim=c(-1,1), col="red", lty=2)
lines(density( saabc.SMA2$post.sample[,2] , from=-1,to=1), xlim=c(-1,1), col="green")
lines(density(rejsampleSMA2$MW2postval[2,], from=-1,to=1), xlim=c(-1,1), col="black")
lines(thetaseq2,postpdfSMA2$marg2, col="purple")
# true param is 0.5 here 
abline(v=0.5, lty=2)

## plot with 100 selected samples from 10^4 initial samples except for E, EV, SA
# postTheta1SMA2K80fullM104.png and postTheta2SMA2K80fullM104.png

##############################################################################
# Plots of selected posterior samples in 2D  with contours, from 2D distances
##############################################################################

dfEVSMA2<-data.frame(rejsampleEVSMA2$unadj.values)
dfESMA2 <-data.frame(rejsampleESMA2$unadj.values)
dfSASMA2 <-data.frame(saabc.SMA2$post.sample)
dfL2SMA2<-data.frame(t(rejsampleSMA2$L2postval))
dfMW2SMA2<-data.frame(t(rejsampleSMA2$MW2postval))

#gllim-MW2 : points and contours 
m <- ggplot(dfMW2SMA2, aes(x = X1, y = X2)) +
  # geom_point(size=.3) + 
  xlim(-2, 2) +
  ylim(-1, 1) + theme(aspect.ratio=0.9)
m<- m +  geom_point(data=dfMW2SMA2, size=.3, color="red") 
m<- m + geom_vline(xintercept=0.7, linetype="dashed", size=.5) + geom_hline(yintercept=0.5, linetype="dashed", size=.5)
m<- m +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
m<- m + geom_density_2d(size = .5, color="black") 
m<- m + xlab("theta1") + ylab("theta2")
# to remove the grey background
m + theme_bw()


## Semi automatic ABC: points and contours
m <- ggplot(dfSASMA2, aes(x = X1, y = X2)) +
  # geom_point(size=.3) + 
  xlim(-2, 2) +
  ylim(-1, 1) + theme(aspect.ratio=0.9)
m<- m +  geom_point(data=dfSASMA2, size=.3, color="red") 
m<- m + geom_vline(xintercept=0.7, linetype="dashed", size=.5) + geom_hline(yintercept=0.5, linetype="dashed", size=.5)
m<- m +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
m<- m + geom_density_2d(size = .5, color="black") 
m<-m + xlab("theta1") + ylab("theta2")
# to remove the grey background
m + theme_bw()

## gllim-L2: points and contours 
m <- ggplot(dfL2SMA2, aes(x = X1, y = X2)) +
  # geom_point(size=.3) + 
  xlim(-2, 2) +
  ylim(-1, 1) + theme(aspect.ratio=0.9)
m<- m +  geom_point(data=dfL2SMA2, size=.3, color="red") 
m<- m + geom_vline(xintercept=0.7, linetype="dashed", size=.5) + geom_hline(yintercept=0.5, linetype="dashed", size=.5)
m<- m +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
m<- m + geom_density_2d(size = .5, color="black") 
m<-m + xlab("theta1") + ylab("theta2")
# to remove the grey background
m + theme_bw()

## gllim-E: points  and contours
m <- ggplot(dfESMA2, aes(x = X1, y = X2)) +
  # geom_point(size=.3) + 
  xlim(-2, 2) +
  ylim(-1, 1) + theme(aspect.ratio=0.9)
m<- m +  geom_point(data=dfESMA2, size=.3, color="red") 
m<- m + geom_vline(xintercept=0.7, linetype="dashed", size=.5) + geom_hline(yintercept=0.5, linetype="dashed", size=.5)
m<- m +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
m<- m + geom_density_2d(size = .5, color="black") 
m<- m + xlab("theta1") + ylab("theta2")
# to remove the grey background
m + theme_bw()

## gllim-EV: points and contours
m <- ggplot(dfEVSMA2, aes(x = X1, y = X2)) +
  # geom_point(size=.3) + 
  xlim(-2, 2) +
  ylim(-1, 1) + theme(aspect.ratio=0.9)
m<- m +  geom_point(data=dfEVSMA2, size=.3, color="red") 
m<- m + geom_vline(xintercept=0.7, linetype="dashed", size=.5) + geom_hline(yintercept=0.5, linetype="dashed", size=.5)
m<- m +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
m<- m + geom_density_2d(size = .5, color="black") 
m<- m + xlab("theta1") + ylab("theta2")
# to remove the grey background
m + theme_bw()


##############################################################################################
##############################################################################################
##
## Example ITD, D=10, L=2
##
##############################################################################################
### Step 1: Setting: target y (observations to be inverted, dimension is D) ;
###         Data sets simulation, one for GLLiM, one for ABC if wanted
# output: ytarget; (thetasimu1, ysimu1) size N; (thetasimu2, ysimu2) size M ;
#  true posterior if possible

load("dataITDdf3.RData")
# recover c("ytargetITD","ysimu1ITD","thetasimu1ITD","ysimu2ITD","thetasimu2ITD","MHvalITDT" )



### Step 2: Fitting GLLiM on (thetasim1, ysim1), model selection (K, constraints on Sigmas) 
### and computing required quantities not directly provided by xLLiM
###########################################################################################################
# xllim  K (chosen manually here) gaussian components  and  (default) isotropic constraints on the Sigma_k
# data and param should be D x N and L x N
# output: a list with the output of xllim + inverse Gammak and Sigmakstar

##### application ITD one pair, dof=3 
K=20
system.time(modgllimITD <- GllimIsoFit(thetasimu1ITD, ysimu1ITD, K))
#  user  system elapsed 
#146.247  53.654 200.478 --> 3 min
save(list = c("modgllimITD"), file = "modgllimITDK20df3simu1.RData")
modgllimITDsimu1<-modgllimITD

# with the larger data set
system.time(modgllimITD <- GllimIsoFit(thetasimu2ITD, ysimu2ITD, K))
#  user  system elapsed 
#1254.550  175.702 1433.473 --> 23 min

save(list = c("modgllimITD"), file = "modgllimITDK20df3simu2.RData")

# one pair in diago df=300
K=20
system.time(modgllimITD <- GllimIsoFit(thetasimu1ITD, ysimu1ITD, K))
#  user  system elapsed 
#146.247  53.654 200.478 --> 3 min
save(list = c("modgllimITD"), file = "modgllimITDK20df300diagsimu1.RData")

# check with simu1 diag with df =1
system.time(modgllimITD <- GllimIsoFit(thetasimu1ITD, ysimu1ITD, K))
#  user  system elapsed 
#146.247  53.654 200.478 --> 3 min
modgllimITDsimu1<-modgllimITD
save(list = c("modgllimITD"), file = "modgllimITDK20df1diagsimu1.RData")

# simu2 diag df =1 
system.time(modgllimITD <- GllimIsoFit(thetasimu2ITD, ysimu2ITD, K))
#  user  system elapsed 
# ici

### other tests and setting
##### New config (-1,0) (1,0) and gllim on same data set as SA
K=20
system.time(modgllimITDs2 <- GllimIsoFit(thetasimu2ITD, ysimu2ITD, K))
#  user  system elapsed 
# 1659.138  187.076 1854.710 -->31 min
save(list = c("modgllimITDs2"), file = "modgllimITDK20s2.RData")

##### New config (-0.5,0) (0.5,0) dofT=300 and gllim on same data set as SA
##### ---> Gaussian noise
K=20
system.time(modgllimITDg <- GllimIsoFit(thetasimu2ITD, ysimu2ITD, K))
#  user  system elapsed 
#2402.278  284.726 2699.027 -->45 min

save(list = c("modgllimITDg"), file = "modgllimITDK20g.RData")

#### For the mixture model
K=20
system.time(modgllimITDgmix <- GllimIsoFit(thetasimu2ITD, ysimu2ITD, K))
#  user  system elapsed 
#1852.820  230.204 2084.032 --> 34 min

# dof=300 sigma=0.01
save(list = c("modgllimITDgmix"), file = "modgllimITDK20gmix.RData")

######## dof=3 sigma=0.01 ########
save(list = c("modgllimITDgmix"), file = "modgllimITDK20gmixdf3.RData")

### Gllim learn with the smaller dataset simu1, dof=3 sigma=0.01 ########
system.time(modgllimITDgmix1 <- GllimIsoFit(thetasimu1ITD, ysimu1ITD, K))
#  user  system elapsed 
#122.031  47.955 171.495 
save(list = c("modgllimITDgmix1"), file = "modgllimITDK20gmixdf3simu1.RData")


##################################################################################
###
### Step 3.1: GLLiM-E and EV using package abc for the ABC scheme 
###
#################################################################################################
# ABC avec sum stat from gllim: 1) espectation
##################################
# summary stat = posterior expectations computed from gllim for the training examples
# ysimu2 is D x M 
# modg is the output of GllimIsoFit
# M1=10^6 : espysimu2 is of size M1 x L 
system.time(espysimu2ITD<-t(gllim_inverse_map(ysimu2ITD,modgllimITD$mod)$x_exp))
#user  system elapsed 
#14.028   1.967  16.109 

##### Plot to visualize the posterior means estimated by Gllim: in the ITD example 
# they should all be closed to 0. Don't plot all the N or M points (too long)
# to plot 2 mics and source locations
dfconfig<-data.frame(matrix(c(-0.5,-0.5,0.5,0.5,0,1), 3,2, byrow=T))

dfExpGllim<-data.frame(espysimu2ITD[1:1000,])
m <- ggplot(dfExpGllim, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfExpGllim, size=1, color=rgb(.1,0,.9,alpha = 0.5), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m

# to check all N or M plot histograms:
hist(espysimu2ITD[,1], nclass=20, xlim=c(-2,2), xlab = "x", main="x")
hist(espysimu2ITD[,2], nclass=50, xlim=c(-2,2), xlab = "y", main="y")

# sum stat = posterior espectation for the target (obs or y test)
# ytarget is 1 x D
# espytarget is 1 x L 
espytargetITD=t(gllim_inverse_map(t(ytargetITD),modgllimITD$mod)$x_exp)

#ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for wasserstein ...
# tol =0.001 select 1000 samples when M1 =10^6]

# thetasimu2 is L x M1
# selected sample of size 1000 x 2 is in rejsampleEITD$unadj.values
# Thrshold 0.001 for 1000 samples
rejsampleEITD <-abc(target=espytargetITD, param=t(thetasimu2ITD), sumstat=espysimu2ITD, tol=.001, method ="rejection")

# Check plot
# plot(rejsampleEITD$unadj.values, xlim=c(-2,2),ylim=c(-2,2), pch=20, cex=.5)
# 

##### Plot to visualize the K mixture estimated by Gllim for ytargetITD
# With package MixSim and function simdataset
# requires  the gllim mixture parameters Pi, Mu, S for ytargetITD

mixparam<-ComputeGllimMixtParam(modgllimITD,t(ytargetITD),0)
samplemixgllim<-simdataset(1000,mixparam$Pi, t(mixparam$Mu),mixparam$S)$X
dfmix<-data.frame(samplemixgllim)

m <- ggplot(dfmix, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfmix, size=1, color=rgb(1,0,0,alpha = 1), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m


##################################
# ABC avec sum stat from gllim: 2) espectation + log variances
##################################

logvarytargetITD<-LogPostVarGllimIso(modgllimITD,t(ytargetITD))

system.time(logvarysimu2ITD<-LogPostVarGllimIso(modgllimITD, ysimu2ITD))
#user  system elapsed 
# 455.056   4.889 461.437 

##### Plot to visualize the 2 posterior log variances estimated by Gllim: in the ITD example 
# Don't plot all the N or M points (too long)

dflogVGllim<-data.frame(t(logvarysimu2ITD)[1:10000,])
m <- ggplot(dflogVGllim, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dflogVGllim, size=1, color=rgb(.1,0,.9,alpha = 0.5), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m

# to check all N or M plot histograms:
hist(logvarysimu2ITD[1,], nclass=20, xlim=c(-2,2), xlab = "var", main="x posterior log variance")
hist(logvarysimu2ITD[2,], nclass=50, xlim=c(-2,2), xlab = "var", main="y posterior log variance")


####### Gllim-EV
#ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for wasserstein ...
# tol =0.001 select 1000 samples when M1 =10^6]

# thetasimu2 is L x M
# selected sample of size L x 1000 is in rejsampleEV$unadj.values
# esplogytarget  1 x 2L 
# esplogysimu2   M x 2L
esplogytargetITD<-c(espytargetITD, t(logvarytargetITD))
esplogysimu2ITD<-cbind(espysimu2ITD, t(logvarysimu2ITD))

# threshold 0.001---> 1000 samples
rejsampleEVITD <-abc(target=esplogytargetITD, param=t(thetasimu2ITD), sumstat=esplogysimu2ITD, tol=.001, method ="rejection")


#plot(rejsampleEVITD$unadj.values, xlim=c(-2,2),ylim=c(-2,2), pch=20, cex=.5)

## SEMI AUTO ABC-- TODO a fair comparison in the usage of data 1 and 2 
mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
# with M1=10^6 and epsilon =0.001
saabc.ITD<-semiauto.abc(obs=ytargetITD,param=t(thetasimu2ITD), sumstats=t(ysimu2ITD),satr=mytf, tol=0.001, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)
## To get plots of parameters vs fitted values, add plot=TRUE but this slow down
## a lot when used with 10^6 data points

# For comparison similar plots with GlliM:
plot(thetasimu2ITD[1,1:1000], espysimu2ITD[1:1000,1],xlim=c(-2,2), ylim=c(-2,2), xlab="True value", ylab="Posterior mean", main="x")
plot(thetasimu2ITD[2,1:1000], espysimu2ITD[1:1000,2],xlim=c(-2,2), ylim=c(-2,2), xlab="True value", ylab="Posterior mean", main="y")


#################################################################################################
###
### Step 3.2: computing distances, L2, MW2 
###
#################################################################################################
# Cases : L=1 or L>1 (to be unified later?) ; L2, MW2 

###########################################################################
##   ITD is a Case L=2 : 1) compute post margins , 2) compute full post     
##########################################################################

### Possibility to reduce the mixtures before computing distances to save time
Wthresh<-0.001

######### 1) Compute posterior margins is much FASTER (NOT DONE for ITD)
## ComputeDistMargin: Compute all z without foreach faster? similar in speed
## Not faster
#system.time(distMW2L2SMA2marg12test<-
#                ComputeDistMarginIso(modgllimSMA2,t(ytargetSMA2),ysimu2SMA2,wthr = 0.001))
##user  system elapsed 
##Timing stopped at: 2786 355.4 3150---> 52 min ???

# To compare in speed with 
numCores <- detectCores() #numCores = 8
registerDoParallel(numCores)

M1<-10^6    ## NOT DONE
system.time(distMW2L2ITDmarg12<- foreach (i=1:M1, .combine=cbind) %dopar% {
  Compute1DistMarginIso(modgllimITD,t(ytargetITD),matrix(ysimu2ITD[,i]),wthr = 0.001)
})
#user  system elapsed 
#15663.188   266.815  2158.846 ---> 36 min

save(list = c("distMW2L2ITDmarg12"), file = "distITDmarg12K20.RData")
### end NOT DONE

############################################################################
#
######## 2) compute distances between full post: 
# avec Compute1DistIso about 4 hours  for M=10^6
# USE PreCompDistIsoPar to accelerate
#
### Possibility to reduce the mixtures before computing distances to save time
#Wthresh<-0.001
# First test with M=10^3 to assess the time
system.time(distMW2L2ITDM103v1<- foreach (i=1:1000, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimITD,t(ytargetITD),matrix(ysimu2ITD[,i]),wthr = 0.001)
})
# user  system elapsed 
# 94.479   7.884  14.946 

#M1<-10^6 dof=3 simu2
system.time(distMW2L2ITD<- foreach (i=1:M1, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimITD,t(ytargetITD),matrix(ysimu2ITD[,i]),wthr = 0.001)
})
#user     system    elapsed 
# 129223.059   2018.199  17831.911 --> 5h

# save all 10^6 distances
save(list = c("distMW2L2ITD"), file = "distITDK20df3simu2.RData")

# dof=3, sigma = 0.01, mixture (gmix)
#M1<-10^6
system.time(distMW2L2ITD<- foreach (i=1:M1, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimITDgmix,t(ytargetITD),matrix(ysimu2ITD[,i]),wthr = 0.001)
})
#user     system    elapsed 
#325764.429   2489.098  42751.873--> 12h !!

# save all 10^6 distances
save(list = c("distMW2L2ITD"), file = "distITDK20mixdf3.RData")

#### TO DO with simu1 and modgllimITDgmix1
# Temp save 
save(list = ls(), file = "env_entierITDmixdf3.RData")

system.time(distMW2L2ITD<- foreach (i=1:M1, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimITDgmix1,t(ytargetITD),matrix(ysimu2ITD[,i]),wthr = 0.001)
})
#user    system   elapsed 
#292902.59  88576.58  65103.85 ---> 18h !!!!

# save all 10^6 distances
save(list = c("distMW2L2ITD"), file = "distITDK20mixdf3simu1.RData")


### Step 4: ABC scheme, default rejection ABC using distance quantiles
#################################################################################################
MW2ITD<-distMW2L2ITD[1,]
L2ITD<-distMW2L2ITD[2,]

rejsampleITD<-RejABCMW2L2(0.001,MW2ITD,L2ITD,thetasimu2ITD)

##############################################################################
# Plots of selected posterior samples in 2D, from 2D distances
##############################################################################

dfEVITD<-data.frame(rejsampleEVITD$unadj.values)
dfEITD <-data.frame(rejsampleEITD$unadj.values)
dfSAITD <-data.frame(saabc.ITD$post.sample)
dfL2ITD<-data.frame(t(rejsampleITD$L2postval))
dfMW2ITD<-data.frame(t(rejsampleITD$MW2postval))


##################################
##### PLOT true post contours with microphones
####################################

# microphones positions and true source location
#dfconfig<-data.frame(matrix(c(-1,0,1,1,0.6,1), 3,2, byrow=T))
# dfconfig<-data.frame(matrix(c(-1,0,1,0,0.6,1), 3,2, byrow=T))
dfconfig<-data.frame(matrix(c(-0.5,0,0.5,0,0,1), 3,2, byrow=T))

###  True posterior unnormalized contours, PostPdfITD is in ITD.R
N_grid=500
ITD_df<-PostPdfITD(N_grid,ytargetITD)$postdf

v <- ggplot(ITD_df) +  xlab("x") + ylab("y")+ theme(aspect.ratio = 1)+ xlim(-2, 2) + ylim(-2, 2) ; 
v + geom_contour( aes(x1, x2, z = z), color="blue", bins=7) + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3) + geom_point(data=dfconfig, mapping=aes(x=X1, y=X2), size=2.5, color="black") 

### Metropolis Hastings sample plot
dfMHITD<-data.frame(MHvalITD)

m <- ggplot(dfMHITD, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfMHITD, size=1, color=rgb(0,0,1,alpha = 1), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m

####################################
####
#### Plots with colors
####
####################################
#a) rgb(.1,0,.9,alpha = 0.5)
#b) rgb(.5,0,.5,alpha = 0.5)
#c) rgb(.9,0,.1,alpha = 0.5)

####### GLLIM-E-ABC  
m <- ggplot(dfEITD, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfEITD, size=1, color=rgb(.9,0,.1,alpha = 0.5), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m

####### GLLIM-EV-ABC  
m <- ggplot(dfEVITD, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfEVITD, size=1, color=rgb(.9,0,.1,alpha = 0.5), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m

####### Semi Automatic ABC  
m <- ggplot(dfSAITD, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfSAITD, size=1, color=rgb(0,0.5,0.5,alpha = 1), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m

####### GLLIM-MW2-ABC 
m <- ggplot(dfMW2ITD, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfMW2ITD, size=1, color=rgb(0,0,.1,alpha = 0.5), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m

####### GLLIM-L2-ABC  
m <- ggplot(dfL2ITD, aes(x = X1, y = X2)) +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfL2ITD, size=1, color=rgb(0.1,0,.9,alpha = 0.5), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m


############## PLOTs to  visualise all computed distances and to compare L2 and MW2
# Coloring with respect to the distance value, for a no-white plot, 10^5 is enough
# REMARK: the plot generation can be so slow that it does not seem to work but it does! 

dfprior<-data.frame(t(thetasimu2ITD[,1:100000]))

# 10^5 L2 or 10^4 there are white holes
dist<-L2ITD[1:100000]
sp2<-ggplot(dfprior, aes(x=X1, y=X2, color=dist)) +
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
sp2<- sp2 + geom_point()
sp2+scale_color_gradient(low="blue", high="red")

# Same for MW2:  less contrasted here
dfprior<-data.frame(t(thetasimu2ITD[,1:100000]))

# 10^5 L2 or 10^4 there are white holes
dist<-MW2ITD[1:100000]
sp2<-ggplot(dfprior, aes(x=X1, y=X2, color=dist)) +
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
sp2<- sp2 + geom_point()
sp2+scale_color_gradient(low="blue", high="red")



##### END Plots for DRAFT #############
# To save a  pdf
#pdf(file = "name.pdf", width = 5.5, height = 4)
#[...]
#dev.off()
####################################################

# Temp save first ITD case (paper)
save(list = ls(), file = "env_entierITDv1.RData")

# save df3 simu2 test1 (1.5 1)
save(list = ls(), file = "env_entierITDdf3test1.RData")

#################################################################################################

##############################################################################################
##
## Real data Example Hapke model inversion  D=10, L=4
## Data sets are simulated with a Gaussian noise with sigma =0.005 (see paper)
##
##############################################################################################
### Step 1: Setting: 
# output: ytarget; (thetasimu1, ysimu1) size N; (thetasimu2, ysimu2) size M ;
# here N = M =10^5 

load("dataHapke.RData")
# recover c("ytargetreal", "ytargetH005","ysimu1H005","thetasimu1H005","ysimu2H005","thetasimu2H005"), file = "dataHapke.RData")

### Step 2: Fitting GLLiM on (thetasim1, ysim1), model selection (K, constraints on Sigmas) 
### and computing required quantities not directly provided by xLLiM
###########################################################################################################
# xllim  K (chosen manually here) gaussian components  and  (default) isotropic constraints on the Sigma_k
# data and param should be D x N and L x N
# output: a list with the output of xllim + inverse Gammak and Sigmakstar

##### application Hapke + noise 0.005
K=40
system.time(modgllimH005 <- GllimIsoFit(thetasimu1H005, ysimu1H005, K))
#user   system  elapsed 
#1179.279  106.380 1289.556 --> 21 min

save(list = c("modgllimH005"), file = "modgllimH005K40.RData")
##################################################################################

### Step 3.1: GLLiM-E and EV using package abc for the ABC scheme 
###
#################################################################################################
# ABC avec sum stat from gllim: 1) espectation
##################################
# summary stat = posterior expectations computed from gllim for the training examples
# ysimu2 is D x M 
# modg is the output of GllimIsoFit
# espysimu2 is of size M x L 
espysimu2H005<-t(gllim_inverse_map(ysimu2H005,modgllimH005$mod)$x_exp)

# sum stat = posterior espectation for the target (obs or y test)
# ytarget is 1 x D
# espytarget is 1 x L 
#
# 1) inversion of the simulated observation
espytargetH005=t(gllim_inverse_map(t(ytargetH005),modgllimH005$mod)$x_exp)
# 2) inversion od the real observation
espytargetreal=t(gllim_inverse_map(t(ytargetreal),modgllimH005$mod)$x_exp)


#ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for wasserstein ...
# tol =0.001 select 100 samples when Nsimu =10^5]

# thetasimu2 is L x M
# selected sample of size L x 100 is in rejsampleESMA1$unadj.values
# Thrshold 0.01 for 1000 samples
#
# simulated observation
rejsampleEH005 <-abc(target=espytargetH005, param=t(thetasimu2H005), sumstat=espysimu2H005, tol=.01, method ="rejection")

# real observation 0.1% quantile, 0.001 and 100 selected
rejsampleEHreal <-abc(target=espytargetreal, param=t(thetasimu2H005), sumstat=espysimu2H005, tol=.001, method ="rejection")


##################################
# ABC avec sum stat from gllim: 2) espectation + log variances
##################################

# simulated observation
logvarytargetH005<-LogPostVarGllimIso(modgllimH005,t(ytargetH005))
# real observation
logvarytargetreal<-LogPostVarGllimIso(modgllimH005,t(ytargetreal))


system.time(logvarysimu2H005<-LogPostVarGllimIso(modgllimH005, ysimu2H005))
#user  system elapsed 
#108.805   3.405 112.920 


####### Gllim-EV
#ABC scheme from abc package, rejection ABC, other are possible using "ridge" or "loclinear" but no big differences observed and 
# rejection is the one used for wasserstein ...
# tol =0.001 select 100 samples when Nsimu =10^5]

# thetasimu2 is L x M
# selected sample of size L x 100 is in rejsampleEV$unadj.values
# esplogytarget  1 x 2L 
# esplogysimu2   M x 2L
esplogytargetH005<-c(espytargetH005, t(logvarytargetH005))
esplogytargetreal<-c(espytargetreal, t(logvarytargetreal))

esplogysimu2H005<-cbind(espysimu2H005, t(logvarysimu2H005))

# threshold 0.01---> 1000 samples
rejsampleEVH005 <-abc(target=esplogytargetH005, param=t(thetasimu2H005), sumstat=esplogysimu2H005, tol=.01, method ="rejection")

# for real obs, 0.1% quantile, 0.001 and 100 selected
rejsampleEVHreal <-abc(target=esplogytargetreal, param=t(thetasimu2H005), sumstat=esplogysimu2H005, tol=.001, method ="rejection")

## SEMI AUTO ABC-- TODO a fair comparison in the usage of data 1 and 2 
mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
# DONT Forget to change: for SMA2 threshold is 0.01
saabc.H005<-semiauto.abc(obs=ytargetH005,param=t(thetasimu2H005), sumstats=t(ysimu2H005),satr=mytf, tol=0.01, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)

# for real obs, 0.1% quantile, 0.001 and 100 selected
saabc.Hreal<-semiauto.abc(obs=ytargetreal,param=t(thetasimu2H005), sumstats=t(ysimu2H005),satr=mytf, tol=0.001, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)

#################################################################################################
###
### Step 3.2: computing distances, L2, MW2 
###
#################################################################################################
###########################################################################
##   Hapke is a Case L=4 : 1) compute post margins , 2) compute full post     
##########################################################################

######### 1) Compute posterior margins is much FASTER BUT RESULTS are quite
### different in this example from the full posterior distances.
#### see 2) below

### Possibility to reduce the mixtures before computing distances to save time
Wthresh<-0.001


## ComputeDistMargin: Compute all z without foreach faster? similar in speed
## Not faster
#system.time(distMW2L2SMA2marg12test<-
#                ComputeDistMarginIso(modgllimSMA2,t(ytargetSMA2),ysimu2SMA2,wthr = 0.001))
##user  system elapsed 
##Timing stopped at: 2786 355.4 3150---> 52 min ???

# To compare in speed with 
numCores <- detectCores() #numCores = 8
registerDoParallel(numCores)

#########################
# Simulated observation
#########################
M<-10^5
system.time(distMW2L2H005marg12<- foreach (i=1:M, .combine=cbind) %dopar% {
  Compute1DistMarginIso(modgllimH005,t(ytargetH005),matrix(ysimu2H005[,i]),wthr = 0.001)
})
#user   system  elapsed 
#7464.092 3259.898 2152.568 --->36 min

save(list = c("distMW2L2H005marg12"), file = "distH005marg12K40.RData")

#########################
# Real observation
#########################
M<-10^5
system.time(distMW2L2Hrealmarg12<- foreach (i=1:M, .combine=cbind) %dopar% {
  Compute1DistMarginIso(modgllimH005,t(ytargetreal),matrix(ysimu2H005[,i]),wthr = 0.001)
})
#user   system  elapsed 
#7305.034 1702.508 1643.209 --> 27 min

save(list = c("distMW2L2Hrealmarg12"), file = "distHrealmarg12K40.RData")


### Step 4: ABC scheme, default rejection ABC using distance quantiles
#################################################################################################

########################
# Simulated observation
########################
MW2H005m1<-distMW2L2H005marg12[1,]
MW2H005m2<-distMW2L2H005marg12[2,]
MW2H005m3<-distMW2L2H005marg12[3,]
MW2H005m4<-distMW2L2H005marg12[4,]
L2H005m1<-distMW2L2H005marg12[5,]
L2H005m2<-distMW2L2H005marg12[6,]
L2H005m3<-distMW2L2H005marg12[7,]
L2H005m4<-distMW2L2H005marg12[8,]


# varying quantile 1% = 0.01 (1000), 0.1% = 0.001 (100), 0.05% = 0.0005 (50)
# 1st margin
rejsampleH005m1<-RejABCMW2L2(0.01,MW2H005m1,L2H005m1,t(thetasimu2H005[1,]))
# 2nd margin
rejsampleH005m2<-RejABCMW2L2(0.01,MW2H005m2,L2H005m2,t(thetasimu2H005[2,]))
# 3rd margin
rejsampleH005m3<-RejABCMW2L2(0.01,MW2H005m3,L2H005m3,t(thetasimu2H005[3,]))
# 4th margin
rejsampleH005m4<-RejABCMW2L2(0.01,MW2H005m4,L2H005m4,t(thetasimu2H005[4,]))

# return(list("MW2postval"= MW2postval,"L2postval"= L2postval ))


########################
# Real observation 0.1% quantile
########################
MW2Hrealm1<-distMW2L2Hrealmarg12[1,]
MW2Hrealm2<-distMW2L2Hrealmarg12[2,]
MW2Hrealm3<-distMW2L2Hrealmarg12[3,]
MW2Hrealm4<-distMW2L2Hrealmarg12[4,]
L2Hrealm1<-distMW2L2Hrealmarg12[5,]
L2Hrealm2<-distMW2L2Hrealmarg12[6,]
L2Hrealm3<-distMW2L2Hrealmarg12[7,]
L2Hrealm4<-distMW2L2Hrealmarg12[8,]


# 1st margin
rejsampleHrealm1<-RejABCMW2L2(0.001,MW2Hrealm1,L2Hrealm1,t(thetasimu2H005[1,]))
# 2nd margin
rejsampleHrealm2<-RejABCMW2L2(0.001,MW2Hrealm2,L2Hrealm2,t(thetasimu2H005[2,]))
# 3rd margin
rejsampleHrealm3<-RejABCMW2L2(0.001,MW2Hrealm3,L2Hrealm3,t(thetasimu2H005[3,]))
# 4th margin
rejsampleHrealm4<-RejABCMW2L2(0.001,MW2Hrealm4,L2Hrealm4,t(thetasimu2H005[4,]))



#### Plots of posterior margins from distances on surrogate posterior margins
#################################################################################################

######################################################
# Simulated obseration inversion
# True param is thetasimu2H005[,57702] = (0.68, 0.04, 0.23, 0.04)
# w 1st margin
plot(density(rejsampleH005m1$L2postval, from=0,to=1), xlim=c(0,1), col="blue", main="w")
lines(density(rejsampleEH005$unadj.values[,1], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVH005$unadj.values[,1], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.H005$post.sample[,1] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleH005m1$MW2postval, from=0,to=1), xlim=c(0,1), col="black")
# true param is 0.68 here 
abline(v=0.68, lty=2)

# theta 2nd margin
plot(density(rejsampleH005m2$L2postval, from=0,to=1), xlim=c(0,1),  col="blue", main="theta")
lines(density(rejsampleEH005$unadj.values[,2], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVH005$unadj.values[,2], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.H005$post.sample[,2] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleH005m2$MW2postval, from=0,to=1), xlim=c(0,1), col="black")
# true param is 0.04 here 
abline(v=0.04, lty=2)

#  b 3rd margin
plot(density(rejsampleH005m3$L2postval, from=0,to=1), xlim=c(0,1), col="blue", main="b")
lines(density(rejsampleEH005$unadj.values[,3], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVH005$unadj.values[,3], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.H005$post.sample[,3] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleH005m3$MW2postval, from=0,to=1), xlim=c(0,1), col="black")
# true param is 0.68 here 
abline(v=0.23, lty=2)

# c 4th margin
plot(density(rejsampleH005m4$L2postval, from=0,to=1), xlim=c(0,1),  col="blue", main="c")
lines(density(rejsampleEH005$unadj.values[,4], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVH005$unadj.values[,4], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.H005$post.sample[,4] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleH005m4$MW2postval, from=0,to=1), xlim=c(0,1), col="black")
# true param is 0.04 here 
abline(v=0.04, lty=2)

######################################################
# Real obseration inversion: True param is not known
# w 1st margin
plot(density(rejsampleHrealm1$L2postval, from=0,to=1), xlim=c(0,1),ylim=c(0,12), col="blue", main="w")
lines(density(rejsampleEHreal$unadj.values[,1], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVHreal$unadj.values[,1], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.Hreal$post.sample[,1] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleHrealm1$MW2postval, from=0,to=1), xlim=c(0,1), col="black")
# param estim
abline(v=0.59, lty=2)

# theta 2nd margin
plot(density(rejsampleHrealm2$L2postval, from=0,to=1), xlim=c(0,1),  col="blue", main="theta")
lines(density(rejsampleEHreal$unadj.values[,2], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVHreal$unadj.values[,2], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.Hreal$post.sample[,2] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleHrealm2$MW2postval, from=0,to=1), xlim=c(0,1), col="black")
#param is 0.42 or 0.15 here 
abline(v=0.42, lty=2)
abline(v=0.15, lty=2)

#  b 3rd margin
plot(density(rejsampleHrealm3$L2postval, from=0,to=1), xlim=c(0,1), col="blue", main="b")
lines(density(rejsampleEHreal$unadj.values[,3], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVHreal$unadj.values[,3], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.Hreal$post.sample[,3] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleHrealm3$MW2postval, from=0,to=1), xlim=c(0,1), col="black")
# param is 0.14 here 
abline(v=0.14, lty=2)

# c 4th margin
plot(density(rejsampleHrealm4$L2postval, from=0,to=1), xlim=c(0,1),  col="blue", main="c")
lines(density(rejsampleEHreal$unadj.values[,4], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVHreal$unadj.values[,4], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.Hreal$post.sample[,4] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleHrealm4$MW2postval, from=0,to=1), xlim=c(0,1), col="black")
# true param is 0.06 here 
abline(v=0.06, lty=2)

##############################################################################
##
## 2) compute distances between full post: for Hapke the post are different ??
## Case in the paper 
##
#Wthresh<-0.001
# USE PreCompDistIsoPar to accelerate

########### Simulated observation
# First test with M=1000
system.time(distMW2L2H005<- foreach (i=1:1000, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimH005,t(ytargetH005),matrix(ysimu2H005[,i]),wthr = 0.001)
})
#user  system elapsed 
#263.472  55.430  58.256 

#M<-10^5
system.time(distMW2L2H005<- foreach (i=1:M, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimH005,t(ytargetH005),matrix(ysimu2H005[,i]),wthr = 0.001)
})
#user     system    elapsed 
#24011.502   884.821  3803.225 ---> 1h

# all 10^5 distances
save(list = c("distMW2L2H005"), file = "distH005fullK40.RData")

########### Real observation

# First test with M=1000
system.time(distMW2L2Hreal<- foreach (i=1:1000, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimH005,t(ytargetreal),matrix(ysimu2H005[,i]),wthr = 0.001)
})
#user  system elapsed 
#320.486  23.338  54.277 

#M<-10^5
system.time(distMW2L2Hreal<- foreach (i=1:M, .combine=cbind) %dopar% {
  Compute1DistIso(modgllimH005,t(ytargetreal),matrix(ysimu2H005[,i]),wthr = 0.001)
})
#user     system    elapsed 
#31122.799   734.598  4820.755 --> 1.3 h 

# all 10^5 distances
save(list = c("distMW2L2Hreal"), file = "distHrealfullK40.RData")

###########################################################################
# Plot of Marginal posteriors  but with samples selected from 2D distances 
##########################################################################

### Simulated  observation  with GT

MW2H005<-distMW2L2H005[1,]
L2H005<-distMW2L2H005[2,]

### several quantile 0.01, 0.001, 0.0005
rejsampleH005Full<-RejABCMW2L2(0.01,MW2H005,L2H005,thetasimu2H005)

# w 1st margin
plot(density(rejsampleH005Full$L2postval[1,], from=0,to=1), xlim=c(0,1), ylim=c(0,4), col="blue", main="w")
lines(density(rejsampleEH005$unadj.values[,1], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVH005$unadj.values[,1], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.H005$post.sample[,1] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleH005Full$MW2postval[1,], from=0,to=1), xlim=c(0,1), col="black")
# 
abline(v=0.68, lty=2)

# theta 2nd margin
plot(density(rejsampleH005Full$L2postval[2,], from=0,to=1), xlim=c(0,1), ylim=c(0,2.3), col="blue", main="theta")
lines(density(rejsampleEH005$unadj.values[,2], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVH005$unadj.values[,2], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.H005$post.sample[,2] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleH005Full$MW2postval[2,], from=0,to=1), xlim=c(0,1), col="black")
#
abline(v=0.04, lty=2)


# b 3rd  margin
plot(density(rejsampleH005Full$L2postval[3,], from=0,to=1), xlim=c(0,1), ylim=c(0,2.3), col="blue", main="b")
lines(density(rejsampleEH005$unadj.values[,3], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVH005$unadj.values[,3], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.H005$post.sample[,3] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleH005Full$MW2postval[3,], from=0,to=1), xlim=c(0,1), col="black")
# 
abline(v=0.23, lty=2)

# c 4th margin
plot(density(rejsampleH005Full$L2postval[4,], from=0,to=1), xlim=c(0,1), ylim=c(0,3.5), col="blue", main="c")
lines(density(rejsampleEH005$unadj.values[,4], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVH005$unadj.values[,4], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.H005$post.sample[,4] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleH005Full$MW2postval[4,], from=0,to=1), xlim=c(0,1), col="black")
#
abline(v=0.04, lty=2)

##################################
### Real observation 

MW2Hreal<-distMW2L2Hreal[1,]
L2Hreal<-distMW2L2Hreal[2,]

rejsampleHrealFull<-RejABCMW2L2(0.001,MW2Hreal,L2Hreal,thetasimu2H005)

# w 1st margin
plot(density(rejsampleHrealFull$L2postval[1,], from=0,to=1), xlim=c(0,1), ylim=c(0,9), col="blue", main="w")
lines(density(rejsampleEHreal$unadj.values[,1], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVHreal$unadj.values[,1], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.Hreal$post.sample[,1] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleHrealFull$MW2postval[1,], from=0,to=1), xlim=c(0,1), col="black")
# 
abline(v=0.59, lty=2)

# theta 2nd margin
plot(density(rejsampleHrealFull$L2postval[2,], from=0,to=1), xlim=c(0,1), ylim=c(0,2.3), col="blue", main="theta")
lines(density(rejsampleEHreal$unadj.values[,2], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVHreal$unadj.values[,2], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.Hreal$post.sample[,2] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleHrealFull$MW2postval[2,], from=0,to=1), xlim=c(0,1), col="black")
#
abline(v=0.15, lty=2)
abline(v=0.42, lty=2)


# b 3rd  margin
plot(density(rejsampleHrealFull$L2postval[3,], from=0,to=1), xlim=c(0,1), ylim=c(0,2.3), col="blue", main="b")
lines(density(rejsampleEHreal$unadj.values[,3], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVHreal$unadj.values[,3], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.Hreal$post.sample[,3] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleHrealFull$MW2postval[3,], from=0,to=1), xlim=c(0,1), col="black")
# 
abline(v=0.14, lty=2)

# c 4th margin
plot(density(rejsampleHrealFull$L2postval[4,], from=0,to=1), xlim=c(0,1), ylim=c(0,3.5), col="blue", main="c")
lines(density(rejsampleEHreal$unadj.values[,4], from=0,to=1), xlim=c(0,1), col="red")
lines(density(rejsampleEVHreal$unadj.values[,4], from=0,to=1), xlim=c(0,1), col="red", lty=2)
lines(density( saabc.Hreal$post.sample[,4] , from=0,to=1), xlim=c(0,1), col="green")
lines(density(rejsampleHrealFull$MW2postval[4,], from=0,to=1), xlim=c(0,1), col="black")
#
abline(v=0.06, lty=2)

##### Reconstruction plots to check estimation, reconstruction made with Julia code
##### for Hapke in Logiciels/Julia2020

## Spurious b-mode for theta= [0.59, 0.42, 0.5, 0.006]
##Signal is:  0.284254 0.279314 0.247754 0.229027 0.221049 0.223845 0.237490 0.265420
##  0.313732 0.359426
## Two other solutions are 
## Theta= [0.59, 0.42, 0.14, 0.006]
## 0.343595 0.323993 0.284616 0.276913 0.277209 0.286020 0.304590 0.335094
## 0.380764 0.421559
## Theta= [0.59, 0.15, 0.14, 0.006]
## 0.338851 0.315779 0.275857 0.268493 0.268899 0.277580 0.295785 0.325652
## 0.365386 0.388834

## plot with ts.plot(cbind(signal1, signal2, etc))
