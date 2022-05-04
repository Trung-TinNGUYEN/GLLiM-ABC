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

# contain obs the observations D x R
#obsfile <- paste0(prefix, "ma2data.d", d, ".n", nobservations, ".RData")
#load(obsfile)
# here:
obs<-matrix(t(ytargetMA2150), ncol=5)

# function to simulate data
target$simulate <- function(theta){
  return(target$robservation(nobservations, theta))
}

# euclidean distance
# and also Euclidean distance
ground_p <- 2
deuclidean <- function(z){
  return(mean(apply(abs(z - obs), 2, function(v) (sum(v^ground_p))^(1/ground_p))^p)^(1/p))
}

# wasserstein distance
wdistance <- get_transport_to_y(obs, p = p)
#



########## prep for distance W2 and L2
# obs = observcations is 30 x 5

## Gllim model used
#modg<-modgllimMA2R5IIDFullK20 
# for TableS2:
modg<-modgllimMA2R5IIDFullK20bis
# 

#### step1 : compute posteriors parameters for y 
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

# REM: y is a DR-dim vector turned into a D x R matrix
#yDR<-matrix(y, ncol=R)
# yDR = obs
## matz is a DRxM matrix turned into an array (D,R,M)
##matzDR<-array(matz,c(D,R,M))

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
gramy <- matrix(NA,nrow=dimay,ncol=dimay) # symmetic
for (ii in 1:dimay) {
  for (jj in ii:dimay) {
    gramy[ii,jj] <- L2scal2normalSMC(postmeany[,ii],
                                  postmeany[,jj],
                                  postcovy[,,ii],
                                  postcovy[,,jj]);
    gramy[jj,ii] <- gramy[ii,jj]
  }
}
qgramy<-quad.form(gramy,Piy)

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
dgllimL2<-function(z){
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
  
  
  #L2dist <- L2normalPreComp(mixy,mixz,qgramy)
  return(L2normalPreCompSMC(mixy,mixz,qgramy))
} # end function distance


####
# distance between summary
summary_obs <- rowMeans(obs)
dsummary <- function(z){
  summary_z <- rowMeans(z)
  return(mean(abs(summary_z - summary_obs)))
}


# common algorithmic parameters
param_algo <- list(nthetas = nparticles, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)


#### For Table S2
filename <- paste0(prefix, "ma2TS2wsmc.d", d, ".n", nobservations, ".wasserstein.RData")
results <- wsmc(wdistance, target, param_algo, maxsimulation = 10^5,savefile = filename)
filename <- paste0(prefix, "ma2TS2wsmc.d", d, ".n", nobservations, ".gllimMW2.RData")
results <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^5,savefile = filename)
filename <- paste0(prefix, "ma2TS2wsmc.d", d, ".n", nobservations, ".gllimL2.RData")
results <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^5,savefile = filename)
#resultsTma2.L2

#### Estimations
## and
#### Plots for FS5
load("xxxx") #put right filename
step<-length(results$thetas_history)
# 2048 distances
dist<-results$distances_history[[step]]
ordist<-order(dist)
# 2048 x 2
thetas<-results$thetas_history[[step]]
thetas100<-thetas[ordist[1:100],]

# Estim
# parameters empirical means from samples, all equal weights...
rowMeans(t(thetas100))
#MW2 0.6182909 0.1764599
# to compute the std on the diagonal of
diag(sqrt(cov(thetas100)) )
# MW2 0.08815274 0.10247397
cor(thetas100)
#[1,] 1.0000000 0.6045888
#[2,] 0.6045888 1.0000000

### L2 
rowMeans(t(thetas100))
#[1] 0.6166277 0.1935900
diag(sqrt(cov(thetas100)) )
#[1] 0.1128262 0.1195882
 cor(thetas100)
#[1,] 1.0000000 0.7270934
#[2,] 0.7270934 1.0000000

### WABC
rowMeans(t(thetas100))
#[1] 0.2423829 0.1171097
 diag(sqrt(cov(thetas100)) )
#[1] 0.2044634 0.2183769
 cor(thetas100)
#[1,] 1.0000000 0.1089243
#[2,] 0.1089243 1.0000000
 

# plot
  df100<-data.frame(thetas100)

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-1,2) + ylim(-0.5,0.75) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=df100, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()





####
filenameE <- paste0(prefix, "ma2wsmc.d", d, ".n", nobservations, ".euclidean.RData")
results <- wsmc(deuclidean, target, param_algo, maxsimulation = 10^6, savefile = filenameE)
load(filename)

filenameWABC <- paste0(prefix, "ma2wsmc.d", d, ".n", nobservations, ".wasserstein.RData")
results <- wsmc(wdistance, target, param_algo, maxsimulation = 10^6, savefile = filenameWABC)
load(filenameWABC)

filenameGMW2 <- paste0(prefix, "ma2wsmc.d", d, ".n", nobservations, ".gllimMW2.RData")
results <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^6, savefile = filenameGMW2)
#load(filenameGMW2)
# 10^5: ok
filenameGMW2105 <- paste0(prefix, "ma2wsmc.d", d, ".n", nobservations, ".gllimMW2105.RData")
results <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^5, savefile = filenameGMW2105)
#

# L2
filenameGL2 <- paste0(prefix, "ma2wsmc.d", d, ".n", nobservations, ".gllimL2.RData")
results <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^6, savefile = filenameGL2)
load(filenameGL2)
# 10^5:ok
filenameGL2 <- paste0(prefix, "ma2wsmc.d", d, ".n", nobservations, ".gllimL2105.RData")
results <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^5, savefile = filenameGL2)


filenameS <- paste0(prefix, "ma2wsmc.d", d, ".n", nobservations, ".summary.RData")
results <- wsmc(dsummary, target, param_algo, maxsimulation = 10^6, savefile = filenameS)
load(filenameS)

# plot check
load(filenameGMW2)
thMW2<- results$thetas_history
plot(density(thMW2[[12]][,1]))
plot(density(thMW2[[12]][,2]))
#
df2048<-data.frame(thMW2[[12]])

v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=0.6, linetype="dashed", size=.5) + geom_hline(yintercept=0.2, linetype="dashed", size=.5)
#v<- v +  geom_vline(xintercept=-0.7, linetype="dashed", size=.5)
v<-v + geom_contour()
v<- v + geom_point(data=df2048, aes(x = X1, y = X2, z=0), size=.3, color="red", alpha=1) 
#v<-v + geom_contour(breaks=c(min(SMA2_df$z), seq(1.24,  max(SMA2_df$z), length.out = 15)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()


##################
####Pas fait
##################

# results <- wsmc_continue(results, savefile = filename, maxtime = 10*60)

# WABC
load(filenameWABC)
thWABC <- results$thetas_history
length(thW2) #20
#dim(thWABC[[17]])

# eucli

load(filename)
theucli <- results$thetas_history
dim(theucli[[17]])
#[1] 2048    2
plot(density(theucli[[17]][,1], from=-1.5,to=0), xlim=c(-1.5,0), ylim=c(0,4), col="black", main="Theta1")

# gllim W2
load(filenameG)
thW2<- results$thetas_history
length(thW2)
#[1] 20
plot(density(thW2[[20]][,1], from=-1.5,to=0), xlim=c(-1.5,0), ylim=c(0,4), col="black", main="Theta1")

# Gllim L2
load(filenameGL2)
thL2<- results$thetas_history
length(thL2) #15

load(filenameS)
thS<- results$thetas_history
length(thS)
#[1] 20
plot(density(thS[[20]][,1], from=-1.5,to=0), xlim=c(-1.5,0), ylim=c(0,4.2), col="black", main="Theta1")


plot(density(results$thetas_history[[3]][,1], from=-1.5,to=0), xlim=c(-1.5,0), ylim=c(0,4), col="black", main="Theta1")


# Gllim mixture margins
gm1bis<-dnorm(thetaseq[1,],postmeany[1,], postcovy[1,1])
gm2bis<-dnorm(thetaseq[2,],postmeany[2,], postcovy[2,2])

######################################
#True Posterior
# yobs is 2iR x 1 needs to be turned into 2x iR
# thetaseq is 2 x nbsample
#PostPdfNLM<-function(thetaseq,yobs,sigmalik,iR){
thetaseq<-matrix(0,2,1000)
thetaseq[1,]<-seq(-2.5,1,length.out=1000)
thetaseq[2,]<-seq(-1.5,1.5,length.out=1000)
trupbis<-PostPdfNLM(thetaseq,obs, sigmalik,100)

#############################################
## PLOTS posterior margins
#############################################

# margin 1 theta1
plot(density(thW2[[20]][,1], from=-1.5,to=0), xlim=c(-1.5,0), ylim=c(0,4.2), col="black", main="Theta1")
lines(density(thS[[20]][,1], from=-1.5,to=0), xlim=c(-1.5,0), col="red")
lines(density(theucli[[17]][,1], from=-1.5,to=0), xlim=c(-1.5,0), col="red", lty=2)
lines(density(thL2[[15]][,1] , from=-1.5,to=0), xlim=c(-1.5,0), col="green")
lines(density(thWABC[[20]][,1], from=-1.5,to=0), xlim=c(-1.5,0), col="blue")
lines(thetaseq[1,],trupbis$marg1, col="purple")
lines(thetaseq[1,],gm1bis, col="green", lty=2)
#lines(density(samplemixgllimBBM[,1], from=0,to=5), xlim=c(0,5), col="green", lty=2)
#
abline(v=-0.71, lty=2)

# margin 2 theta2
plot(density(thW2[[20]][,2], from=-0.5,to=0.75), xlim=c(-0.5,0.75), ylim=c(0,4), col="black", main="Theta2")
lines(density(thS[[20]][,2], from=-0.5,to=0.75), xlim=c(-0.5,0.75), col="red")
lines(density(theucli[[17]][,2], from=-0.5,to=0.75), xlim=c(-0.5,0.75), col="red", lty=2)
lines(density(thL2[[15]][,2] , from=-0.5,to=0.75), xlim=c(-0.5,0.75), col="green")
lines(density(thWABC[[20]][,2], from=-0.5,to=0.75), xlim=c(-0.5,0.75), col="blue")
lines(thetaseq[2,],trupbis$marg2, col="purple")
lines(thetaseq[2,],gm2bis, col="green", lty=2)
#lines(density(samplemixgllimBBM[,2], from=0,to=5), xlim=c(0,5), col="green", lty=2)
#
abline(v=0.09, lty=2)


########TIMES
### eg w_to_post_summary
###############


#
#
# load(filename)
# plot_threshold_time(results) + geom_point()
# mle <- rowMeans(obs)
# plot_bivariate(results, 1, 2, from = 10) + geom_vline(xintercept = mle[1]) + geom_hline(yintercept = mle[2])
# plot_marginal(results, 1, from = 10)

# library(microbenchmark)
# microbenchmark(deuclidean(target$simulate(true_theta)), times = 1000)
