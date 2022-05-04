library(winference)
registerDoParallel(cores = detectCores())
#rm(list = ls())
setmytheme()

#my_colors <- get_my_colors()
# or
gllim_colors <- get_gllim_colors()

set.seed(11)

doRun <- FALSE
max_time <- 30*60
d <- 2  # D 


target <- get_bbm(d)
#target$parameters$tau <- 5
nobservations <- 100  # R or 
nparticles <- 2048
p <- 1
prefix <- ""

# contain obs the observations D x R
#obsfile <- paste0(prefix, "ma2data.d", d, ".n", nobservations, ".RData")
#load(obsfile)
# here:
# to get the same observations as in Fig S1 and S2 reduced to 100 iid (500 initially)
load("~/florence/local/THESES/FabienBOUX/DATA/GITaout2021/GLLiM-ABC/data-et-al4examples/dataBBM.RData")
#rm(bbmparams) # not needed for smc
#rm(bbmdata)
obs<-matrix(bbmobs[1:200,], ncol=nobservations)


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
load("~/florence/local/THESES/FabienBOUX/DATA/ABC-GLLiM-Package-April21-MyTest/TEST-MA2-mai2021/modgllimBBMIIDFullK100n100.RData")
modg<-modgllimBBMIIDFullK100n100 
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

wthr<-0.001

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
gramy <- matrix(NA,nrow=dimay,ncol=dimay) # symmetric
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
# mixture_rmixmod(nclust=5) # default
# here nclust set to 2 otherwise induces warnings and errors
param_algo <- list(nthetas = nparticles, nmoves = 1, proposal = mixture_rmixmod(nclust=2),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)

####
filenameE <- paste0(prefix, "bbmwsmc.d", d, ".n", nobservations, ".euclidean.RData")
results <- wsmc(deuclidean, target, param_algo, maxsimulation = 10^5, savefile = filenameE)
results_euclidean <- results
load(filename)

filenameWABC <- paste0(prefix, "bbmwsmc.d", d, ".n", nobservations, ".wasserstein.RData")
results <- wsmc(wdistance, target, param_algo, maxsimulation = 10^5, savefile = filenameWABC)
results_wasserstein <- results
load(filenameWABC)
# 10^6
filenameWABC <- paste0(prefix, "bbmwsmc.d", d, ".n", nobservations, ".wasserstein106.RData")
results <- wsmc(wdistance, target, param_algo, maxsimulation = 10^6, savefile = filenameWABC)
results_wasserstein106 <- results

# BBM MW2 10^5
filenameGMW2 <- paste0(prefix, "bbmwsmc.d", d, ".n", nobservations, ".gllimMW2.RData")
results <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^5, savefile = filenameGMW2)
results_gllimMW2 <- results
#load(filenameGMW2)

# BBM MW2 10^6
filenameGMW2 <- paste0(prefix, "bbmwsmc.d", d, ".n", nobservations, ".gllimMW2106.RData")
results <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^6, savefile = filenameGMW2)
results_gllimMW2106 <- results

# BBM L2 10^5 :ok  with nclust=2
filenameGL2 <- paste0(prefix, "bbmwsmc.d", d, ".n", nobservations, ".gllimL2.RData")
results <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^5, savefile = filenameGL2)
results_gllimL2 <- results
load(filenameGL2)
# L2 10^6 : des erreurs et warning??
# re done : pb a partir du step 10! (143511 dist computed): stop at 11
# try with set.seed(12): idem a step 4
# ----> put nclust to 2: warnings and errors at step 10 but goes on until 17
filenameGL2 <- paste0(prefix, "bbmwsmc.d", d, ".n", nobservations, ".gllimL2106.RData")
results <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^6, savefile = filenameGL2)
results_gllimL2106 <- results

filenameS <- paste0(prefix, "bbmwsmc.d", d, ".n", nobservations, ".summary.RData")
results <- wsmc(dsummary, target, param_algo, maxsimulation = 10^5, savefile = filenameS)
results_summary <- results
load(filenameS)

# plot check
load(filenameGMW2)
thMW2<- results$thetas_history
plot(density(thMW2[[12]][,1]))
plot(density(thMW2[[12]][,2]))
#
df2048<-data.frame(thMW2[[12]])

# 50 best
titi<-results$distances_history[[12]]
toto<-order(titi)
thMW250<-thMW2[[12]][toto[1:50],]
plot(density(thMW250[,1]))

### To compute MSE et al for Table Sup S1: 
### reduce samples to 50 best for WABC, SMCMW2 , SMC-L2
### x2 10^5 and 10^6

th2048<- results_wasserstein$thetas_history
step<-length(th2048)
titi<-results_wasserstein$distances_history[[step]]
toto<-order(titi)
thWABC50M105<-th2048[[step]][toto[1:50],]
#
th2048<- results_wasserstein106$thetas_history
step<-length(th2048)
titi<-results_wasserstein106$distances_history[[step]]
toto<-order(titi)
thWABC50M106<-th2048[[step]][toto[1:50],]

th2048<- results_gllimMW2$thetas_history
step<-length(th2048)
titi<-results_gllimMW2$distances_history[[step]]
toto<-order(titi)
thMW250M105<-th2048[[step]][toto[1:50],]
#
th2048<- results_gllimMW2106$thetas_history
step<-length(th2048)
titi<-results_gllimMW2106$distances_history[[step]]
toto<-order(titi)
thMW250M106<-th2048[[step]][toto[1:50],]

th2048<- results_gllimL2$thetas_history
step<-length(th2048)
titi<-results_gllimL2$distances_history[[step]]
toto<-order(titi)
thL250M105<-th2048[[step]][toto[1:50],]
#
th2048<- results_gllimL2106$thetas_history
step<-length(th2048)
titi<-results_gllimL2106$distances_history[[step]]
toto<-order(titi)
thL250M106<-th2048[[step]][toto[1:50],]

#### compute means and MSE 
#######################################################
# Empirical means and rmse for each ABC samples etc.
#######################################################

## Empirical means for each ABC samples etc....
colMeans(thWABC50M105)
#[1] 1.8710045 1.7899847 0.7550683 0.9391245 1.7668801
colMeans(thWABC50M106)
#[1] 1.5458223 1.7316641 0.7687267 0.8211230 1.5656387
 colMeans(thMW250M105)
#[1] 1.1485440 1.8489203 0.6409566 0.5241444 1.4166453
 colMeans(thMW250M106)
#[1] 1.3164512 1.6555310 0.6842534 0.7076463 1.4949772
 colMeans(thL250M105)
#[1] 2.133420 2.581240 2.443232 2.344856 3.074158
 colMeans(thL250M106)
#[1] 1.874628 2.947641 3.715076 3.062133 3.429933


## # compute RMSE
Mm<-50 #pts
truemat<-matrix(1,5,Mm)
msemat<-(t(thMW250M105) -truemat)^2
rmsem<-sqrt(rowMeans(msemat))

# W2 10^5 : 0.5749088 1.0718358 0.5625080 0.5584576 0.5021307
# W2 10^6: 0.5080650 0.8516947 0.3991505 0.4208917 0.5473526
# L2 10^5  1.393789 1.985199 1.996444 1.735478 2.216891
# L2 10^6: 1.437933 2.279484 2.908325 2.147152 2.512803
# WABC 10^5: 1.1579165 1.2008127 0.5022536 0.5082945 0.8180683
# WABC 10^6: 0.7726928 0.8793367 0.3725209 0.3650109 0.6210949


save(list = ls(), file = "env_entier15mars2022.RData")


#############################################################################
# PLOTS for paper 
#############################################################################
wsmc.euclidean.df <- wsmc_to_dataframe(results_euclidean) %>% filter(step == length(results_euclidean$thetas_history))
wsmc.summary.df <- wsmc_to_dataframe(results_summary) %>% filter(step == length(results_summary$thetas_history))
wsmc.wasserstein.df <- wsmc_to_dataframe(results_wasserstein) %>% filter(step == length(results_wasserstein$thetas_history))
wsmc.gllimMW2.df <- wsmc_to_dataframe(results_gllimMW2) %>% filter(step == length(results_gllimMW2$thetas_history))
wsmc.gllimL2.df <- wsmc_to_dataframe(results_gllimL2) %>% filter(step == length(results_gllimL2$thetas_history))
#
wsmc.gllimMW2106.df <- wsmc_to_dataframe(results_gllimMW2106) %>% filter(step == length(results_gllimMW2106$thetas_history))
wsmc.wasserstein106.df <- wsmc_to_dataframe(results_wasserstein106) %>% filter(step == length(results_wasserstein106$thetas_history))
wsmc.gllimL2106.df <- wsmc_to_dataframe(results_gllimL2106) %>% filter(step == length(results_gllimL2106$thetas_history))

g1 <- ggplot(wsmc.wasserstein.df, aes(x = Theta1)) + geom_density(aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5, linetype = "dashed") + theme(aspect.ratio=.5)
g1<- g1 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=.5)
#g1 <- ggplot(wsmc.summary.df, aes(x = Theta1)) + geom_density(aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
##g1 <- ggplot(data.frame(posterior_sample), aes(x = X1)) + geom_density(aes(y = ..density.., fill = " ", colour = "Posterior"), alpha = 0.5)
###g1<-ggplot(data.frame(thetaseq[1,], trupbis$marg1), aes(x=thetaseq.1..., y=trupbis.marg1)) + geom_point()
##g1 <- ggplot(data.frame(pos), aes(x = X1)) + geom_density(aes(y = ..density.., fill = " ", colour = "Posterior"), alpha = 0.5)
##g1 <- g1 + geom_density(data=wsmc.summary.df, aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
#g1 <- g1 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g1 <- g1 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5, linetype = "dashed")
g1 <- g1 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "L2"), alpha = 0.5, linetype = "dashed")
g1 <- g1 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5)
g1 <- g1 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "L2"), alpha = 0.5)
g1 <- g1 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
## N=100 good
## g1 <- g1 + geom_density(data=wsmc.gllimW2100.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
## g1 <- g1 + geom_density(data=wsmc.gllimL2100.df, aes(y = ..density.., fill = "", colour = "Rej. Summary"), alpha = 0.5)
###g1 <- g1 + geom_density(data=rej.wasserstein.df, aes(y = ..density.., fill = "Rej. Summary", colour = "Rej. Summary"), alpha = 0.5)
#g1 <- g1 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g1<- g1 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
#g1 <- g1 + scale_color_manual(name = "", values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors)
g1 <- g1 + scale_color_manual(values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors) + theme(legend.position = "none",title =element_text(size=18, face='bold'))
g1 <- g1 + ggtitle("Theta1") + xlab("") + ylab("Density") + theme(plot.title = element_text(hjust = 0.5))
# for old style comment:
#g1 <- g1 + xlab(expression(theta[1]))
#g1 <- g1 + geom_label(data = data.frame(x = c(4,4,4), y = c( 0.6, 0.9, 1.2),
#                                  method = c( "WABC", "Gllim-MW2", "Gllim-L2")),
#                                  aes(x = x, y = y, colour = c( "marginal","MW2","L2"), label = method), size = 8) + theme(legend.position = "none")
##g1 <- g1 + geom_label(data = data.frame(x = c(3.3, 2.9, 3, 1.6, 3,4,4.5,4.5), y = c(0.4, 0.55, 0.7, 0.85, 1, 0.85,0.55,1),
##                                        method = c( "Summary", "WABC", "Euclidean", "Gllim-MW2", "Gllim-L2","Gllim-MW2b", "WABCb", "Gllim-L2b")),
##                      aes(x = x, y = y, colour = c( "Summary", "Wasserstein", "Euclidean","marginal","Rej. Summary","Posterior","W + constraint","Swap"), label = method), size = 8) + theme(legend.position = "none")
g1
ggsave(filename = paste0(prefix, "bbm.posterior1.pdf"), plot = g1, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "bbm.posterior1.png"), plot = g1, width = 7.5, height = 7, dpi = 150)
# old style
ggsave(filename = paste0(prefix, "bbm.posterior1v1.pdf"), plot = g1, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "bbm.posterior1v1.png"), plot = g1, width = 7.5, height = 7, dpi = 150)
# for paper:
ggsave(filename = paste0(prefix, "bbm.smc.posterior1.pdf"), plot = g1, width = 7.5, height = 7)

 
# Theta 2 style for rev1 
g2 <- ggplot(wsmc.wasserstein.df, aes(x = Theta2)) + geom_density(aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5, linetype = "dashed")+ theme(aspect.ratio=.5)
g2<- g2 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=.5)
#g2 <- ggplot(wsmc.summary.df, aes(x = Theta2)) + geom_density(aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
#g2 <- g2 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5, linetype = "dashed")
g2 <- g2 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "L2"), alpha = 0.5, linetype = "dashed")
g2 <- g2 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "L2"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
#g2 <- g2 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g2<- g2 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
g2 <- g2 + scale_color_manual(values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors) + theme(legend.position = "none",title =element_text(size=18, face='bold'))
g2 <- g2 + ggtitle("Theta2") + xlab("") + ylab("Density") + theme(plot.title = element_text(hjust = 0.5))
# for old style see next plot
g2
ggsave(filename = paste0(prefix, "bbm.smc.posterior2.pdf"), plot = g2, width = 7.5, height = 7)


# Theta 2 New Style (winference style)
g2 <- ggplot(wsmc.wasserstein.df, aes(x = Theta2)) + geom_density(aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
#g2 <- ggplot(wsmc.summary.df, aes(x = Theta2)) + geom_density(aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
#g2 <- g2 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "Rej. Summary"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "Posterior"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "Swap"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "W + constraint"), alpha = 0.5)
#g2 <- g2 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g2<- g2 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
g2 <- g2 + scale_color_manual(name = "", values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors)
g2 <- g2 + xlab(expression(theta[2]))
g2 <- g2 + geom_label(data = data.frame(x = c( 2.5, 1.5, 3,4, 4.1,4.1), y = c( 0.65, 0.95, 1.1, 0.95, 0.65,1.1),
method = c( "WABC", "Gllim-MW2", "Gllim-L2","Gllim-MW2b", "WABCb", "Gllim-L2b")),
aes(x = x, y = y, colour = c( "Wasserstein","marginal","Rej. Summary","Posterior","W + constraint","Swap"), label = method), size = 8) + theme(legend.position = "none")
#g2 <- g2 + geom_label(data = data.frame(x = c(4, 2.5, 3, 1.5, 3,4, 4.1,4.1), y = c(0.5, 0.65, 0.8, 0.95, 1.1, 0.95, 0.65,1.1),
 #                                       method = c( "Summary", "WABC", "Euclidean", "Gllim-MW2", "Gllim-L2", "Gllim-MW2b", "WABCb", "Gllim-L2b")),
 #                     aes(x = x, y = y, colour = c( "Summary", "Wasserstein", "Euclidean","marginal","Rej. Summary","Posterior", "W + constraint","Swap"), label = method), size = 8) + theme(legend.position = "none")
g2
ggsave(filename = paste0(prefix, "bbm.posterior2.pdf"), plot = g2, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "bbm.posterior2.png"), plot = g2, width = 7.5, height = 7, dpi = 150)

# Theta 3 for Rev1
g3 <- ggplot(wsmc.wasserstein.df, aes(x = Theta3)) + geom_density(aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5, linetype = "dashed")+ theme(aspect.ratio=.5)
g3<- g3 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=.5)
#g3 <- g3 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g3 <- g3 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5, linetype = "dashed")
g3 <- g3 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "L2"), alpha = 0.5, linetype = "dashed")
g3 <- g3 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5)
g3 <- g3 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "L2"), alpha = 0.5)
g3 <- g3 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
#g3 <- g3 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g3<- g3 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
g3 <- g3 + scale_color_manual(values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors) + theme(legend.position = "none",title =element_text(size=18, face='bold'))
g3 <- g3 + ggtitle("Theta3") + xlab("") + ylab("Density") + theme(plot.title = element_text(hjust = 0.5))
#
g3
ggsave(filename = paste0(prefix, "bbm.smc.posterior3.pdf"), plot = g3, width = 7.5, height = 7)

## other plot

g3 <- ggplot(wsmc.wasserstein.df, aes(x = Theta3)) + geom_density(aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
#g3 <- ggplot(wsmc.summary.df, aes(x = Theta3)) + geom_density(aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
#g3 <- g3 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g3 <- g3 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
g3 <- g3 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "Rej. Summary"), alpha = 0.5)
g3 <- g3 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "Posterior"), alpha = 0.5)
g3 <- g3 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
g3 <- g3 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "W + constraint"), alpha = 0.5)
#g3 <- g3 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g3<- g3 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
g3 <- g3 + scale_color_manual(name = "", values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors)
g3 <- g3 + xlab(expression(theta[3]))
g3 <- g3 + geom_label(data = data.frame(x = c(2.5, 2.4, 1.8, 4,4.2,4), y = c(0.95, 0.65,  1.1, 0.45,0.95, 1.1),
                                        method = c( "WABC",  "Gllim-MW2", "Gllim-L2", "Gllim-MW2b", "WABCb", "Gllim-L2b")),
                      aes(x = x, y = y, colour = c( "Wasserstein","marginal","Rej. Summary","Posterior", "W + constraint","Summary"), label = method), size = 8) + theme(legend.position = "none")
#g3 <- g3 + geom_label(data = data.frame(x = c(3, 2.5, 3, 1.8, 3,4,4), y = c(0.5, 0.65, 0.8, 0.95, 1.1,0.95, 0.65),
#                                        method = c( "Summary", "WABC", "Euclidean", "Gllim-MW2", "Gllim-L2", "Gllim-MW2b", "WABCb", "Gllim-L2b")),
#                      aes(x = x, y = y, colour = c( "Summary", "Wasserstein", "Euclidean","marginal","Rej. Summary","Posterior", "W + constraint","Summary"), label = method), size = 8) + theme(legend.position = "none")
g3
ggsave(filename = paste0(prefix, "bbm.posterior3.pdf"), plot = g3, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "bbm.posterior3.png"), plot = g3, width = 7.5, height = 7, dpi = 150)


# Theta 4 for rev 1 
g4 <- ggplot(wsmc.wasserstein.df, aes(x = Theta4)) + geom_density(aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5, linetype = "dashed")+ theme(aspect.ratio=.5)
g4<- g4 + theme_bw() + theme( panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=.5)
#g4 <- ggplot(wsmc.summary.df, aes(x = Theta4)) + geom_density(aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
#g4 <- g4 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g4 <- g4 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5, linetype = "dashed")
g4 <- g4 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "L2"), alpha = 0.5, linetype = "dashed")
g4 <- g4 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5)
g4 <- g4 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "L2"), alpha = 0.5)
g4 <- g4 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
#g4 <- g4 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g4<- g4 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
g4 <- g4 + scale_color_manual(values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors) + theme(legend.position = "none",title =element_text(size=18, face='bold'))
g4 <- g4 + ggtitle("Theta4") + xlab("") + ylab("Density") + theme(plot.title = element_text(hjust = 0.5))
#
g4
ggsave(filename = paste0(prefix, "bbm.smc.posterior4.pdf"), plot = g4, width = 7.5, height = 7)


# other style

g4 <- ggplot(wsmc.wasserstein.df, aes(x = Theta4)) + geom_density(aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
#g4 <- ggplot(wsmc.summary.df, aes(x = Theta4)) + geom_density(aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
#g4 <- g4 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g4 <- g4 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5)
g4 <- g4 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "L2"), alpha = 0.5)
g4 <- g4 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5)
g4 <- g4 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "L2"), alpha = 0.5)
g4 <- g4 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
#g4 <- g4 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g4<- g4 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
g4 <- g4 + scale_color_manual(name = "", values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors)
g4 <- g4 + xlab(expression(theta[4]))
g4 <- g4 + geom_label(data = data.frame(x = c( 2.5,  1.7, 3,4,4,4), y = c( 0.85,  1.1, 0.65,1.1,0.85, 0.45),
                                        method = c( "WABC", "Gllim-MW2", "Gllim-L2", "Gllim-MW2b", "WABCb", "Gllim-L2b")),
                      aes(x = x, y = y, colour = c(  "Wasserstein", "marginal","Rej. Summary","Posterior","W + constraint","Summary"), label = method), size = 8) + theme(legend.position = "none")
#g4 <- g4 + geom_label(data = data.frame(x = c(2.9, 2.5, 3, 1.7, 3,4,4), y = c(1.3, 0.65, 0.85, 1.1, 0.45,1.1,0.65),
#                                       method = c( "Summary", "WABC", "Euclidean", "Gllim-MW2", "Gllim-L2", "Gllim-MW2b", "WABCb", "Gllim-L2b")),
#                      aes(x = x, y = y, colour = c( "Summary", "Wasserstein", "Euclidean","marginal","Rej. Summary","Posterior","W + constraint","Swap"), label = method), size = 8) + theme(legend.position = "none")
g4
ggsave(filename = paste0(prefix, "bbm.posterior4.pdf"), plot = g4, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "bbm.posterior4.png"), plot = g4, width = 7.5, height = 7, dpi = 150)


# Theta 5 for Rev 1
g5 <- ggplot(wsmc.wasserstein.df, aes(x = Theta5)) + geom_density(aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5, linetype="dashed") + theme(aspect.ratio=.5)
g5<- g5 + theme_bw() + theme( panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=.5)
#g5 <- ggplot(wsmc.summary.df, aes(x = Theta5)) + geom_density(aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
#g5 <- g5 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g5 <- g5 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5, linetype="dashed")
g5 <- g5 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "L2"), alpha = 0.5, linetype="dashed")
g5 <- g5 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "MW2"), alpha = 0.5)
g5 <- g5 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "L2"), alpha = 0.5)
g5 <- g5 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
#g5 <- g5 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g5<- g5 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
g5 <- g5 + scale_color_manual(values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors) + theme(legend.position = "none",title =element_text(size=18, face='bold'))
g5 <- g5 + ggtitle("Theta5") + xlab("") + ylab("Density") + theme(plot.title = element_text(hjust = 0.5))
#
g5
ggsave(filename = paste0(prefix, "bbm.smc.posterior5.pdf"), plot = g5, width = 7.5, height = 7)

#Other plot
g5 <- ggplot(wsmc.wasserstein.df, aes(x = Theta5)) + geom_density(aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
#g5 <- ggplot(wsmc.summary.df, aes(x = Theta5)) + geom_density(aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
#g5 <- g5 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g5 <- g5 + geom_density(data=wsmc.gllimMW2.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
g5 <- g5 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "Rej. Summary"), alpha = 0.5)
g5 <- g5 + geom_density(data=wsmc.gllimMW2106.df, aes(y = ..density.., fill = " ", colour = "Posterior"), alpha = 0.5)
g5 <- g5 + geom_density(data=wsmc.gllimL2106.df, aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
g5 <- g5 + geom_density(data=wsmc.wasserstein106.df, aes(y = ..density.., fill = " ", colour = "W + constraint"), alpha = 0.5)
#g5 <- g5 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g5<- g5 +  geom_vline(xintercept=1, linetype="dashed", size=.5)
g5 <- g5 + scale_color_manual(name = "", values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors)
g5 <- g5 + xlab(expression(theta[5]))
g5 <- g5 + geom_label(data = data.frame(x = c(2.8, 2,  2.8, 4.1,4.2,4.4), y = c( 1.5,  2, 0.78, 2,1.5,0.78),
                                        method = c(  "WABC",  "Gllim-MW2", "Gllim-L2", "Gllim-MW2b", "WABCb", "Gllim-L2b")),
                      aes(x = x, y = y, colour = c(  "Wasserstein", "marginal","Rej. Summary","Posterior","W + constraint","Summary"), label = method), size = 8) + theme(legend.position = "none")
#g5 <- g5 + geom_label(data = data.frame(x = c(4, 3, 3, 1.9, 3,4.1,4.5), y = c(0.5, 1.25, 1, 1.5, 0.78,1.5, 1.25),
#                                        method = c( "Summary", "WABC", "Euclidean", "Gllim-MW2", "Gllim-L2", "Gllim-MW2b", "WABCb", "Gllim-L2b")),
#                      aes(x = x, y = y, colour = c( "Summary", "Wasserstein", "Euclidean","marginal","Rej. Summary","Posterior","W + constraint","Swap"), label = method), size = 8) + theme(legend.position = "none")
g5
ggsave(filename = paste0(prefix, "bbm.posterior5.pdf"), plot = g5, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "bbm.posterior5.png"), plot = g5, width = 7.5, height = 7, dpi = 150)




#OLD 
g2 <- ggplot(data.frame(posterior_sample), aes(x = X2)) + geom_density(aes(y = ..density.., fill = " ", colour = "Posterior"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.summary.df, aes(y = ..density.., fill = " ", colour = "Summary"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = " ", colour = "Euclidean"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.gllimW2.df, aes(y = ..density.., fill = " ", colour = "marginal"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.gllimL2.df, aes(y = ..density.., fill = "", colour = "Rej. Summary"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = " ", colour = "Wasserstein"), alpha = 0.5)
g2 <- g2 + scale_color_manual(name = "", values = gllim_colors) + scale_fill_manual(name = "", values = gllim_colors)
g2 <- g2 + geom_label(data = data.frame(x = c(0.43,-0.25, -0.23, 0.55, 0.5,-0.25), y = c(3.8, 3, 1.8, 0.7, 2.6, 3.6),
                                        method = c("Summary", "Posterior", "WABC","Euclidean","Gllim-MW2","Gllim-L2")),
                      aes(x = x, y = y, colour = c("Summary", "Posterior","Euclidean","Wasserstein","marginal","Rej. Summary"), label = method), size = 8) + theme(legend.position = "none")
g2 <- g2 + xlab(expression(theta[2]))
g2
ggsave(filename = paste0(prefix, "bbm.posterior2.pdf"), plot = g2, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "bbm.posterior2.png"), plot = g2, width = 7.5, height = 7, dpi = 150)



###### OLD

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
