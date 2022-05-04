#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)

ny<-10
target <- get_itdmix(ny)
# number of  iid observations (R=1)
##### see itdmix_generate_data.R to generate observations...
nobservations <- 1
d <- 10  # D 
p<-1

nparticles <- 2048

# obs is the observations D x 1
#obsfile <- paste0(prefix, "ma2data.d", d, ".n", nobservations, ".RData")
#load(obsfile)
# here:
# load file dataITDmixdf3.RData in ABC-GLLIMxxxxMyTest...
# use observation: ytargetITD is 1 x 10
##### use an already Learned GLLiM model: modgllimITDgmix1K38

obs<-t(ytargetITD)

# function to simulate data
target$simulate <- function(theta){
  return(target$robservation(theta,target$parameters,ny))
}

# distances :

# wasserstein distance: NA here ???
wdistance <- get_transport_to_y(obs, p = p)
#

########## prep for SMC with  distance W2 and L2
# obs = observations is 10 x 1

## Gllim model used
modg<-modgllimITDgmix1K38

# for SMC
param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)
##t = proc.time()
#filename <- paste0(prefix, "gandk3wsmc.n", nobservations, ".RData")
#results <- wsmc(compute_d, target, param_algo, savefile = filename, maxsim = 2.4e6)
##t = proc.time() - t
#load(filename)
##results <- wsmc_continue(results, savefile = filename, maxtime = 14*60*60)


#### SMC-MW2 : Attention specific to D=1 or not ??
# obs = observations is a vector not a matri x 1 x 250

#obsmat<-matrix(obs, ncol=nobservations)


#### step1 : compute posteriors parameters for y 
# notational shortcut: gllim direct parameters
modeleg<-modg$mod
Aa<-modeleg$A
Sigmaa<-modeleg$Sigma
ca<-modeleg$c
ba<-modeleg$b
#pia<-modeleg$pi
covstara<-modg$covstar
#logdetVR<-modg$logdetVR
invGammaa<-modg$invGamma
#invSigmaa<-modg$invSigma

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
# xllim inverse model
postweighty<-gllim_inverse_map(obs,modeleg)$alpha
Piy=t(postweighty)
##Kdim<-length(Piy)
# in case of very low weights
seuil<-min(wthr,sort(Piy, decreasing=TRUE)[3])
leftky<-seq(1,K)[Piy>seuil]
#leftky<-seq(1,Kdim)[Piy>wthr]
# dimay: number of Gaussians whose weight is above the threshold
dimay<-sum(Piy>seuil)
Piy<-Piy[Piy>seuil]

postcovy<-array(0,c(L,L,dimay))
postmeany<-matrix(0,L,dimay)

postcovy[,,1:dimay]=covstara[,,leftky]

Asybs<-NULL
for (k in leftky){
  # ISOTROPIC case
  Asybs<-cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(obs-ba[,k]) +invGammaa[,,k]%*%ca[,k]) )
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
  leftky<-c(leftky,leftky)
} 


tracecost<-matrix(0,dimay,K)
for (ii in 1:dimay) {
  sigma1<-postcovy[,,ii]
  for (jj in 1:K){
    sigma2<-covstara[,,jj]
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
  ##ICI
  z<-t(z)
  postweightz<-gllim_inverse_map(z,modeleg)$alpha
  Piz=t(postweightz)
  # in case of very low weights
  seuil<-min(wthr,sort(Piz, decreasing=TRUE)[3])
  leftkz<-seq(1,K)[Piz>seuil]
  #leftkz<-seq(1,Kdim)[Piz>wthr]
  dimaz<-sum(Piz>seuil)
  Piz<-Piz[Piz>seuil]
  
  postcovz<-array(0,c(L,L,dimaz))
  postmeanz<-matrix(0,L,dimaz)
  
  postcovz[,,1:dimaz]=covstara[,,leftkz]
  
  Asybs<-NULL
  for (k in leftkz){
    # ISOTROPIC case
    Asybs<-cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(z-ba[,k]) +invGammaa[,,k]%*%ca[,k]) )
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
  
  ##
  mixy<-list("Mu"=postmeany, "S"=postcovy, "Pi"=Piy)
  mixz<-list("Mu"=postmeanz, "S"=postcovz, "Pi"=Piz)
  
  # tracecost has to be reduced to leftky x leftkz before the call
  #tracecostz<-matrix(0,dimay,dimaz)
  # t() if tracecost is 1 x xxxx 
  #if (dim(tracecost)[1]==1){
  #tracecostz<-t(tracecost[,leftkz])}
  # else 
  {tracecostz<-tracecost[,leftkz]}
  
  #MW2dist <- was2mixPreComp(mixy,mixz,tracecostz)
  return(was2mixPreComp(mixy,mixz,tracecostz) )
} # end function distance MW2

# 
# 10^5
filename <- paste0(prefix, "itdmixK38wsmc.d", d, ".n", nobservations, ".gllimMW2.RData")
results <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^5,savefile = filename)
#total # distances calculated: 124336 (for this step: 57398)
#total time spent: 2302.212 seconds (for this step: 833.747 seconds)

# 10^5 nclust=8 dans algo_param for 8 branches
## Pas mieux que nclust=5....
param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(nclust=8),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)
#
filename <- paste0(prefix, "itdmixK38105nc8wsmc.d", d, ".n", nobservations, ".gllimMW2.RData")
results <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^5,savefile = filename)
#step 5... running until # distances calculated >= 1e+05 
#acceptance rates: 47.31445 %, threshold = 0.3926264 , min. dist. = 0.05725834 
#total # distances calculated: 113213 (for this step: 50384)
#total time spent: 1775.671 seconds (for this step: 806.595 seconds)

# 10^6 nclust=5 (default) dans algo_param 
param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)
#
filename <- paste0(prefix, "itdmixK38106wsmc.d", d, ".n", nobservations, ".gllimMW2.RData")
results <- wsmc(dgllimMW2, target, param_algo, maxsimulation = 10^6,savefile = filename)
#
#step 9... running until # distances calculated >= 1e+06 
#acceptance rates: 45.70312 %, threshold = 0.1561745 , min. dist. = 0.04070304 
#total # distances calculated: 1021032 (for this step: 404991)
#total time spent: 4209.796 seconds (for this step: 1652.771 seconds)

## Plot for Revision

step<-length(results$thetas_history)
# 2048 distances
dist<-results$distances_history[[step]]
ordist<-order(dist)
# 2048 x 2
thetas<-results$thetas_history[[step]]
thetas1000<-thetas[ordist[1:1000],]

dfMW2ITDSMC<-data.frame(thetas1000)

dfconfig<-data.frame(matrix(c(-0.5,0,0.5,0,0,-0.5,0,0.5,1.5,1), 5,2, byrow=T))


####### GLLIM-MW2-ABC-SMC 
m <- ggplot(dfMW2ITDSMC, aes(x = X1, y = X2)) + theme_gray() +
  #  geom_point(size=.3) + 
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
#m<- m + geom_density_2d(size = .5, aes(color = ..level..)) # + scale_color_gradient(low = "blue", high = "red")
m<- m +  geom_point(data=dfMW2ITDSMC, size=1, color=rgb(0,0,.1,alpha = 0.5), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m


# To save a  pdf
#pdf(file = "name.pdf", width = 5.5, height = 4)
#[...]
#dev.off()

#### Plots for FS5
#load("xxxx") #put right filename
step<-length(results$thetas_history)
# 2048 distances
#dist<-results$distances_history[[step]]
#ordist<-order(dist)
# 2048 x 2
thetas<-results$thetas_history[[step]]
#thetas100<-thetas[ordist[1:100],]

#ICI

# load(filename)
# wsmc.df <- wsmc_to_dataframe(results)
# nsteps <- max(wsmc.df$step)
#
# # plot_bivariate_polygon(results, 1, 2)
# # plot_bivariate_polygon(results, 3, 4)
#
# library(gridExtra)
# grid.arrange(plot_marginal_time(results, 1),
#   plot_marginal_time(results, 2),
#   plot_marginal_time(results, 3),
#   plot_marginal_time(results, 4), nrow = 2)




step<-length(results$thetas_history)
# 2048 distances
#dist<-results$distances_history[[step]]
#ordist<-order(dist)
# 2048 x 2
thetas<-results$thetas_history[[step]]



