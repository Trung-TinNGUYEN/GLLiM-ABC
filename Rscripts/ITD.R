
#############################################################################################
####  Example ITD + absolute value + Student noise model with dof=1 , D=10, L=2
####  goal: find the location of a source localised at 2D coordinates x12 from a so-called 
####        "ITD" of dimension D (see details in the paper)
##############################################################################################
#### Functions specific to this example, to create data sets, computation of the 
#### true posterior, all saved in a dataITD.Rdata file in the end


#################################################################################################
###
### Step 1: Setting: target y (observations to be inverted, dimension is D) ;
###         Data sets simulation, one for GLLiM, one for ABC if wanted
###
#################################################################################################
# input: true parameters, simulated parameters etc...
# output: ytarget; (thetasimu1, ysimu1) size N; (thetasimu2, ysimu2) size M ;
#  true posterior computed numerically if possible


### To simulate a vector of dim D with |ITD| in the mean and a Student distribution
# ie simulation of D ITD with the same source 2D position given by x12 
# with a Student noise ie isotropic scale matrix with sigmaT=0.01 (little noise)
# and dofT=3 
# micros positions, configuration setting, eg. micr1=c(-0.5,0) and micr2=c(0.5,0)
# OUTPUT: a 1 x D matrix

#############################
# ITD model  with micr1 =c(-0.5,0) micr2=c(0.5,0), dofT=3
#############################
simuITDT<-function(x12,m1,m2,sigmaT,dofT,ny=10){
  micr11=m1[1]
  micr12=m1[2]
  micr21=m2[1]
  micr22=m2[2]
  x1=x12[1]
  x2=x12[2]
  d1=sqrt((x1-micr11)^2 + (x2-micr12)^2)
  d2=sqrt((x1-micr21)^2 + (x2-micr22)^2)
  Mutest = rep(abs(d1-d2),ny)
  Sigmatest<-diag(sigmaT,ny)
  #ysimu= Mutest+ rmvt(1, sigma=Sigmatest, df=1)
  ysimu=rmvt(1, delta=Mutest, sigma=Sigmatest, df=dofT)
  ysimu
}


################ loglikelihood x prior of yobs for the ITD model, 
# needed to run the MH algo of Geyer's package mcmc
# below for a Student noise dof=dofT 
logunpostITDT<-function(yobs, x12, m1, m2, sigmaT, dofT){
  x1=x12[1]
  x2=x12[2]
  micr11=m1[1]
  micr12=m1[2]
  micr21=m2[1]
  micr22=m2[2]
  ny=length(yobs)
  
  d1=sqrt((x1-micr11)^2 + (x2-micr12)^2)
  d2=sqrt((x1-micr21)^2 + (x2-micr22)^2)
  Mut = rep(abs(d1-d2), ny)
  Sigmat=diag(sigmaT,ny)
  if ((x1< -2) | (x1>2) |(x2< -2) | (x2>2) )
    # attention delta = mode not mean
    return(-Inf) else return(dmvt(yobs, delta=Mut, sigma=Sigmat, df=dofT, log=T))}
#return(-Inf) else return(dmvt(yobs, mean=Mut, sigma=Sigmat, df=1, log=T))}
##################################

########################################################
# Numerical computation of the true posterior for a ITD model
# when the observation is y.
# thetaseq are the theta values at which the pdf is computed and that
# depends on the chosen prior,

# eg y<-ytargetITD, N_grid=500

# unnormalized likelihood for contour plots
# TODO: unified later with function for MH algo
LikeITDTun<-function(theta,yobs,sigmaT,dofT, m1, m2){
ny <- length(yobs)
itd <- abs(sqrt(sum((theta-m1)^2))-sqrt(sum((theta-m2)^2)))
(1+1/(dofT*sigmaT)*sum((yobs-itd)^2))^(-(dofT+ny)/2)
}

PostPdfITD<-function(N_grid, y){
  
  x1 <- seq(-2,2, length = N_grid)
  x2 <- seq(-2,2, length = N_grid)
  grid <- expand.grid(x1 = x1, x2 = x2)
  z <- apply(grid,1, LikeITDTun, yobs=y,sigmaT=sigmaT,dofT=dofT, m1=micr1, m2=micr2 )
  
  # Full posterior (L=2) unormalized, ok for contour plots
  full_df <- cbind(grid, z)
  # ITD_df
  list("postdf"= full_df)
}

########################################################################
## target observation D=10 for true location in (1.5,1)
# Microphones configuration and default values in the function
micr1=c(-0.5,0)
micr2=c(0.5,0)
# Noise 
sigmaT=0.01
dofT=3

ytargetITD<-simuITDT(c(1.5,1), micr1,micr2,sigmaT,dofT,ny=10)
#ytargetITDtemp<-ytargetITD
#ytargetITD<-simuITDT(c(1,1), micr1,micr2,sigmaT,dofT,ny=10)
#ytargetITD<-simuITDT(c(0.6,1), micr1,micr2,sigmaT,dofT,ny=10)
# output: dim(ytargetITDT) =  1 x 10


## Example of use and  plot to  check
## Mics  and source locations (5)
dfconfig<-data.frame(matrix(c(-0.5,0,0.5,0,1.5,1), 3,2, byrow=T))
###  True posterior unnormalized contours

N_grid=500
ITD_df<-PostPdfITD(N_grid,ytargetITD)$postdf

v <- ggplot(ITD_df) +  xlab("x") + ylab("y")+ theme(aspect.ratio = 1)+ xlim(-2, 2) + ylim(-2, 2) ; 
v + geom_contour( aes(x1, x2, z = z), color="blue", bins=7) + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3) + geom_point(data=dfconfig, mapping=aes(x=X1, y=X2), size=2.5, color="black") 



###########################################################################################################
# Simulation of a first training set for GLLIM 
###########################################################################################################
#
# training set for GLLIM, to learn a parametric representation of the posteriors
# simulate N  positions  x1 and  x2 uniformly in  [-2,2]^2 and the associated ITD vectors
#
# For GLLIM data and param should be D x N and L x N

N = 10^5
# simulate N source locations (x1,x2)  uniformly in [-2,2]^2 
thetasimu1ITD <- matrix(runif(2*N,-2,2), 2, N)
# Output: dim is   L x N (L=2) 

# simulate as many corresponding  observations in dimension 10 (ny) with
# default values (make sure they coincide with what you want)
# using apply is faster even with 8 cores ??
system.time(ysimu1ITD<-apply(thetasimu1ITD, 2, simuITDT, m1=micr1,m2=micr2,sigmaT=sigmaT,dofT=dofT))
# output: dim is D x N  (D=10)
#user  system elapsed 
#39.210   0.521  39.884 

##########################################################################################################
# a second training samples for ABC  M1=10^6
##########################################################################################################
#
M1=10^6
thetasimu2ITD<-matrix( runif(2*M1,-2,2), 2, M1)
# simulate as many corresponding ITDs in dimension 10 (ny)
system.time(ysimu2ITD<-apply(thetasimu2ITD,2,simuITDT, m1=micr1,m2=micr2,sigmaT=sigmaT,dofT=dofT))
#user  system elapsed 
#390.520   3.480 395.033 

###################################################################################
#### Metropolis Hasting algorithm
#### Simulation of a sample of the "true" posterior approximated using a MH algorithm
###################################################################################
# Update February 2023:before that date the metrop function for the MH algorithm was not run correctly. Bug found by 
# Henrik Häggström: in the metrop function, blen has to be set to 1, which
# is the default value. Check details in help("metrop")
# scale maybe changed to 1 too, the higher "scale" the more moving is the chain, to be set
# in conjunction with nspac to compensate higher rejection. If scale is too small eg 0.1
# the chain has more trouble exploring oll the branches. 

# For burnin: first run of metrop
<<<<<<< HEAD
outmetropITDT<-metrop(logunpostITDTmix, c(0,0), nbatch=10^5, blen=1, yobs=ytargetITD, m1=micr1, m2=micr2, m1p=micr1p, m2p=micr2p,
=======
outmetropITDT<-metrop(logunpostITDT, c(0,0), nbatch=10^5, blen=1, yobs=ytargetITD, m1=micr1, m2=micr2, 
>>>>>>> florence
                      sigmaT=sigmaT,dofT=dofT, scale=1)

# sample after burnin:  nbatch =1000 produces 1000 points but not necessary different due to the rejection scheme
# to get 1000 different values in the sample, set nspac to a larger value, eg 100 or 1000 (more time consuming)
# Remark: if nspac=1 or 10 the number of unique values in the sample is likely to be less than 1000
# 
<<<<<<< HEAD
outITDT <- metrop(outmetropITDT, nspac=1000, blen = 1, nbatch = 1000, yobs=ytargetITD, m1=micr1, m2=micr2,m1p=micr1p, m2p=micr2p, sigmaT=sigmaT, dofT=dofT, scale=1 )
=======
outITDT <- metrop(outmetropITDT, nspac=1000, blen = 1, nbatch = 1000, yobs=ytargetITD, m1=micr1, m2=micr2, sigmaT=sigmaT, dofT=dofT, scale=1 )
>>>>>>> florence

############## The MH generated sample:
MHvalITD<-outITDT$batch


# Just a plot to check
# plot cex for the size of points, for a square plotpar(pty="s")
# plot(MHvalITD,xlim=c(-2,2),ylim=c(-2,2), pch=20, cex=0.5) 
dfMHITD<-data.frame(MHvalITD)
m <- ggplot(dfMHITD, aes(x = X1, y = X2)) +
  theme(aspect.ratio=1) +
  xlim(-2, 2) +
  ylim(-2, 2)
m<- m +  geom_point(data=dfMHITD, size=1, color=rgb(.1,0,.9,alpha = 0.5), shape=1) + xlab("x") + ylab("y")
m<- m + geom_point(data=dfconfig, size=2.5, color="black") + geom_abline(slope=0 , intercept=0, linetype="dashed", size=.3)
m


##################################################################################
# To save the data for later use: create a file dataITD.Rdata
save(list = c("ytargetITD","ysimu1ITD","thetasimu1ITD","ysimu2ITD","thetasimu2ITD","MHvalITD"), file = "dataITDdf3.RData")
##################################################################################



