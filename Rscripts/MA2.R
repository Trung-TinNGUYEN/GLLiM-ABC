##############################################################################################
#### 
#### Example of MA(2) D=150 , L=2 
####                   
##############################################################################################
#### Functions specific to this example, to create data sets, computation of the 
#### true posterior, all saved in a dataMA2.Rdata file in the end AND data3MA2.RData (with 10 AC for BSL)

#################################################################################################
###
### Step 1: Setting: target y (observations to be inverted, dimension is D) ;
###         Data set simulation, same for GLLiM and ABC 
###
#################################################################################################
# input: true parameters, simulated parameters etc...
# output: ytargetMA2150,ysimu2MA2150 etc... (thetasimu2MA2, ysimu3MA2150) size M
#  true posterior computed numerically if possible

#################################################
## To simulate a vector of dim ny from a MA2 
## theta is a vector of length 2 c(val1, val2)
## Remark: ny is at least 3 
#################################################
MA2<-function(theta,ny=150){
  theta1<-theta[1]
  theta2<-theta[2]
  if (ny>3) {firstrow<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2, rep(0,ny-3))}
  else {firstrow<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2)}
  Sigmat<-toeplitz(firstrow)
  ysim<-rmvnorm(1,rep(0,ny),Sigmat)
  ysim
}

###################################################
## Likelihood  x the prior uniform on the triangle
###################################################
LikeMA2<-function(theta,yobs){
  ny=length(yobs)
  theta1=theta[1]
  theta2=theta[2]
  if (ny>3) {firstrow<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2, rep(0,ny-3))}
  else {firstrow<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2)}
  Sigmat<-toeplitz(firstrow)
  # to account for the prior here the triangle
  if ((theta1<=-2) | (theta1>=2) | ((theta1+theta2)<=-1) | ((theta1-theta2) >=1))
    # attention delta = mode not mean
    return(0) else return(dmvnorm(yobs,rep(0,ny),Sigmat))
}

###########################################################################
## Likelihood for R iid observations yR of the MA(2) 
## x the prior uniform on the triangle
## yobsR is of length RD (the R iid vectors are concatanated in one vector)
###########################################################################
LikeMA2R<-function(theta,yobsR,R){
  ny=length(yobsR)/R  #D
  theta1=theta[1]
  theta2=theta[2]
  if (ny>3) {firstrow<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2, rep(0,ny-3))}
  else {firstrow<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2)}
  Sigmat<-toeplitz(firstrow)
  # to account for the prior here the triangle
  if ((theta1<=-2) | (theta1>=2) | ((theta1+theta2)<=-1) | ((theta1-theta2) >=1))
    # attention delta = mode not mean
    likeR<-0 else{ 
      likeR<-1
      for(r in 1:R){
        likeR<-likeR * dmvnorm(yobsR[(1+ny*(r-1)):(r*ny)],rep(0,ny),Sigmat)
      }
    } 
  likeR
}

######################################################################
# Numerical computation of the true posterior for a MA2 model
# when the observation is y.
# thetaseq are the theta values at which the pdf is computed and that
# depends on the chosen prior, eg the triangle
######################################################################

# plot  normalise of the TRUE POST
thetaseq1<-seq(-2,2,length=400)
#thetaseq2<-seq(-1,1,length=400)
thetaseq2<-seq(-1,1,length=200) # no bumps 

# eg y<-ytargetMA2

PostPdfMA2<-function(thetaseq1,thetaseq2, y){
  # numerical normalisation for plot of the true marginal posterior
  # unormalized for the joint post
  
  thetaN1<-length(thetaseq1)
  thetaN2<-length(thetaseq2)
  # unormalized posterior: Likelihood x triangle prior on a grid 
  grid <- expand.grid(x1 = thetaseq1, x2 = thetaseq2)
  z <- apply(grid,1, LikeMA2, yobs=y)*10^4
  
  ## Marginal posteriors
  Like1<-rep(0,thetaN1)
  Like2<-rep(0,thetaN2)
  for(i in 1:thetaN1){
    Liketemp1<-NULL
    # Liketemp2<-NULL
    for(j in 0:(thetaN2-1)){
      Liketemp1<-c(Liketemp1,z[i+j*thetaN1])
      # Liketemp2<-c(Liketemp2,z[(j+1)+(i-1)*thetaN])
      Like2[j+1]<-integrate.xy(thetaseq1, z[(1+j*thetaN1):((j+1)*thetaN1)])
    }
    Like1[i]<-integrate.xy(thetaseq2, Liketemp1)
    # Like1[i]<-normc1
    # normc2<-integrate.xy(thetaseq1, Liketemp2)
    # Like2[i]<-normc2
  }
  margpdf1<-Like1/(integrate.xy(thetaseq1, Like1))
  margpdf2<-Like2/(integrate.xy(thetaseq2, Like2))

  # Full posterior (L=2) unormalized, ok for contour plots
  full_df <- cbind(grid, z)
  # MA2_df
  list("marg1"=margpdf1,"marg2"=margpdf2,"postdf"= full_df)
}

#################
## Example of use 
## MA2_df<-PostPdfMA2(thetaseq1,thetaseq2, ytargetMA2)$postdf
##
##################
### PLOT of contours of true posterior  for y=ytargetMA2 
#v + geom_contour(bins=50)
v <- ggplot(MA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
v<- v + geom_vline(xintercept=1, linetype="dashed", size=.5) + geom_hline(yintercept=0.6, linetype="dashed", size=.5)
v<- v +  geom_vline(xintercept=-1, linetype="dashed", size=.5)
#v<- v + geom_point(data=dftrue, size=2.5, color="black")
v<-v + geom_contour(breaks=c(min(z), seq(1.24,  max(z), length.out = 5)), color="black")
v<- v + xlab("theta1") + ylab("theta2")
# to get rid of the grey background
v + theme_bw()
#v + geom_contour(breaks=c(  0.9, 1.1, 1.2, seq(1.23,  max(z), length.out = 10)))
#v + geom_contour(bins=20)
#v + geom_density_2d(size = .5, aes(color = ..level..)) 
#

###########################################################
### For true posterior marginals and joint in the IID case 
## for IID GLLiM R>1 here R=10
#y<-ytargetRMA2[1,11:20]
#z2<-apply(grid,1, LikeMA2, yobs=y)*10^10
#z<-z*z2
#y<-ytargetR3MA2[1,21:30]
#z3<-apply(grid,1, LikeMA2, yobs=y)*10^10
#z<-z*z3
#for(r in 1:7){
#  y<-ytargetR10MA2[1,(21+r*10):(30+r*10)]
#  z3<-apply(grid,1, LikeMA2,yobs=y)*10^10
#  z<-z*z3
#}

#yR is the concatenation of R MA2 simulation eg. for R=10
#yR<-NULL
#for(r in 1:10){yR<-c(yR,MA2(c(0.7,0.5), 10))}
###########################################################

PostPdfMA2R<-function(thetaseq1,thetaseq2, yR,R){
  # numerical normalisation for plot of the true marginal posterior
  # unormalized for the joint post
  
  thetaN1<-length(thetaseq1)
  thetaN2<-length(thetaseq2)
  # unormalized posterior: Likelihood x triangle prior on a grid 
  grid <- expand.grid(x1 = thetaseq1, x2 = thetaseq2)
  z <- apply(grid,1, LikeMA2R, yobs=yR, R=R)*10^10
  
  ## Marginal posteriors
  Like1<-rep(0,thetaN1)
  Like2<-rep(0,thetaN2)
  for(i in 1:thetaN1){
    Liketemp1<-NULL
    # Liketemp2<-NULL
    for(j in 0:(thetaN2-1)){
      Liketemp1<-c(Liketemp1,z[i+j*thetaN1])
      # Liketemp2<-c(Liketemp2,z[(j+1)+(i-1)*thetaN])
      Like2[j+1]<-integrate.xy(thetaseq1, z[(1+j*thetaN1):((j+1)*thetaN1)])
    }
    Like1[i]<-integrate.xy(thetaseq2, Liketemp1)
    # Like1[i]<-normc1
    # normc2<-integrate.xy(thetaseq1, Liketemp2)
    # Like2[i]<-normc2
  }
  margpdf1<-Like1/(integrate.xy(thetaseq1, Like1))
  margpdf2<-Like2/(integrate.xy(thetaseq2, Like2))
  
  # Full posterior (L=2) unormalized, ok for contour plots
  full_df <- cbind(grid, z)
  # MA2_df
  list("marg1"=margpdf1,"marg2"=margpdf2,"postdf"= full_df)
}

# Example of use
#postpdfMA2R<-PostPdfMA2R(thetaseq1,thetaseq2, ytargetR50MA2, 10)


########################################################################
##
## target observation D=150 for true theta1=0.6 and theta2=0.2
##
########################################################################
ytargetMA2150<-MA2(c(0.6,0.2),ny=150)


########################################################################
# Simulation of a first training set for ABC and GLLIM 
# NOT USED at first

N =10^5
# simulate N theta1 theta2  uniformly in  triangle
thetasimu1MA2 <- runif_in_triangle(N, c(0,-1), c(-2,1), c(2,1));
# Output: dim is  N x L (L=2) 

# simulate as many corresponding  observations in dimension 150 (ny)
# using apply is faster than doParallel even with 8 cores ??
system.time(ysimu1MA2<-apply(thetasimu1MA2, 1, MA2))
# output: dim is D x N  (D=150)
#user  system elapsed 
#28.987   0.580  29.612  

# For xllim data and param should be D x N and L x N
thetasimu1MA2 <- t(thetasimu1MA2)
# if doParallel is used then: ysimu1MA2<-t(ysimu1MA2)


########################################################################
##
## a second training samples for ABC... USED
##
########################################################################

M=10^5
thetasimu2MA2 <- runif_in_triangle(M, c(0,-1), c(-2,1), c(2,1));
# output is M x L 
# simulate as many corresponding  observations in dimension 150 (ny)
ysimu2MA2150<-apply(thetasimu2MA2, 1, MA2) # long
# output: dim is D x M  (D=150)

thetasimu2MA2 <- t(thetasimu2MA2)


########################################################################
# a Third training samples for ABC...
M3=10^6
thetasimu3MA2 <- runif_in_triangle(M3, c(0,-1), c(-2,1), c(2,1));
# output is M3 x L 
# simulate as many corresponding  observations in dimension 150 (ny)
ysimu3MA2150<-apply(thetasimu3MA2, 1, MA2) # takes more than 3h 
# output: dim is D x M3  (D=150)

thetasimu3MA2 <- t(thetasimu3MA2)

### 10 first lag Autocovariances in ysimu3ACF10MA2150
ysimu3ACF10MA2150<-apply(ysimu3MA2150, 2, fnSum_acf, Ka=10)
# dim is 10 x 10^6


########################################################################
# Test with autocovariances as summary statistics
########################################################################
# requires function fnSum_acf.R
# Here we use Ka= 10 autocovariances, ie 10 different lags
# target 
ytargetACF10MA2150<-fnSum_acf(ytargetMA2150,10)
# Learning and ABC set, from the same simulations as before
ysimu2ACF10MA2150<-apply(ysimu2MA2150, 2, fnSum_acf, Ka=10)
# dim is 10 x 10^5



##################################################################################
##
## To save the data for later use: create a file dataMA2.Rdata
save(list = c("ytargetMA2150","ysimu2MA2150","thetasimu2MA2"), file = "dataMA2.RData")
##################################################################################

##################################################################################
# To save the data for later use: create a file data3MA2.Rdata
save(list = c("ytargetMA2150","ysimu3MA2150","ysimu3ACF10MA2150","thetasimu3MA2"), file = "data3MA2.RData")
######################################################################################




