##############################################################################################
####  Example Sum of  2 MA(2) with opposite theta1, D=10, L=2
##############################################################################################
#### Functions specific to this example, to create data sets, computation of the 
#### true posterior, all saved in a dataSMA2.Rdata file in the end

#################################################################################################
###
### Step 1: Setting: target y (observations to be inverted, dimension is D) ;
###         Data sets simulation, one for GLLiM, one for ABC if wanted
###
#################################################################################################
# input: true parameters, simulated parameters etc...
# output: ytarget; (thetasimu1, ysimu1) size N; (thetasimu2, ysimu2) size M ;
#  true posterior computed numerically if possible

### To simulate a vector of dim ny from a sum of SMA2 a sort of MA1 not quite
# theta is a vector of length 2 c(val1, val2)
SMA2<-function(theta,ny=10){
  theta1<-theta[1]
  theta2<-theta[2]
  firstrow1<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2, rep(0,ny-3))
  Sigmat1<-toeplitz(firstrow1)
  firstrow2<-c(theta1^2+theta2^2+1, -theta1*theta2-theta1, theta2, rep(0,ny-3))
  Sigmat2<-toeplitz(firstrow2)
  #Sigmat2<-0
  Sigmat<-Sigmat1+Sigmat2
  ysim<-rmvnorm(1,rep(0,ny),Sigmat)
  ysim
}

## Likelihood of the sum x the prior uniform on the triangle
LikeSMA2<-function(theta,yobs){
  ny=length(yobs)
  theta1=theta[1]
  theta2=theta[2]
  firstrow1<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2, rep(0,ny-3))
  Sigmat1<-toeplitz(firstrow1)
  firstrow2<-c(theta1^2+theta2^2+1, -theta1*theta2-theta1, theta2, rep(0,ny-3))
  Sigmat2<-toeplitz(firstrow2)
  #Sigmat2<-0
  Sigmat<-Sigmat1+Sigmat2
  # to account for the prior here the triangle
  if ((theta1<=-2) | (theta1>=2) | ((theta1+theta2)<=-1) | ((theta1-theta2) >=1))
    # attention delta = mode not mean
    return(0) else return(dmvnorm(yobs,rep(0,ny),Sigmat))
}

## Likelihood for R iid observations yR of the SMA(2) 
## x the prior uniform on the triangle
## yobsR is of length RD (the R iid vectors are concatanated in one vector)
LikeSMA2R<-function(theta,yobsR,R){
  ny=length(yobsR)/R  #D
  theta1=theta[1]
  theta2=theta[2]
  firstrow1<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2, rep(0,ny-3))
  Sigmat1<-toeplitz(firstrow1)
  firstrow2<-c(theta1^2+theta2^2+1, -theta1*theta2-theta1, theta2, rep(0,ny-3))
  Sigmat2<-toeplitz(firstrow2)
  #Sigmat2<-0
  Sigmat<-Sigmat1+Sigmat2
  # to account for the prior here the triangle
  if ((theta1<=-2) | (theta1>=2) | ((theta1+theta2)<=-1) | ((theta1-theta2) >=1))
    # attention delta = mode not mean
    return(0) else{ 
      likeR<-1
      for(r in 1:R){
      likeR<-LikeR * dmvnorm(yobsR[(1+ny*(r-1)):(r*ny)],rep(0,ny),Sigmat)
       }
}}

########################################################
# Numerical computation of the true posterior for a SMA2 model
# when the observation is y.
# thetaseq are the theta values at which the pdf is computed and that
# depends on the chosen prior, eg the triangle

# plot  normalise of the TRUE POST
thetaseq1<-seq(-2,2,length=400)
#thetaseq2<-seq(-1,1,length=400)
thetaseq2<-seq(-1,1,length=200) # no bumps 

# eg y<-ytargetSMA2

PostPdfSMA2<-function(thetaseq1,thetaseq2, y){
  # numerical normalisation for plot of the true marginal posterior
  # unormalized for the joint post
  
  thetaN1<-length(thetaseq1)
  thetaN2<-length(thetaseq2)
  # unormalized posterior: Likelihood x triangle prior on a grid 
  grid <- expand.grid(x1 = thetaseq1, x2 = thetaseq2)
  z <- apply(grid,1, LikeSMA2, yobs=y)*10^4
  
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
  # SMA2_df
  list("marg1"=margpdf1,"marg2"=margpdf2,"postdf"= full_df)
}

## Example of use 
# SMA2_df<-PostPdfSMA2(thetaseq1,thetaseq2, ytargetSMA2)$postdf
#
### PLOT of contours of true posterior  for y=ytargetSMA2 (jan 2021)
#v + geom_contour(bins=50)
v <- ggplot(SMA2_df, aes(x1, x2, z = z)) + xlim(-2,2) + ylim(-1,1) +theme(aspect.ratio=.9)
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


### For true posterior marginals and joint in the IID case 
## for IID GLLiM R>1 here R=10
#y<-ytargetRSMA2[1,11:20]
#z2<-apply(grid,1, LikeSMA2, yobs=y)*10^10
#z<-z*z2
#y<-ytargetR3SMA2[1,21:30]
#z3<-apply(grid,1, LikeSMA2, yobs=y)*10^10
#z<-z*z3
#for(r in 1:7){
#  y<-ytargetR10SMA2[1,(21+r*10):(30+r*10)]
#  z3<-apply(grid,1, LikeSMA2,yobs=y)*10^10
#  z<-z*z3
#}

#yR is the concatenation of R SMA2 simulation eg. for R=10
#yR<-NULL
#for(r in 1:10){yR<-c(yR,SMA2(c(0.7,0.5), 10))}

PostPdfSMA2R<-function(thetaseq1,thetaseq2, yR,R){
  # numerical normalisation for plot of the true marginal posterior
  # unormalized for the joint post
  
  thetaN<-length(thetaseq1)
  # unormalized posterior: Likelihood x triangle prior on a grid 
  grid <- expand.grid(x1 = thetaseq1, x2 = thetaseq2)
  z <- apply(grid,1, LikeSMA2R, yobs=yR, R=R)*10^10
  
  ## Marginal posteriors
  Like1<-rep(0,thetaN)
  Like2<-rep(0,thetaN)
  for(i in 1:thetaN){
    Liketemp1<-NULL
    Liketemp2<-NULL
    for(j in 0:(thetaN-1)){
      Liketemp1<-c(Liketemp1,z[i+j*thetaN])
      Liketemp2<-c(Liketemp2,z[(j+1)+(i-1)*thetaN])
    }
    normc1<-integrate.xy(thetaseq2, Liketemp1)
    Like1[i]<-normc1
    normc2<-integrate.xy(thetaseq1, Liketemp2)
    Like2[i]<-normc2
  }
  margpdf1<-Like1/(integrate.xy(thetaseq1, Like1))
  margpdf2<-Like2/(integrate.xy(thetaseq2, Like2))
  
  # Full posterior (L=2) unormalized, ok for contour plots
  full_df <- cbind(grid, z)
  # SMA2_df
  list("marg1"=margpdf1,"marg2"=margpdf2,"postdf"= full_df)
}



########################################################################
## target observation D=10 for true theta1=0.7 and theta2=0.5
ytargetSMA2<-SMA2(c(0.7,0.5),ny=10)

## iid case R=10
ytargetSMA2R10<-NULL
for(r in 1:10){ytargetSMA2R<-c(ytargetSMA2R,SMA2(c(0.7,0.5), 10))}

########################################################################
# Simulation of a first training set for ABC and GLLIM 

N =10^5
# simulate N theta1 theta2  uniformly in  triangle
thetasimu1SMA2 <- runif_in_triangle(N, c(0,-1), c(-2,1), c(2,1));
# Output: dim is  N x L (L=2) 

# simulate as many corresponding  observations in dimension 10 (ny)
# using apply is faster even with 8 cores ??
system.time(ysimu1SMA2<-apply(thetasimu1SMA2, 1, SMA2))
# output: dim is D x N  (D=10)
#user  system elapsed 
#28.987   0.580  29.612 

# using doP
#numCores <- detectCores()
#registerDoParallel(numCores)
#system.time(ysimu1SMA2<- foreach (i=1:N, .combine=rbind) %dopar% {
#  SMA2(thetasimu1SMA2[i,], 10)
#})
#output is N x D
##user  system elapsed 
##54.551   2.300  57.015 

# For xllim data and param should be D x N and L x N
thetasimu1SMA2 <- t(thetasimu1SMA2)
# if doP is used then: ysimu1SMA2<-t(ysimu1SMA2)


########################################################################
# a second training samples for ABC...
M=10^5
thetasimu2SMA2 <- runif_in_triangle(M, c(0,-1), c(-2,1), c(2,1));
# output is M x L 
# simulate as many corresponding  observations in dimension 10 (ny)
ysimu2SMA2<-apply(thetasimu2SMA2, 1, SMA2)
# output: dim is D x M  (D=10)

thetasimu2SMA2 <- t(thetasimu2SMA2)

########################################################################
# a third training samples for ABC... BIG in memory
M1=10^6
thetasimu3SMA2 <- runif_in_triangle(M1, c(0,-1), c(-2,1), c(2,1));
# simulate as many corresponding  observations in dimension 10 (ny)
ysimu3SMA2<-apply(thetasimu3SMA2, 1, SMA2)
thetasimu3SMA2 <- t(thetasimu3SMA2)

##################################################################################
# To save the data for later use: create a file dataSMA2.Rdata
save(list = c("ytargetSMA2","ysimu1SMA2","thetasimu1SMA2","ysimu2SMA2","thetasimu2SMA2","ysimu3SMA2","thetasimu3SMA2"), file = "dataSMA2.RData")
##################################################################################

