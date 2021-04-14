#############################################################################################
####  Real data Planetary Science Example D=10 geometries, L=4 parameters
####  goal:  (see details in the paper)
##############################################################################################
#### Functions specific to this example, to create data sets, computation of the 
#### true posterior, all saved in a dataHapke.Rdata file at the end
#### 

#################################################################################################
###
### Step 1: Setting: target y (observations to be inverted, dimension is D) ;
###         Data sets simulation, one for GLLiM, one for ABC if wanted
###          are read from files provided by running the Hapke model
#################################################################################################
# training set for GLLIM, to learn a parametric representation of the posteriors
# The Hapke model + gaussian noise (see paper)
# Parameters: N=10^5 and L=4
thetasimu1H005<-as.matrix(read.csv2("../data-et-al4examples/Hapke005/Xgllim.csv", sep=",", header=FALSE))
# this command is reading "strings" and not numerics so that > xtrain[1,] looks like
#V1          V2          V3 
#"0.21357" "1.153e-05"   "0.85163" 
# therefore the following command:
thetasimu1H005<-apply(thetasimu1H005,2,as.numeric)
# output is N x 4

# geometries
ysimu1H005<-as.matrix(read.csv2("../data-et-al4examples/Hapke005/Ygllim.csv", sep=",", header=FALSE))
ysimu1H005<-apply(ysimu1H005,2,as.numeric)
# N x 10

# For GLLIM data and param should be D x N and L x N
thetasimu1H005<-t(thetasimu1H005)
ysimu1H005<-t(ysimu1H005)

#################################################
# Another training set for ABC 
# The Hapke model + gaussian noise (see paper)
# Parameters: N=10^5 and L=4

thetasimu2H005<-as.matrix(read.csv2("../data-et-al4examples/Hapke005/Xabc.csv", sep=",", header=FALSE))
thetasimu2H005<-apply(thetasimu2H005,2,as.numeric)
ysimu2H005<-as.matrix(read.csv2("../data-et-al4examples/Hapke005/Yabc.csv", sep=",", header=FALSE))
ysimu2H005<-apply(ysimu2H005,2,as.numeric)
thetasimu2H005<-t(thetasimu2H005)
ysimu2H005<-t(ysimu2H005)

################################################
# observations to be inverted
# 1) real observation : ytargetreal
# 2) simulated observation with highest correlation with ytargetreal: ytargetH005

####### 1)  real observation to be inverted : 1 x D 
ytargetreal<-as.matrix(read.csv2("../data-et-al4examples/Hapke005/Obs.csv", sep=",", header=FALSE))
ytargetreal<-apply(ytargetreal,2,as.numeric)
ytargetreal<-t(ytargetreal)

###### Compute Correlation or MSE (not R2) to find the vector in the simulated data set
# that is the closest to the observed vector
###
#
rsq <- function (x, y) cor(x, y) ^ 2
#
mse<-function(x,y){
  sum((x-y)^2)}

#yobs dim 1 x 10
#ysimu dim N x 10 
# find the index of the closest ysimu
MatchObs<-function(yobs, ysimu){
  nsimu<-dim(ysimu)[1]
  rsqlist<-NULL
  for (i in 1:nsimu){
    rsqlist<-c(rsqlist, rsq(yobs[1,], ysimu[i,]))
    # for the MSE change into
    #rsqlist<-c(rsqlist, mse(yobs[1,], ysimu[i,]))
  }
  index<-which.max(rsqlist)
  #index<-which.min(rsqlist)
  list(mind=index, allmse=rsqlist)
}

# example of use:  MatchObs(ytargetreal, t(ysimu2H005))
toto<-MatchObs(ytargetreal, t(ysimu2H005))
toto$mind 
#[1] 57702

# to look at the more correlated, 2nd most correlated etc...
# ATTENTION order is reversed if MSE is used
# titi<-rank(toto$allmse)
#seq(1:10^5)[titi==10^5]
#[1] 57702
#> seq(1:10^5)[titi==(10^5-1)]
#[1] 90006
#> seq(1:10^5)[titi==(10^5-2)]
#[1] 4108

# Plot to check real observation and selected signal: 1 x D
ts.plot(ts(ysimu2H005[,57702]),ts(ytargetreal[1,]), xlab="Geometry", col=c(2,1), lty=c(1,2))

##### 2) Highest correlated : 1 x D
ytargetH005<-t(ysimu2H005[,57702])

##################################################################################
# To save the data for later use: create a file dataHapke.Rdata
save(list = c("ytargetreal", "ytargetH005","ysimu1H005","thetasimu1H005","ysimu2H005","thetasimu2H005"), file = "dataHapke.RData")
##################################################################################

