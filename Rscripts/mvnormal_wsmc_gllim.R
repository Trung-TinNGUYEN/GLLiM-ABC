library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()

set.seed(11)

doRun <- FALSE
max_time <- 30*60
d <- 2
target <- get_multivariate_normal(d)
target$parameters$tau <- 5
nobservations <- 100
nparticles <- 2048
p <- 1
prefix <- ""

obsfile <- paste0(prefix, "mvnormaldata.d", d, ".n", nobservations, ".RData")
load(obsfile)

# function to simulate data
target$simulate <- function(theta){
  return(target$robservation(nobservations, theta, target$parameters))
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
# obs = observcations is 2 x 100
## Gllim model used
modg<-modgllimNLMIIDFullK1N105
# test with modg<-modgllimNLMIIDFullK1N102

#### step1 : compute posteriors parameters for y and z (mixtures)
# notational shortcut: gllim direct parameters
Aa<-modg$A
Sigmaa<-modg$Sigma
ca<-modg$c
ba<-modg$b
covstarR<-modg$covstarR
invGammaa<-modg$invGamma
invSigmaa<-modg$invSigma
D= dim(ba)[1]
L=dim(ca)[1]
# K<-1


# y is a DR-dim vector turned into a D x R matrix
#yDR<-matrix(y, ncol=R)
# yDR = obs

postcovy<-covstarR

# for K=1, the trace part is 0 in W2
# and the b part in the means cancels also because they are the same

# FULL case without the constant part...
# L x 1  (y)
postmeany<-t(t(covstarR))%*%t(Aa)%*%invSigmaa%*%t(t(rowSums(obs)))

# W2
dgllim<-function(z){
  # L x 1
  postmeanz<-t(t(covstarR))%*%t(Aa)%*%invSigmaa%*%t(t(rowSums(z)))
  
 return(sum((postmeany - postmeanz)^2))
}

# L2
#first constant part in L2 : can be made simpler...
postmean0<-rep(0,L)
L2cst<-2*L2scal2normal(postmean0,postmean0,postcovy,postcovy);

dgllimL2<-function(z){
  # L x 1
  postmeanz<-t(t(covstarR))%*%t(Aa)%*%invSigmaa%*%t(t(rowSums(z)))
  L2dist <- L2cst-2*L2scal2normal(postmeanz[,1],postmeany[,1],postcovy,postcovy)
  return(L2dist)
}
  
# distance between summary
summary_obs <- rowMeans(obs)
dsummary <- function(z){
  summary_z <- rowMeans(z)
  return(mean(abs(summary_z - summary_obs)))
}

# common algorithmic parameters
param_algo <- list(nthetas = nparticles, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)

filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".euclidean.RData")
results <- wsmc(deuclidean, target, param_algo, maxsimulation = 10^6, savefile = filename)
load(filename)

filenameWABC <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".wasserstein.RData")
results <- wsmc(wdistance, target, param_algo, maxsimulation = 10^6, savefile = filenameWABC)
load(filenameWABC)

filenameG <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".gllim.RData")
results <- wsmc(dgllim, target, param_algo, maxsimulation = 10^6, savefile = filenameG)
load(filenameG)

filenameG100 <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".gllim100.RData")
results <- wsmc(dgllim, target, param_algo, maxsimulation = 10^6, savefile = filenameG)
load(filenameG100)


# L2
filenameGL2 <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".gllimL2.RData")
results <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^6, savefile = filenameGL2)
load(filenameGL2)
# bis
filenameGL2b <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".gllimL2.RData")
results <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^6, savefile = filenameGL2b)
load(filenameGL2b)

filenameGL2100 <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".gllimL2100.RData")
results <- wsmc(dgllimL2, target, param_algo, maxsimulation = 10^6, savefile = filenameGL2)
load(filenameGL2100)


filenameS <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".summary.RData")
results <- wsmc(dsummary, target, param_algo, maxsimulation = 10^6, savefile = filenameS)
load(filenameS)



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
