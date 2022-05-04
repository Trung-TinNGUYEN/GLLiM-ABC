# The observations are  D x R from a MA(2) model
#'@rdname get_ma2
#'@title MA(2) model
#'@description This function returns a list representing
#' a MA(2) model with length 150 # specified by the user as an argument.
#'@return The list contains:
#' rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), 
#' robservation (create synthetic data sets),
#' parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_ma2 <- function(len){
  target <- list()
  # unif on triangle:
  target$rprior <- function(nparticles,parameters){
    return(runif_in_triangle(nparticles, c(0,-1), c(-2,1), c(2,1)))
    #return(fast_rmvnorm(nparticles, rep(parameters$mu_0, dimension), diag(parameters$tau^2, dimension, dimension)))
    # particles <- matrix(nrow = nparticles, ncol = dimension)
    # for (id in 1:dimension){
    #   particles[,id] <- rnorm(nparticles, mean = parameters$mu_0, sd = parameters$tau)
    # }
    # return(particles)
  }
  # evaluate the log-density of the prior, for  particle thetas= (th1,th2)
  target$dprior <- function(thetas,parameters){
    # logdensities <- rep(0, nrow(thetas))
    theta1<-thetas[,1]
    theta2<-thetas[,2]
    if ((theta1<=-2) | (theta1>=2) | ((theta1+theta2)<=-1) | ((theta1-theta2) >=1))
      # attention delta = mode not mean
      logdens<- -Inf else{ 
      logdens<- 0}
     return(logdens)
    #return(fast_dmvnorm(thetas, rep(parameters$mu_0, dimension), diag(parameters$tau^2, dimension, dimension)))
    # for (id in 1:dimension){
    # logdensities <- logdensities + dnorm(thetas[,id],  mean = parameters$mu_0, sd = parameters$tau, log = TRUE)
    # }
    # return(logdensities)
  }
  
  # generate random variables used to compute a synthetic dataset
  target$generate_randomness <- function(nobservations){
    return(list())
    #return(fast_rmvnorm(nobservations, rep(0, dimension), diag(1, dimension, dimension)))
  }
  
  #S <- diag(1, dimension, dimension)
  #if (dimension > 1){
  #  for (i in 1:(dimension-1)){
   #   S[i,i+1] <- S[i+1,i] <- 0.5
  #  }
  #}
  target$parameters <- list()
  #target$parameters <- list(S = S, mu_0 = 0, tau = 10)
  ny<-len # 150
  
  # function to compute a dataset for each theta value
  target$robservation <- function(nobservations, theta, parameters, ...){
    obs <- matrix(0, nobservations,30)
    theta1<-theta[1]
    theta2<-theta[2]
    if (ny>3) {firstrow<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2, rep(0,ny-3))}
    else {firstrow<-c(theta1^2+theta2^2+1, theta1*theta2+theta1, theta2)}
    Sigmat<-toeplitz(firstrow)
      # D x R=nobs 
    obs<-matrix(rmvnorm(1,rep(0,ny),Sigmat), ncol=nobservations) 
    return(obs) # D=30 x nobs=5
    
    #return(t(fast_rmvnorm(nobservations, theta, parameters$S)))
    }
    
  
  #target$loglikelihood <- function(thetaparticles, observations, parameters){
  #  logdensities <- rep(0, nrow(thetaparticles))
   # for (i in 1:nrow(thetaparticles)){
   #   logdensities[i] <- sum(fast_dmvnorm(t(observations), thetaparticles[i,], parameters$S))
   # }
   # return(logdensities)
 # }
  target$thetadim <- 2
  target$ydim <- 30
  target$parameter_names <- paste0("X", 1:2)
  
  return(target)
}

