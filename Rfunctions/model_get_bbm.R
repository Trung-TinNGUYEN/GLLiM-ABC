# The observations are  D x R from a Bivariate Beta model
#'@rdname get_bbm
#'@title BBM model
#'@description This function returns a list representing
#' a BBM model 
#'@return The list contains:
#' rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), 
#' robservation (create synthetic data sets),
#' parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_bbm <- function(dimension){
  target <- list()
  # unif on [0,5]^5 
  target$rprior <- function(nparticles,parameters){
    # output is nparticles x L (L=5) 
    return(matrix(runif(nparticles*5, 0,5), ncol=5) )
    #return(fast_rmvnorm(nparticles, rep(parameters$mu_0, dimension), diag(parameters$tau^2, dimension, dimension)))
    # particles <- matrix(nrow = nparticles, ncol = dimension)
    # for (id in 1:dimension){
    #   particles[,id] <- rnorm(nparticles, mean = parameters$mu_0, sd = parameters$tau)
    # }
    # return(particles)
  }
  # evaluate the log-density of the prior, for  particle thetas= (th1,th2..th5)
  target$dprior <- function(thetas,parameters){
    if ((sum(thetas<=5) + sum(thetas>=0))==10 )
      # attention delta = mode not mean
      logdens<- 0 else{ 
        logdens<- -Inf}
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
 # ny<-len # 150
  
  # function to compute a dataset for each theta value : nobs x D =2
  target$robservation <- function(nobservations, theta, parameters, ...){
    obs <- matrix(0, nobservations,2)
        U1 <- rgamma(nobservations,theta[1],1)
        U2 <- rgamma(nobservations,theta[2],1)
        U3 <- rgamma(nobservations,theta[3],1)
        U4 <- rgamma(nobservations,theta[4],1)
        U5 <- rgamma(nobservations,theta[5],1)
        V1 <- (U1+U3)/(U5+U4)
        V2 <- (U2+U4)/(U5+U3)
        Z1 <- V1/(V1+1)
        Z2 <- V2/(V2+1)
        simu_data <- t(cbind(Z1,Z2))
      return(simu_data)   # D=2 x nobs=100
    }

  #target$loglikelihood <- function(thetaparticles, observations, parameters){
  #  logdensities <- rep(0, nrow(thetaparticles))
  # for (i in 1:nrow(thetaparticles)){
  #   logdensities[i] <- sum(fast_dmvnorm(t(observations), thetaparticles[i,], parameters$S))
  # }
  # return(logdensities)
  # }
  target$thetadim <- 5
  target$ydim <- 2
  target$parameter_names <- paste0("Theta", 1:5)
  
  return(target)
}

