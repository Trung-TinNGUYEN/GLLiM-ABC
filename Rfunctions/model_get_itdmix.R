# The observations are  D x R from a ITD mixture model
#'@rdname get_itdmix
#'@title multiple Hyperboloid model
#'@description This function returns a list representing
#' a ITD mixture model with length 10 # specified by the user as an argument.
#'@return The list contains:
#' rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), 
#' robservation (create synthetic data sets),
#' parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_itdmix <- function(ny){
  # unif on [-2,2]^2
  rprior <- function(nparticles,parameters){
    return(matrix(runif(nparticles*2, min = -2, max = 2), ncol = 2))
  }
  # evaluate the log-density of the prior, for each particle
  dprior <- function(thetaparticles,parameters){
    densities <- rep(0, nrow(thetaparticles))
    for (i in 1:nrow(thetaparticles)){
      if (any(thetaparticles[i,] > 2) || any(thetaparticles[i,] < -2)){
        densities[i] <- -Inf
      }
    }
    return(densities)
  }
  
  ##################################################################################
  # To simulate a ny-dim ITD observation from source in x12 (One component in the mixture)
  simuITDT<-function(x12,m1,m2,sigmaT,dofT,ny){
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
  ################################################################################
  # To simulate a 2 component ITD Mixture: 1 vector ysimu from source x12
  #parameters: m1,m2,m1p,m2p,sigmaT,dofT
  robservation<-function(theta,parameters,ny){
    m1<-c(-0.5,0)
    m2<-c(0.5, 0)
    m1p<-c(0,-0.5)
    m2p<-c(0, 0.5)
    sigmaT<-0.01
    dofT<-3
    usimu<-runif(1,0,1)
    if(runif(1,0,1)>0.5)
      ysimu<-simuITDT(theta,m1,m2,sigmaT,dofT, ny)
    else ysimu<-simuITDT(theta,m1p,m2p,sigmaT,dofT,ny)
   return(ysimu)
  }

  # parameters <- list()
  #
  parameters <- list()
  model <- list(rprior = rprior,
                dprior = dprior,
                robservation = robservation,
                parameter_names = c("Theta1", "Theta2"),
                parameters = parameters,
                thetadim = 2, ydim = 10)
  return(model)
}




