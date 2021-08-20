##############################################################################################
####  Example Bivariate Beta model, D=2, L=5, R=500 (reduced in the code to 100)
##############################################################################################
#### Functions specific to this example, to create data sets, computation of the 
#### true posterior, all saved in a dataBBM.Rdata file in the end
#### Same principle/structure as previous example eg SMA2.R etc.
#### code modified from Hien's to get the data as DxiR x Nsimu matrix
##############################################################################################


# iR=500  D=2 and L=5 nsimu=10^5 (=nb_iters)
# iR is n the number of iid 2-dim simulation for the same theta

BBM_dataset_simu <- function(iR, nsimu) {
  # Initialize parameter values
  rand_Theta <- array(NA, 5)
  #simu_data <- array(NA, c(nsimu,iR,2))
  # stored_Theta <- array(NA, c(nsimu,5))
  simu_data <- matrix(NA, nsimu, 2*iR)  # here D=2
  stored_Theta <- matrix(NA, nsimu, 5) # here L=5
  
  for (k in 1:nsimu) {
    rand_Theta <- runif(5, min = 0, max = 5)
    U1 <- rgamma(iR,rand_Theta[1],1)
    U2 <- rgamma(iR,rand_Theta[2],1)
    U3 <- rgamma(iR,rand_Theta[3],1)
    U4 <- rgamma(iR,rand_Theta[4],1)
    U5 <- rgamma(iR,rand_Theta[5],1)
    V1 <- (U1+U3)/(U5+U4)
    V2 <- (U2+U4)/(U5+U3)
    Z1 <- V1/(V1+1)
    Z2 <- V2/(V2+1)
    simu_data[k,] <- as.vector(t(cbind(Z1,Z2)))
    # 
    # for (i in 1:nrow(true_data)) {
    #     simu_data[k,i,1] <- rbeta(1, rand_Theta[1]+rand_Theta[3], 
    #                               rand_Theta[4]+rand_Theta[5])
    #     simu_data[k,i,2] <- rbeta(1, rand_Theta[2]+rand_Theta[4], 
    #                               rand_Theta[3]+rand_Theta[5])            
    # }
    stored_Theta[k,] <- rand_Theta
  }
  
  list(Data = simu_data, Theta = stored_Theta)
}

# True param (1,1,1,1,1) rgmma(iR,theta,1)
iR=500
U1 <- rgamma(iR,1,1)
U2 <- rgamma(iR,1,1)
U3 <- rgamma(iR,1,1)
U4 <- rgamma(iR,1,1)
U5 <- rgamma(iR,1,1)
V1 <- (U1+U3)/(U5+U4)
V2 <- (U2+U4)/(U5+U3)
Z1 <- V1/(V1+1)
Z2 <- V2/(V2+1)
bbmobs<-as.vector(t(cbind(Z1,Z2)))  # 1x 1000

###################################### 
# simulation of 10^5 BBM learning set
Nsimu = 10^5  # (N)
learnset<-BBM_dataset_simu(iR=500,Nsimu)
bbmdata<-learnset$Data  # N x 2*iR (Nx 1000)
bbmparams<-learnset$Theta  # N x L (N x5)
rm(learnset)
# transpose for xllim
bbmdata<-t(bbmdata) # 1000 x 10^5
bbmparams<-t(bbmparams) # 5 x 10^5
bbmobs<-as.matrix(bbmobs)  # 1000 x 1

##################################################################################
# To save the data for later use: create a file dataBBM.Rdata
save(list = c("bbmobs","bbmdata","bbmparams"), file = "dataBBM.RData")
##################################################################################

