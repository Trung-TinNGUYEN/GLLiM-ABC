#############################################################################
#
# Script to compute RMSE over 10 repetitions of the same experiment for the
# Bivariate Beta model  to get emprical means and MSE as in the
# ES paper 
# simulation of more observations (10) from theta=(1,1,1,1,1)
# requires function MA2 which is in MA2.R file
# REM: it's better to run it online and not as a function not to loose all
# computations in case it crashes...
# 
############################################################################


## Simulation of the 10 ytarget observations
iR=100
repet=10
truemat<-matrix(1,5,50)  # epsilon for 50 pts, 5 param
ymatobsBBM<-matrix(0,2*iR,repet)  # 200x 10
qymatobsBBM<-matrix(0,14,repet)  # 14 x 10

for (ni in 1:repet){
  U1 <- rgamma(iR,1,1)
  U2 <- rgamma(iR,1,1)
  U3 <- rgamma(iR,1,1)
  U4 <- rgamma(iR,1,1)
  U5 <- rgamma(iR,1,1)
  V1 <- (U1+U3)/(U5+U4)
  V2 <- (U2+U4)/(U5+U3)
  Z1 <- V1/(V1+1)
  Z2 <- V2/(V2+1)
  ymatobsBBM[,ni]<-as.vector(t(cbind(Z1,Z2)))  # 1x 200
  # add quantiles too
  qprob=seq(0,1,by=1/6)
  # for the reduced target: bbmobs[1:200,]
  #QuantTransf<-function(dataRDN,qprob){
   # Nsim<-dim(dataRDN)[[2]] # N ou 1
   # ql<-length(qprob)
   # quantmat<-matrix(0,2*ql,Nsim)
   # for (i in 1:Nsim){
    #  temp<-matrix(dataRDN[,i],nrow=2)  #length 2R turned into a 2 x R matrix
    #  quantmat[,i]<-c(quantile(temp[1,], qprob), quantile(temp[2,],qprob))
    #}
    #quantmat
  #}
  # ex 
  qymatobsBBM[,ni]<-QuantTransf(t(t(ymatobsBBM[,ni])),qprob) #14 x1

}

####################################################################
# Do not run with foreach, use a for loop
#system.time(listestim<- foreach (i=1:M, .combine=cbind) %dopar% { 
# 
# 11 methods to compare, 10 estimators: 5 means and 5 RMSE
####################################################################
listestimBBM<-array(0, c(11,10,repet))

### FOR over number of repetitions
for(i in 1:repet){
  #for(i in 3:repet){
  ytarget<-t(ymatobsBBM[,i])  # 200 x1 
  ytargetquant<-t(qymatobsBBM[,i]) # 14 x 1 
  # L2 and MW2 distance computation with th 10^-4 to save time
  distMW2L2<-PreCompDistFullParIID(modgllimBBMIIDFullK100n100,R=100,ytarget,bbmdata[1:200,],wthr = 10^(-4))
  # for 50 selected points 
  rejsampledist<-RejABCMW2L2(0.0005,distMW2L2[1,],distMW2L2[2,],bbmparams)
  
  #### comppute gllim posterior expectation and log variances in the IID case
  #ExpVarMatz<-ExpLogPostVarGllimFullIID(modgllimBBMIIDFullK100n100,bbmdata[1:200,],R=100) 
  # Do not change with target: ExpVarMatzBBMFullK100n100
  
  # t(t()) to keep the matrix format...
  ExpVarTarget<-ExpLogPostVarGllimFullIID(modgllimBBMIIDFullK100n100,t(ytarget),100) 
  #
  L=5
  espysimu<-ExpVarMatzBBMFullK100n100[,1:L]
  # espytarget is 1 x L 
  espytarget<-ExpVarTarget[1:L]
  #
  rejsampleE <-abc(target=espytarget, param=t(bbmparams), sumstat=espysimu, tol=.0005, method ="rejection")
  ### exp + log var 
  rejsampleEV <-abc(target=ExpVarTarget, param=t(bbmparams), sumstat=ExpVarMatzBBMFullK100n100, tol=.0005, method ="rejection")
  
  #mytf<-list(function(x){cbind(x,x^2,x^3,x^4)})
  saabc.BBMr<-semiauto.abc(obs=ytargetquant,param=t(bbmparams), sumstats=t(quantbbmdata),satr=mytf, tol=0.0005, overlap=TRUE, saprop=1, abcprop = 1,method="rejection", final.dens = TRUE)
  
  # gllim mixture for yobs
  ## 2D plot GLLiM mixture for target 50 points
  mixp<-ComputeGllimFullMixtParamIID(modgllimBBMIIDFullK100n100,R=100,ytarget,0)
  samplemixg<-simdataset(100,mixp$Pi, t(mixp$Mu),mixp$S)$X
  
  ################ Quantiles K=40
  mixparamQ<-ComputeGllimFullMixtParam(modgllimQuantFullK40,t(ytargetquant),0)
  samplemixgQ<-simdataset(50,mixparamQ$Pi, t(mixparamQ$Mu),mixparamQ$S)$X
  
  # Compute MW2 et L2 dist
  distMW2L2Q<-PreCompDistFullPar(modgllimQuantFullK40,t(ytargetquant),quantbbmdata,wthr = 0.0001)
  
  rejsampledistQ<-RejABCMW2L2(0.0005,distMW2L2Q[1,],distMW2L2Q[2,],bbmparams)
  # 
  # E et EV
  # espyQuant is of size N x L do not change with target
  #espyQ<-t(gllim_inverse_map(quantbbmdata,modgllimQuantFullK40$mod)$x_exp)
  
  # sum stat = posterior espectation for the target (obs or y test)
  # ytarget is 1 x D
  # espytarget is 1 x L 
  #
  # 1) inversion of the simulated observation
  espytargetQ=t(gllim_inverse_map(t(ytargetquant),modgllimQuantFullK40$mod)$x_exp)
  # 
  rejsampleEQ <-abc(target=espytargetQ, param=t(bbmparams), sumstat=espyQuant, tol=.0005, method ="rejection")
  
  ### EV  
  ##
  # simulated observation
  logvarytargetQ<-LogPostVarGllimFull(modgllimQuantFullK40,t(ytargetquant))
    
  ####### Gllim-EV
  # thetasimu2 is L x M
  # selected sample of size L x 100 is in rejsampleEV$unadj.values
  # esplogytarget  1 x 2L 
  # esplogysimu2   M x 2L
  esplogytargetQ<-c(espytargetQ, t(logvarytargetQ))
  esplogyQuant<-cbind(espyQuant, t(logvaryQuant))
  
  # threshold 0.0005---> 50 samples
  rejsampleEVQ <-abc(target=esplogytargetQ, param=t(bbmparams), sumstat=esplogyQuant, tol=.0005, method ="rejection")
  

    
  #############################################################
  # computing various estimators
  #############################################################

  estim<-matrix(0,11,10) # methods x estim 5 means, 5 RMSE
  #methods order: SA, GM, E, EV, L2, MW2, qGM, qE, qEV, qL2, qMW2
  # means L=5
  L=5
  estim[1,1:L]<-colMeans(saabc.BBMr$post.sample)
  estim[2,1:L]<-colMeans(samplemixg[1:50,])
  estim[3,1:L]<-colMeans(rejsampleE$unadj.values)
  estim[4,1:L]<-colMeans(rejsampleEV$unadj.values)
  estim[5,1:L]<-rowMeans(rejsampledist$L2postval)
  estim[6,1:L]<-rowMeans(rejsampledist$MW2postval)
  # with quantiles
  estim[7,1:L]<-colMeans(samplemixgQ)
  estim[8,1:L]<-colMeans(rejsampleEQ$unadj.values)
  estim[9,1:L]<-colMeans(rejsampleEVQ$unadj.values)
  estim[10,1:L]<-rowMeans(rejsampledistQ$L2postval)
  estim[11,1:L]<-rowMeans(rejsampledistQ$MW2postval)
  # RMSE
  estim[1,(L+1):(2*L)]<-sqrt(rowMeans((t(saabc.BBMr$post.sample) -truemat)^2))
  estim[2,(L+1):(2*L)]<-sqrt(rowMeans((t(samplemixg[1:50,]) -truemat)^2))
  estim[3,(L+1):(2*L)]<-sqrt(rowMeans((t(rejsampleE$unadj.values) -truemat)^2))
  estim[4,(L+1):(2*L)]<-sqrt(rowMeans((t(rejsampleEV$unadj.values) -truemat)^2))
  estim[5,(L+1):(2*L)]<-sqrt(rowMeans((rejsampledist$L2postval -truemat)^2))
  estim[6,(L+1):(2*L)]<-sqrt(rowMeans((rejsampledist$MW2postval -truemat)^2))
  # with quantiles
  estim[7,(L+1):(2*L)]<-sqrt(rowMeans((t(samplemixgQ) -truemat)^2))
  estim[8,(L+1):(2*L)]<-sqrt(rowMeans((t(rejsampleEQ$unadj.values) -truemat)^2))
  estim[9,(L+1):(2*L)]<-sqrt(rowMeans((t(rejsampleEVQ$unadj.values) -truemat)^2))
  estim[10,(L+1):(2*L)]<-sqrt(rowMeans((rejsampledistQ$L2postval -truemat)^2))
  estim[11,(L+1):(2*L)]<-sqrt(rowMeans((rejsampledistQ$MW2postval -truemat)^2))
  #
  
  listestimBBM[,,i]<-estim
} ##### end for loop

#) # end system.time
# listestim
#}  
#     user   system  elapsed 
# 2053.367   45.293  680.174 --> 11min pour M=2
# For repet=10
save(list = c("listestimBBM"), file = "listestim10BBM.RData")

#################################################
# Compute averages over repet
#################################################
avermse<-apply(listestimBBM, c(1,2), mean)
avermse3<-signif(avermse,4)

################################################
# make a latex table
################################################
#install.packages("lazyWeave")
library(lazyWeave)

orig_option <- getOption("lazyReportFormat")
options(lazyReportFormat="latex")
lazy.file.start(docClass="report", packages=c("pslatex", "palatino", "avant"),
                title="Report Name", author="Your Name")
#* Return the original option setting
options(lazyReportFormat=orig_option)

#titi<-lazy.table(toto[,,1])
#lazy.write(
 # lazy.file.start(),
 # titi,
  #lazy.file.end(),
 # OutFile="titi.tex")
#unlink("titi.tex")

#
tablatex<-lazy.table(avermse3, cborder=seq(0,10))
# REM : no borders and add then afterwards

#






