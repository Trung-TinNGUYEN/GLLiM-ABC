ComputeDistMarginIso<-function(modg,y,matz,imarg=0,wthr=0){
#################################################################################
# Case L> 1 and ALL distances (both MW2 and L2) from the y posterior to each of the z posterior for
# z in matz are computed    
# Compute both the MW2 and L2 distances between the GLLiM posterior for y and 
# the Gllim posterior for z 
# For faster computation when K is large, the mixture components with very small 
# weights can be removed as they are unlikely to impact the distance values
#
# Particularity: 
# this function computes distances between posterior MARGINS so that 
# the L=1 case can be applied which is much faster
#################################################################################
### Input:
#  - modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
#  - y and matz are resp D x 1 (target) and D x M (ABC sample) 
#  - imarg is the list of index margins to be computed ie c(1,2) 
#     or c(1,2,..., L)  by default if imarg set to 0
#  - wthr= weight threshold eg 0.001 recommended: we keep the mixture weights 
#    higher than wthr. Default is 0 meaning no threshold
### Output: 
# a matrix 2I x M  with the MW2 and L2 distance values for the chosen
#  margins (if there are I of them), eg. if all margins are computed: 2L x M where
# the first L rows contaim the M MW2 distances for margin 1 to L from z to y and
# last L rows contain the  L2 distances
#################################################################################
# Remark: although quantities for y are recomputed at each new z, this version is 
# faster! probably this is negligeable wrt the time spent in computing distances.
# ALSO the version for a matrix of z M x D is faster than a foreach of
# the 1 distance at a time version 
  
  #### step1 : compute posteriors parameters for y and z (mixtures)
  # notational shortcut: gllim direct parameters
  modeleg<-modg$mod
  #Gammaa<-modeleg$Gamma
  Aa<-modeleg$A
  Sigmaa<-modeleg$Sigma
  ca<-modeleg$c
  ba<-modeleg$b
  covstara<-modg$covstar
  invGammaa<-modg$invGamma
  
  L=dim(ca)[1]
  if(imarg==0) imarg<-seq(1,L)
  
  # compute inverse model for y and z : a D x (M+1) matrix
  
  matyz<-cbind(y,matz)
  # xllim inverse model
  # post mixture weights: dim is M+1 x K 
  postweightyz<-gllim_inverse_map(matyz,modeleg)$alpha
  dimMplus1<-dim(postweightyz)[1]
  
  #### the y post mixture needs to be computed separtly
  Piy=t(postweightyz)[,1]
  leftky<-seq(1,K)[Piy>wthr]
  dimay<-sum(Piy>wthr)
  Piy<-Piy[Piy>wthr]
  
  postcovy<-array(0,c(L,L,dimay))
  postmeany<-matrix(0,L,dimay)
  
  postcovy[,,1:dimay]=covstara[,,leftky]
  
  Asybs<-NULL
  for (k in leftky){
    # ISOTROPIC case
    Asybs<-cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(matyz[,1]-ba[,k]) +invGammaa[,,k]%*%ca[,k]) )
  }
  # L x dimay  (y)
  postmeany=Asybs
  
  # problem with networkflow when K=1
  # to avoid we duplicate the mixture in this case
  if (dimay==1) {
    postmeany<-cbind(postmeany,postmeany)
    temp=array(0,c(L,L,2))
    temp[,,1]<-postcovy
    temp[,,2]<-postcovy
    postcovy<-temp
    Piy=c(0.5,0.5)
    dimay<-2
  }
  ### end case of y 
  
  ### for all other z 
  listdist<-NULL
  
  for (ny in 2:dimMplus1){
    Piz=t(postweightyz)[,ny]
    leftkz<-seq(1,K)[Piz>wthr]
    dimaz<-sum(Piz>wthr)
    Piz<-Piz[Piz>wthr]
    
    postcovz<-array(0,c(L,L,dimaz))
    postmeanz<-matrix(0,L,dimaz)
    
    postcovz[,,1:dimaz]=covstara[,,leftkz]
    
    Asybs<-NULL
    for (k in leftkz){
      # ISOTROPIC case
      Asybs<-cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(matyz[,ny]-ba[,k]) +invGammaa[,,k]%*%ca[,k]) )
    }
    # L x dimaz  (y)
    postmeanz=Asybs
    
    if (dimaz==1) {
      postmeanz<-cbind(postmeanz,postmeanz)
      temp=array(0,c(L,L,2))
      temp[,,1]<-postcovz
      temp[,,2]<-postcovz
      postcovz<-temp
      Piz=c(0.5,0.5)
      dimaz<-2
    } 
    
    ## here we compute the mixture for margin in imarg
    MW2dist<-NULL
    L2dist<-NULL
    for(im in imarg){
      postmeanzmarg<-postmeanz[im,]
      postcovzmarg<-postcovz[im,im,]
      postmeanymarg<-postmeany[im,]
      postcovymarg<-postcovy[im,im,]
      # trick for dim 1
      postcovzmarg<-array(postcovzmarg, dim = c(1,1,dimaz))
      postmeanzmarg<-array(postmeanzmarg, dim = c(1,dimaz))
      postcovymarg<-array(postcovymarg, dim = c(1,1,dimay))
      postmeanymarg<-array(postmeanymarg, dim = c(1,dimay))
      
      mixy<-list("Mu"=postmeanymarg, "S"=postcovymarg, "Pi"=Piy)
      mixz<-list("Mu"=postmeanzmarg, "S"=postcovzmarg, "Pi"=Piz)
      
      #### step2 : Compute MW2 and L2 between imarg posterior margin given y and post given z
      MW2dist <- c(MW2dist, was2mixL1(mixy,mixz))
      L2dist <- c(L2dist, L2normalL1(mixy,mixz))
    }# end for imarg
    
    # list dist is 2L x M if all margins are computed
    listdist<-cbind(listdist, c(MW2dist,L2dist)) 
  } # end for ny
  
  return(listdist)
}
