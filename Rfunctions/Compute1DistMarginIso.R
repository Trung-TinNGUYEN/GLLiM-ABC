Compute1DistMarginIso<-function(modg,y,z,imarg=0,wthr=0){
### in the L>1 case, it is faster to compute distances between posterior margins
### and then to apply the L=1 case!!!
### imarg is the list of index margins to be computed ie c(1,2) 
### or c(1,2, L)  by default if imarg set to 0
## wtr= weigh trheshold eg 0.001 we keep the weights higher
## modg is the output of GllimFit, contains the gllim model computed by xllim and other quantites that do not
## depend on y or z
# modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
## sumstaty includes quantities, ie posterior mixture, posterior expectation,
## posterior variance for the target y. They are provided by CompuSumStat function
# sumstaty<-CompSumstat(modg,y,wthr=0)

# z is D x1
### Faster  although quantitie for y are recomputed at each new z, 
### probably this is negligeable wrt the time spent in computing distances...
 
# modg: xllim direct model ($mod) + Inv Gammak, Sigmakstar (precomputed)
    
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
    
    # compute inverse model for y and z : a D x 2 matrix
    matyz<-cbind(y,z)
    # xllim inverse model
    gllimpredyz<-gllim_inverse_map(matyz,modeleg)
    # post mixture weights: dim is 2 x K 
    postweightyz<-gllimpredyz$alpha
    Piy=t(postweightyz)[,1]
    leftky<-seq(1,K)[Piy>wthr]
    dimay<-sum(Piy>wthr)
    Piy<-Piy[Piy>wthr]
    
    Piz=t(postweightyz)[,2]
    leftkz<-seq(1,K)[Piz>wthr]
    dimaz<-sum(Piz>wthr)
    Piz<-Piz[Piz>wthr]
    
    # quantites that depend on y (z)
    #Aks y + bks: 2 because y and z 
    postcovy<-array(0,c(L,L,dimay))
    postmeany<-matrix(0,L,dimay)
    postcovz<-array(0,c(L,L,dimaz))
    postmeanz<-matrix(0,L,dimaz)
    
    
    # cov dep pas de y or z
    postcovy[,,1:dimay]=covstara[,,leftky]
    postcovz[,,1:dimaz]=covstara[,,leftkz]
    # dim L x K : the mixtures means
    # for y 
    Asybs<-NULL
    for (k in leftky){
      # ISOTROPIC case
      Asybs<-cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(matyz[,1]-ba[,k]) +invGammaa[,,k]%*%ca[,k]) )
    }
    # L x dimay  (y)
    postmeany=Asybs
    
    # for z 
    Asybs<-NULL
    for (k in leftkz){
      # ISOTROPIC case
      Asybs <- cbind(Asybs,covstara[,,k]%*%(t(Aa[,,k]/Sigmaa[1,1,k])%*%(matyz[,2]-ba[,k]) +invGammaa[,,k]%*%ca[,k])) 
    }
    # L x dimaz  ( z)
    postmeanz=Asybs

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
      postmeanymarg<-postmeany[im,]
      postcovymarg<-postcovy[im,im,]
      postmeanzmarg<-postmeanz[im,]
      postcovzmarg<-postcovz[im,im,]
      
      # trick for dim 1
      postcovymarg<-array(postcovymarg, dim = c(1,1,dimay))
      postcovzmarg<-array(postcovzmarg, dim = c(1,1,dimaz))
      postmeanymarg<-array(postmeanymarg, dim = c(1,dimay))
      postmeanzmarg<-array(postmeanzmarg, dim = c(1,dimaz))
      
      mixy<-list("Mu"=postmeanymarg, "S"=postcovymarg, "Pi"=Piy)
      mixz<-list("Mu"=postmeanzmarg, "S"=postcovzmarg, "Pi"=Piz)
      
      #### step2 : Compute MW2 and L2 between imarg posterior margin given y and post given z
      MW2dist <- c(MW2dist, was2mixL1(mixy,mixz))
      L2dist <- c(L2dist, L2normalL1(mixy,mixz))
    }
    # end for imarg
    c(MW2dist,L2dist)
  }
  
 