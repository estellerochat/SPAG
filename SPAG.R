# Univariate SPAG
spag<-function(gene,env){
#Inputs:
  # for n individuals
  # gene is a dataframe with n rows and 3 columns
  #one row per individuals
  #columns: x, y, and the presence/absence of the genotype (or allele) of interest (coded: 0,1, or NA)
  # env is the raster layer of the environmental predictor
#Output: 
  # raster layer with the same extent and spatial resolution as env, 
  # showing the probability to find the genotype of interest
  
  ENV=raster::extract(env,cbind(gene$x,gene$y))
  gene=gene[,-which(colnames(gene) %in% c("x","y"))]
  me=glm(gene~ENV,family='binomial',na.action=na.omit)
  coeff=t(me$coefficients)
  PROB<-calc(env,fun=function(x){exp(coeff[1]+x %*% coeff[-1])/(1+exp(coeff[1]+x%*%coeff[-1]))})
  return(PROB)
}

# Intersection SPAG
I_spag<-function(gene,env){
#Inputs:
  # for n individuals and y genotypes(alleles) to be intersected
  # gene is a dataframe with n rows and (y+2) columns
  #one row per individuals
  # columns: x, y, and the presence/absence of the y genotypes (or allele) of interest (coded 0,1, or NA)
  # env is a raster stack with y layers for the y environmental predictors (in the same order as in the columns of gene)
#Output: 
  # raster layer with the same extent and spatial resolution as env, 
  # showing the probability to simultaneoulsy find all the genotypes provided in gene
  
  coordlist=cbind(gene$x,gene$y) #list of coordinates of the individuals
  colnames(coordlist)=c("x","y")
  gene=gene[,-which(colnames(gene) %in% c("x","y"))]
  
  #Compute all univariate SPAGs
  UnivSPAGs=raster()
  for(i in 1:ncol(gene)){
    UnivSPAGs=stack(UnivSPAGs,spag(cbind(coordlist,gene[,i,drop=F]),env[[i]]))
  }
  names(UnivSPAGs)=colnames(gene)
  
  #Compute the intersection-SPAG
  PROBmulti=UnivSPAGs[[1]] 
  for(i in 2:ncol(gene)){
    #Intermediate result: compute a SPAG with the U
    PROBint=spag(cbind(coordlist,gene[,i,drop=F]),stack(env[[i]],PROBmulti[[i-1]]))
    PROBmulti=stack(PROBmulti,PROBint*PROBmulti[[i-1]])
  }
  return(PROBmulti[[length(names(PROBmulti))]])
}

# Union SPAG
U_spag<-function(gene,env){
#Inputs:
  # for n individuals and y genotypes(alleles) to compute the union
  # gene is a dataframe with n rows and (y+2) columns
  #one row per individuals
  #columns: x, y, and the presence/absence of the y genotypes (or allele) of interest (coded 0,1, or NA)
  # env is a raster stack with y layers for the y environmental predictors (in the same order as in the columns of gene)
#Output: 
  # raster layer with the same extent and spatial resolution as env, 
  # showing the probability to find at least one of the genotypes provided in gene
  
  coordlist=cbind(gene$x,gene$y) #list of coordinates of the individuals
  colnames(coordlist)=c("x","y")
  gene=gene[,-which(colnames(gene) %in% c("x","y"))]
  
  #Compute all univariate SPAGs
  UnivSPAGs=raster()
  for(i in 1:ncol(gene)){
    UnivSPAGs=stack(UnivSPAGs,spag(cbind(coordlist,gene[,i,drop=F]),env[[i]]))
  }
  names(UnivSPAGs)=colnames(gene)
  
  #Compute the union-SPAG
  PROBmulti=sum(UnivSPAGs,na.rm=T)
  l=2
  while(ncol(gene)>l){
    GROUPES=combn(c(1:ncol(gene)),l)
    for(c in (1:ncol(GROUPES))){
      gpe=GROUPES[,c]
      PROBmulti=PROBmulti+(-1)^(l-1)*I_spag(cbind(coordlist,gene[,gpe,drop=F]),env[[gpe]])
    }
    l=l+1
  }
  PROBmulti=PROBmulti+(-1)^(l-1)*I_spag(cbind(coordlist,gene),env)
  return(PROBmulti)
}

# K-Percentage SPAG
K_spag<-function(gene,env,K){
#Inputs:
  # for n individuals and y genotypes(alleles) to compute the union
  # gene is a dataframe with n rows and (y+2) columns
  #one row per individuals
  #columns: x, y, and the presence/absence of the y genotypes (or allele) of interest (coded 0,1, or NA)
  # env is a raster stack with y layers for the y environmental predictors (in the same order as in the columns of gene)
  # K is the percentage to use for the K-percentage SPAG
#Output: 
  # raster layer with the same extent and spatial resolution as env, 
  # showing the probability to find at least one of the genotypes provided in gene
  
  coordlist=cbind(gene$x,gene$y) #list of coordinates of the individuals
  colnames(coordlist)=c("x","y")
  gene=gene[,-which(colnames(gene) %in% c("x","y"))]
  N=ncol(gene)
  NK=ceiling(K*N)
  
  #Compute the Kpercentage-SPAG
  PROBmulti=env[[1]]*0
  GROUPES=t(as.matrix(combn(1:N,NK)))
  gpefinal=cbind(rep(1,nrow(GROUPES)),GROUPES,matrix(0,nrow(GROUPES),N-NK))
  l=2
  while(nrow(GROUPES)>=l){
    BIGGROUPES=as.matrix(combn(1:nrow(GROUPES),l))
    for(c in (1:ncol(BIGGROUPES))){
      biggpecol=BIGGROUPES[,c]
      biggpe=unique(as.vector(GROUPES[biggpecol,]))
      res=matrix(0,1,N+1)
      res[1:(length(biggpe)+1)]=c((-1)^(l-1),biggpe)
      gpefinal=rbind(gpefinal,res)
    }
    l=l+1
  }
  gpefinal[,2:(N+1)]=t(apply(gpefinal[,2:(N+1)],1,sort))
  gpefinal=data.frame(gpefinal)
  colnames(gpefinal)=c("multi",paste("G",1:N,sep=""))
  sqltext=paste("SELECT sum(multi) as multi, ",paste("G",1:N,sep="",collapse=",")," FROM gpefinal GROUP BY ",paste("G",1:N,sep="",collapse=","),sep="")
  gpefinal=sqldf(sqltext)
  
  for(i in 1:nrow(gpefinal)){
    genelist=gpefinal[i,2:ncol(gpefinal)][gpefinal[i,2:ncol(gpefinal)]!=0]
    PROBmulti=PROBmulti+gpefinal$multi[i]*I_spag(cbind(coordlist,gene[,genelist,drop=F]),env[[genelist]])
  }
  return(PROBmulti) 
}
