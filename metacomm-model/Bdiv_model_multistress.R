#disp is the proportion of biomass that disperses in each time step - recommend using values between 0 and 0.001
#for Int_type choose the type of interactions: 
#"NoInt" - no interactions; "Comp" - competition; "Mixed" -competition, mutualism, and parasitsm; "Trophic" - tri trophic food web

Bdiv_model<-function(disp=0.01,Int_type="Comp",nSpecies=120,Stress = 1,Stress2=1){
  Tmax<-1000+1000*Stress #length of simulation 1000 is enough to reach equilibrium
  
  nCom_rows<-20 #number of rows in landscape matrix
  nCom_cols<-20 #number of columns in landscape matrix
  nCom<-nCom_cols*nCom_rows #total number of patches in landscape matrix
  
  #number of species in each trophic level
  if(Int_type=="Trophic"){
    nPrey<-nSpecies*0.5
    nHerb<-nSpecies*0.3
    nPred<-nSpecies*0.2
  } else{
    nPrey<-nSpecies
    nHerb<-0
    nPred<-0
  }
  
  #Identity vectors for each trophic level
  preyV<-1:nPrey
  herbV<-(nPrey+1):(nPrey+nHerb)
  pred<-(nSpecies-nPred+1):(nSpecies)
  
  #environmental gradients
  Env1<-matrix(seq(1,nCom_rows),nCom_rows,nCom_rows)
  Env2<-matrix(rep(seq(1,nCom_cols),each=nCom_cols),nCom_cols,nCom_cols)
  maxStress1<-0 #the amount of increase in environmental variable 1
  maxStress2<-0 #the amount of increase in environmental variable 2
  
  #Stress
  Stressor<-c(rep(0,1000),seq(0,Stress,length=(Stress*1000)))
  Stressor1<-rep(1,nCom)
  Stressor2<-rep(0,nCom)
  Stressor2[sample(nCom,nCom/2,replace=F)]<-Stress2
    
  #species environmental niches
  Opt1<-c(seq(1-10,max(Env1)+10,length=nPrey), seq(1,max(Env1),length=nHerb),seq(1,max(Env1),length=nPred))
  T_Norm<-apply(t(Opt1),2,dnorm,sd=50,x=seq(1,maxStress1+max(Env1)))*300
  A1<-(T_Norm-max(T_Norm))
  
  
  Opt2<-c(seq(1-10,max(Env1)+10,length=nPrey)[sample(nPrey,replace = F)], seq(1,max(Env2),length=nHerb)[sample(nHerb,replace = F)],seq(1,max(Env2),length=nPred)[sample(nPred,replace = F)])
  T_Norm<-apply(t(Opt2),2,dnorm,sd=50,x=seq(1,maxStress2+max(Env2)))*300
  A2<-(T_Norm-max(T_Norm))
  
  Stress_tolerance<-sample(seq(-0.01,0,length=1000),size = nSpecies, replace=T)
  Stress_tolerance2<-sample(seq(-0.01,0,length=1000),size = nSpecies, replace=T)
  
  traits<-cbind(Opt1,Opt2,Stress_tolerance,Stress_tolerance2)
  row.names(traits)<-1:nSpecies
  
  #interaction matricies####
  #competitive
  weight=1/nSpecies*10 #weight interaction strength
  
  if(Int_type=="Comp" | Int_type=="NoInt"){
    b11=-.15 #mean interspecific competition strength
    bdiag1=-.2 #intraspecific competition strength
    BB=b11*matrix(runif(nPrey*nPrey),nPrey,nPrey)
    diag(BB)<-bdiag1
    BB=weight*BB
    B<-BB
    if(Int_type=="NoInt"){B=diag(diag(BB))}
  } else {if(Int_type=="Mixed"){
    b11=-.15 #mean interspecific competition strength
    bdiag1=-.2 #intraspecific competition strength
    BB=b11*matrix(runif(nPrey*nPrey),nPrey,nPrey)
    diag(BB)<-bdiag1
    BB=weight*BB
    BI<-BB
    
    BB<-matrix(-1,nSpecies,nSpecies)
    int.n<-sum(BB[upper.tri(BB)])*-1
    BB[upper.tri(BB)][sample(int.n, replace=F)<=(0.35*int.n)]<-0.5
    BB[lower.tri(BB)][t(BB)[lower.tri(BB)]>0][sample(0.35*int.n, replace=F)<(0.10*int.n)]<-0.5
    B<-BB*-BI
  } else {if(Int_type == "Trophic"){
    b11=-0.1
    b12=-0.3
    b21=0.1
    b23=-.1
    b32=.08
    bdiag1=-.2
    bdiag2=-.15
    
    B11=b11*matrix(runif(nPrey*nPrey),nPrey,nPrey)
    B12=b12*matrix(runif(nPrey*nHerb),nPrey,nHerb)
    B13=matrix(0,nPrey,nPred)
    B21=b21*matrix(runif(nHerb*nPrey),nHerb,nPrey)
    B22=matrix(0,nHerb,nHerb)
    B23=b23*matrix(runif(nHerb*nPred),nHerb,nPred)
    B31=matrix(0,nPred,nPrey)
    B32=b32*matrix(runif(nPred*nHerb),nPred,nHerb)
    B33=matrix(0,nPred,nPred)
    BB=rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
    diag(BB)<-bdiag1
    diag(BB[(nPrey+nHerb+1):nSpecies,(nPrey+nHerb+1):nSpecies])<-bdiag2
    BB=weight*BB
    B<-BB
  }}}
  
  C1<-c(rep(0.05,nPrey),rep(0,nSpecies-nPrey))
  
  #dispersal####
  disp_mat<-matrix(1/(nCom-1),nCom,nCom)
  diag(disp_mat)<-0
  
  #model####
  X=array(NA,dim=c(nCom,nSpecies,Tmax))
  X[,,1]<-10
  
  for(l in 1:(Tmax-1)){
    X[,,l+1]<-X[,,l]*exp(rep(C1,nCom)+X[,,l]%*%B+A1[Env1,]+A2[Env2,]+Stressor1*Stressor[l]*rep(Stress_tolerance,each=nCom)+Stressor2*Stressor[l]*rep(Stress_tolerance2,each=nCom))+(disp_mat%*%X[,,l])*disp-X[,,l]*disp
    X[,,l+1][(X[,,l+1]<10^-2.5)]<-0
  }
  Abundance_ts<-round(X[,,],digits = 3)*1000
  colnames(Abundance_ts)<-1:nSpecies
  
  Stress.df<-data.frame(Stressor1=Stressor1,Stressor2=Stressor2)
  
  return(list(Abundance_ts=Abundance_ts,traits=traits, Stressor=Stressor,Stress.df=Stress.df))
}