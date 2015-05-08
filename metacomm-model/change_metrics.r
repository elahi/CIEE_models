Temporal_Jaccard<-function(pre,post){
  if(class(pre)!="matrix"){beta<-vegdist(rbind(pre,post),method = "jaccard")} else 
  {beta<-rep(NA,nrow(pre))
   for(i in 1:nrow(pre)){
     beta[i]<-vegdist(rbind(pre[i,],post[i,]),method = "jaccard")}}
  return(beta)
}

Temporal_Bray<-function(pre,post){
  if(class(pre)!="matrix"){beta<-vegdist(rbind(pre,post),method = "bray")} else 
  {beta<-rep(NA,nrow(pre))
   for(i in 1:nrow(pre)){
     beta[i]<-vegdist(rbind(pre[i,],post[i,]),method = "bray")}}
  return(beta)
}


change_metrics<-function(x,means=T){
  SR_ts<-apply(x>0,3,rowSums)
  Com_abund_ts<-apply(x,3,rowSums)
  Shan_ts<-apply(x,3,renyi,hill=T,scales=c(1))
  Simp_ts<-apply(x,3,renyi,hill=T,scales=c(2))
  
  log_SR<-log(SR_ts[,-1]/SR_ts[,1])
  log_abund<-log(Com_abund_ts[,-1]/Com_abund_ts[,1])
  log_shan<-log(Shan_ts[,-1]/Shan_ts[,1])
  log_simp<-log(Simp_ts[,-1]/Simp_ts[,1])
  
  Jaccard<-log_simp
  for(i in 2:dim(x)[3]){
    Jaccard[,i-1]<-Temporal_Jaccard(pre = x[,,1],post = x[,,i])
  }
  
  Bray_Curtis<-log_simp
  for(i in 2:dim(x)[3]){
    Bray_Curtis[,i-1]<-Temporal_Bray(pre = x[,,1],post = x[,,i])
  }
  
  if(means==T){Change.df<-data.frame(log_SR=colMeans(log_SR),log_Abund=colMeans(log_abund),log_Shannon=colMeans(log_shan),log_Simpson=colMeans(log_simp),Jaccard=colMeans(Jaccard),Bray_Curtis=colMeans(Bray_Curtis),Stress=Stressor_sub[-1])
  } else{Change.df<-data.frame(log_SR=c(log_SR),log_Abund=c(log_abund),log_Shannon=c(log_shan),log_Simpson=c(log_simp),Jaccard=c(Jaccard),Bray_Curtis=c(Bray_Curtis),Stress=rep(Stressor_sub[-1],each=dim(x)[1]))}
  return(Change.df)
}

