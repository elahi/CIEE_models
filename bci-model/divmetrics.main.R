# Functions for UBC Working Group on assessing change in biodiversity manuscript
# modified by Sarah Supp, from Fangliang He, based on He and Legendre 20002, and He 2012
# May 4-7, 2015

divmetrics.main = function(data0,data03,data05,data08) {
  # calculate diversity change using: 
  # log SR, log Shannon, log Simpson, slope SR, slope abundance, slope Shannon, slope Simpson, 
  # Jaccard, Bray-Curtis
  
  require(reshape)
  require(reshape2)
  
  rich00=numeric()
  rich03=numeric()
  rich05=numeric()
  rich08=numeric()
  
  Hshannon00=numeric()
  Hshannon03=numeric() 
  Hshannon05=numeric()
  Hshannon08=numeric()
  
  Hsimpson00=numeric()
  Hsimpson03=numeric()
  Hsimpson05=numeric()
  Hsimpson08=numeric()
  
  bc03=numeric()
  bc05=numeric()
  bc08=numeric()
  
  data0=data0[,-1]
  data03=data03[,-1]
  data05=data05[,-1]
  data08=data08[,-1]
  
  ncell=dim(data0)[1]
  
  print(ncell)
  
  for (i in 1:ncell) {
    sp00=data0[i,]
    sp03=data03[i,]
    sp05=data05[i,]
    sp08=data08[i,]
    
    # richness effect size
    rich00[i]=length(sp00[sp00>0])
    rich03[i]=length(sp03[sp03>0])
    rich05[i]=length(sp05[sp05>0])
    rich08[i]=length(sp08[sp08>0])
    
    # Hill's Shannon
    p00=sp00/sum(sp00)
    p03=sp03/sum(sp03)
    p05=sp05/sum(sp05)
    p08=sp08/sum(sp08)
    
    Hshannon00[i]=exp(-sum(p00[p00>0]*log(p00[p00>0])))
    Hshannon03[i]=exp(-sum(p03[p03>0]*log(p03[p03>0]))) 
    Hshannon05[i]=exp(-sum(p05[p05>0]*log(p05[p05>0])))
    Hshannon08[i]=exp(-sum(p08[p08>0]*log(p08[p08>0])))
    
    # Hill's Simpson
    Hsimpson00[i]=1/sum(p00^2)
    Hsimpson03[i]=1/sum(p03^2) 
    Hsimpson05[i]=1/sum(p05^2)
    Hsimpson08[i]=1/sum(p08^2)
    
    # Bray-Curtis
    bc03[i]=sum(abs(sp00-sp03))/sum(sp00+sp03)
    bc05[i]=sum(abs(sp00-sp05))/sum(sp00+sp05)
    bc08[i]=sum(abs(sp00-sp08))/sum(sp00+sp08)
    
    print(i)
  }
  
  #dataframe with everything
  metric = data.frame(stress=c(rep(0,ncell),rep(0.3,ncell),rep(0.5,ncell),rep(0.8,ncell)),
                               Richness=c(rich00, rich03, rich05, rich08),
                               Hshannon=c(Hshannon00,Hshannon03,Hshannon05,Hshannon08),
                               Hsimpson=c(Hsimpson00,Hsimpson03,Hsimpson05,Hsimpson08))
  
  metric.melt = melt(metric, id.vars="stress")
  names(metric.melt) = c("stress", "metric", "value")
  
  fx = data.frame(stress=c(rep(0.3,ncell),rep(0.5,ncell),rep(0.8,ncell)),
                           rich.abs=c(rich03-rich00,rich05-rich00,rich08-rich00),
                           Hshannon.abs=c(Hshannon03-Hshannon00,Hshannon05-Hshannon00,Hshannon08-Hshannon00),
                           Hsimpson.abs=c(Hsimpson03-Hsimpson00,Hsimpson05-Hsimpson00,Hsimpson08-Hsimpson00),
                           bc=c(bc03,bc05,bc08),
                           rich.es=c(log(rich03/rich00),log(rich05/rich00),log(rich08/rich00)),
                           Hshannon.es=c(log(Hshannon03/Hshannon00),log(Hshannon05/Hshannon00),log(Hshannon08/Hshannon00)),
                           Hsimpson.es=c(log(Hsimpson03/Hsimpson00),log(Hsimpson05/Hsimpson00),log(Hsimpson08/Hsimpson00)))
  
  fx.melt = melt(fx, id.vars="stress")
  names(fx.melt) = c("stress", "metric", "value")
  
  vals = rbind(fx.melt, metric.melt)
  
  return(vals)
}
  


# #the metrics
# Richness.out=data.frame(stress2, Richness=c(rich00, rich03, rich05, rich08))
# Hshannon.out=data.frame(stress2,Hshannon=c(Hshannon00,Hshannon03,Hshannon05,Hshannon08))
# Hsimpson.out=data.frame(stress2,Hsimpson=c(Hsimpson00,Hsimpson03,Hsimpson05,Hsimpson08))
# bc.out=data.frame(stress1,bc=c(bc03,bc05,bc08))
# 
# #effect size of metrics
# rich.es.out=data.frame(stress1, rich.es=c(log(rich03/rich00),log(rich05/rich00),log(rich08/rich00)))
# Hshannon.es.out=data.frame(stress1, Hshannon.es=c(log(Hshannon03/Hshannon00),log(Hshannon05/Hshannon00),log(Hshannon08/Hshannon00)))
# Hsimpson.es.out=data.frame(stress1, Hsimpson.es=c(log(Hsimpson03/Hsimpson00),log(Hsimpson05/Hsimpson00),log(Hsimpson08/Hsimpson00)))
# 
# #absolute difference in metrics
# rich.abs.out=data.frame(stress1, rich.abs=c(rich00-rich03,rich00-rich05,rich00-rich08))
# Hshannon.abs.out=data.frame(stress1, Hshannon.abs=c(Hshannon00-Hshannon03,Hshannon00-Hshannon05,Hshannon00-Hshannon08))
# Hsimpson.abs.out=data.frame(stress1, Hsimpson.abs=c(Hsimpson00-Hsimpson03,Hsimpson00-Hsimpson05,Hsimpson00-Hsimpson08))

#   par(mfrow=c(1,4))
#   
#   boxplot(Richness~stress2, data=Richness.out, main="Richness")
#   boxplot(Hshannon~stress2, data=Hshannon.out,main="Hill Shannon")
#   boxplot(Hsimpson~stress2, data=Hsimpson.out, main="Hill Simpson")
#   boxplot(bc~stress1, data=bc.out, main="Bray-Curtis")
#   
#   par(mfrow=c(2,3))
#   boxplot(rich.abs~stress1, data=rich.abs.out,main="Richness difference")
#   boxplot(Hshannon.abs~stress1, data=Hshannon.abs.out,main="Hill Shannon difference")
#   boxplot(Hsimpson.abs~stress1, data=Hsimpson.abs.out,main="Hill Simpson difference")
#   
#   boxplot(rich.es~stress1, data=rich.es.out,main="Richness effect size")
#   boxplot(Hshannon.es~stress1, data=Hshannon.es.out,main="Hill Shannon effect size")
#   boxplot(Hsimpson.es~stress1, data=Hsimpson.es.out,main="Hill Simpson effect size")
  
