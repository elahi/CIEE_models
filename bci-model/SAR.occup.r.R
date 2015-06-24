# Functions for UBC Working Group on assessing change in biodiversity manuscript
# modified by Sarah Supp, from Fangliang He, based on He and Legendre 20002, and He 2012
# May 4-7, 2015

sar.occup.r = function (abund,p) {
  # plot the effect of habitat loss on # of species
  # this is for the UBC working group on assessing change in biodiversity
  # using He 2012 Ecology paper model
  # May 5, 2015
  
  # p  is the occupancy of a list of species before disturbance, e.g., occup.bci25$occup
  # pc is the occupancy of a list of species after c% abundance being lost after disturbance 
  # but both p and pc have to be converted to proportional data
  
  sar=numeric()
  
  p=p/500000     # 500000 = xmax*ymax from bci coordinates?
  nsp=length(p)  #number of species
  
  #c1=seq(0,1,length=10) #levels for abundance reduction
  c1 = c(0, 0.3, 0.5, 0.8)
  nc1=length(c1)
  
  area=seq(1,50,length=500)
  narea=length(area)
  
  par(mfrow=c(2,2))
  
  # plot SAR for different degrees of habitat loss/abundance reduction
  plot(c(1,50),c(50,300),xlab="Area Sampled",ylab="Number of species",ylim=c(50,300),type="n",log="xy")
  
  for (ii in 1:nc1) {
    for (i in 1:narea) {
      sar.mid=0
      for (j in 1:nsp) { sar.mid=sar.mid+(1-area[i]/50)^(abund[j]*(1-c1[ii])) } 
      sar[i]=nsp-sar.mid
    }
    lines(area[-500],sar[-500],col=ii,lty=ii,lwd=2)
  }
  
  # plot the relationship between abundance reduction and the number of species lost
  plot(c(0,50),c(0,300),xlab="Area Sampled",ylab="Number of species lost",ylim=c(1,300),type="n")
  
  sploss=numeric()
  for (ii in 1:nc1) {
    for (i in 1:narea) {
      sploss.mid=0
      for (j in 1:nsp) { sploss.mid=sploss.mid+(area[i]/50)^(abund[j]-(abund[j]*(c1[ii]))) }
      sploss[i]=sploss.mid
    }
    lines(area,sploss,col=ii,lty=ii,lwd=2)
  }
  
  # plot the relationship between abundance reduction and the number of species retained
  sp.retain=numeric()
  
  # area2=seq(0.05,49,length=20)
  plot(c(0,50),c(0,300),xlab="Area Sampled",ylab="Number of species remaining",ylim=c(1,300),type="n")
  
  for (ii in 1:nc1) {
    for (i in 1:narea) {
      sp.retain.mid=0
      for (j in 1:nsp) { sp.retain.mid=sp.retain.mid+(area[i]/50)^(abund[j]*(1-c1[ii])) }
      sp.retain[i]=nsp-sp.retain.mid
    }
    lines(area,sp.retain,col=ii,lty=ii,lwd=2)
  }
  
  # plot the relationship between effect size and habitat loss
  
  plot(c(0,50),c(-1,1),xlab="Area Sampled",ylab="Effect size",ylim=c(-1,1),type="n")
  
  for (ii in 1:nc1) {
    for (i in 1:narea) {
      sp.retain.mid=0
      for (j in 1:nsp) { sp.retain.mid=sp.retain.mid+(1-area[i]/50)^(abund[j]*(1-c1[ii])) }
      sp.retain[i]=nsp-sp.retain.mid
    }
    if(ii==1){ baseline=sp.retain }
    if(ii>1){ effect.size=log(sp.retain/baseline) 
    lines(area,effect.size,col=ii,lty=ii, lwd=2)
    abline(h=0,lt=2,col=2) 
  }
}

} # END Program
