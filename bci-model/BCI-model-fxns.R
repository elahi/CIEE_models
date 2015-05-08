# Functions for UBC Working Group on assessing change in biodiversity manuscript
# modified from Fangliang He, based on He and Legendre 20002, and He 2012
# May 4-7, 2015


#------------------------------- Functions for sampling and manipulating BCI dataset
count.fn = function(size,nxcell,nycell,x,y){
  # 1. count the no of points (trees) in a cell
  
  # Divide the plot into a grid system with cell size = size
  npt=numeric()           # no of points in each cell
  z=numeric()             # presence in each cell
  
  
  ncell=0                 # total no of cells
  
  xlo=-size
  
  for (i in 1:nxcell){
    xlo=xlo+size
    xup=xlo+size
    
    ylo=-size
    
    for (j in 1:nycell){
      ylo=ylo+size
      yup=ylo+size
      
      ncell=ncell+1
      
      npt[ncell]=length(x[(x>xlo&x<=xup)&(y>ylo&y<=yup)])
    }
  }
  
  return(npt)
}


matrix.fn = function(z) {
  # This function does WHAT FILL IN THE BLANK
  norow=dim(z)[1]+1
  nocol=dim(z)[2]+1
  
  x=matrix(0,nrow=norow,ncol=nocol)
  
  for(i in 2:norow){
    for(j in 2:nocol){
      
      x[i,j]=z[i-1,j-1]
    }
  }
  return(c(t(x))) 
}


plot.figure.r = function (data1, data2, cc) {
  #Plots FILL IN THE BLANK
  p = data1$occup/5e+05
  pc = data2$occup/5e+05
  plot(p, pc, xlab = "", ylab = "", cex = 0.7, ylim = c(0, 
                                                        1), xlim = c(0, 1))
  curve(1 - (1 - x)^(1 - cc), add = T, col = "red", lwd = 1)
  curve((1 - cc) * x, add = T, col = "blue", lwd = 1)
}


plot.figure3.r = function (data1, data2, data3, cc) {
  #Plots FILL IN THE BLANK
  p = data1$occup/5e+05
  pc.rand = data2$occup/5e+05
  pc.syst = data3$occup/5e+05
  pc = (pc.rand + pc.syst)/2
  plot(p, pc, xlab = "", ylab = "", cex = 0.7, ylim = c(0, 
                                                        1), xlim = c(0, 1))
  curve(1 - (1 - x)^(1 - cc), add = T, col = "red", lwd = 1)
  curve((1 - cc) * x, add = T, col = "blue", lwd = 1)
}


plot.r = function (data1, data2, data3, cc) {
  # Plots FILL IN THE BLANK
  pobs = data1$occup/5e+05
  pc = data2$occup/5e+05
  pc2 = data3$occup/5e+05
  pc.mean = (pc + pc2)/2
  prd = 1 - (1 - pobs)^(1 - cc)
  print(cor(pc.mean, prd))
  plot(pobs, pc.mean, xlim = c(0, 1), ylim = c(0, 1), xlab = "", 
       ylab = "", cex = 0.8)
  curve(1 - (1 - x)^(1 - cc), col = "red", add = T)
  curve((1 - cc) * x, col = "blue", add = T)
}


plotxy.fn = function(x,y,xmax,ymax){
  #Plot FILL IN THE BLANK
  par(mfrow=c(2,1),mai=c(0.45,1.1,0.2,1.1))
  plot(x,y,pch=16,cex=0.25,xlim=c(0,xmax),ylim=c(0,ymax),xlab="",ylab="",col="blue")
  lines(c(0,0),c(0,ymax),col="red")
  lines(c(0,xmax),c(ymax,ymax),col="red")
  lines(c(0,xmax),c(0,0),col="red")
  lines(c(xmax,xmax),c(0,ymax),col="red")
}


presence.fn = function(size,nxcell,nycell,x,y,abund,xmax,ymax){
  # convert cell count into presencen/absence data.
  
  xx=seq(0,xmax,len=nxcell+1)
  yy=seq(0,ymax,len=nycell+1)
  
  npt=count.fn(size,nxcell,nycell,x,y)  	# call function count.fn for observed pattern
  
  z=ifelse(npt>0,1,0)				# convert abundance npt into presence/absence
  
  noccup=sum(z)
  
  #  z=matrix(z,ncol=nxcell)
  npt=matrix(npt,ncol=nxcell)
  
  plot(x,y,xlim=c(0,xmax),ylim=c(0,ymax),xlab="",ylab="",type="n")
  #  image(xx,yy,t(z),breaks=c(0,0.5,1),col=c("white","blue"),xlab="",ylab="",add=T)	# map presence/absence map
  image(xx,yy,t(npt),col=topo.colors(12),xlab="",ylab="",add=T)				# map the number of trees (tree count)
  
  for (i in 1:(nxcell+1)){
    lines(c(xx[i],xx[i]),c(0,ymax),col=2)
  }
  for (i in 1:(nycell+1)){
    lines(c(0,xmax),c(yy[i],yy[i]),col=2) 
  }
  return(list(npt=as.vector(npt),noccup=noccup))
}


reduce.fn = function(xy.dat,cc){
  # randomly removal of c proportion of trees
  # c=0.8 (critically endangered)
  # c=0.5 (endangered)
  # c=0.2 (vulnerable)
  
  abund=length(xy.dat$x)
  
  if(abund<5) {
    x=runif(1,0,1000)
    y=runif(1,0,500)
  }
  
  else{
    ntree=round((1-cc)*abund)  	# reduce cc number of trees. This is opposite to "reduce2.fn"
    iseq=sample(1:abund,ntree)		# randomly sample ntree from sequence 1:abund
    
    x=xy.dat$x[iseq]
    y=xy.dat$y[iseq]
  }
  return(data.frame(x=x,y=y))
}


reduce2.fn = function(xy.dat,cc){
  # aggregated removal of c proportion of trees
  # c=0.8 (critically endangered)
  # c=0.5 (endangered)
  # c=0.2 (vulnerable)
  
  abund=length(xy.dat$x)
  
  if(abund<5) {
    x=runif(1,0,1000)
    y=runif(1,0,500)
  }
  
  else{
    x=xy.dat$x
    
    norder=order(x)
    
    xy.dat=xy.dat[norder,]  		#sort xy.dat according to x-axis
    
    ntree=round(cc*abund)			# reduce cc number of trees
    iseq=(ntree+1):abund			# retain those trees with x-axis larger than the n-th tree.
    
    xy.dat=xy.dat[iseq,]			# retain both x and y axes
    x=xy.dat$x
    y=xy.dat$y
  }
  return(data.frame(x=x,y=y))
}
