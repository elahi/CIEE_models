# Functions for UBC Working Group on assessing change in biodiversity manuscript
# modified by Sarah Supp, from Fangliang He, based on He and Legendre 20002, and He 2012
# May 4-7, 2015

gridplot.main = function(ctfs.dat,size,cc,plotsize=c(1000,500)) {
  # To calculate and compare the diversity metrics used for the UBC biodiversity change paper
  # May 4-7, 2015
  #
  # This program is based on the program used to calculate IUCN Ecology 2012 paper
  # the reduction of population size is cc = 0.8 (critically endangered), 0.5 (endangered) and 0.3 (vulnerable)
  #
  # The program has two reduce.fn. One is "reduce.fn" which randomly removes trees
  # The second program is "reduce2.fn" which is "aggregated removal"
  # 
  # ctfs.dat - BCI stem mapping plot data, e.g., "bci82.dat" which is the 1982 census data after removing all the NAs from "bci.full1".
  # size   - lattice size in meters, that is scale
  # cc - reduction of population size
  # this is main program for reading species x, y coordinates for each species
  
  
  abund=numeric()    # no of individuals of species
  noccup=numeric()		# no of occupied cells for a species
  
  sp=ctfs.dat$sp
  splist=unique(sp)
  nsp=length(splist)		# no of species
  
  x=ctfs.dat$gx
  y=ctfs.dat$gy
  
  xmax=1000
  ymax=500
  
  nxcell=xmax/size		# no of cells along x-axis
  nycell=ymax/size		# no of cells along y-axis
  
  ntree.dat=data.frame(abund=rep(-99,nxcell*nycell))
  
  #
  for (i in 1:nsp) {
    
    xx=x[sp==splist[i]]
    yy=y[sp==splist[i]]
    
    xy.dat=data.frame(x=xx,y=yy)
    
    xy0.dat=reduce.fn(xy.dat,cc)		# random removal
    #  xy0.dat=reduce2.fn(xy.dat,cc)	# aggregated removal from left to right side of the plot #TODO: Need to bias removal by abundance (e.g. common vs. rare)
    
    x0=xy0.dat$x
    y0=xy0.dat$y
    
    abund[i]=length(x0)
    
    plotxy.fn(x0,y0,xmax,ymax)
    
    # call program presence.fn which converts the points into presence/absence data
    zz=presence.fn(size,nxcell,nycell,x0,y0,abund[i],xmax,ymax)
    
    
    print(i)
    noccup[i]=zz$noccup		# no of occupied cells
    ntree.dat=data.frame(ntree.dat,ntree=zz$npt)
  }
  
  return(ntree.dat)
  # return(data.frame(abund=abund,occup=noccup*size*size))
}

