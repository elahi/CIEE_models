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
  # The program has three reduce.fn. One is "reduce.fn" which randomly removes trees
  # The second program is "reducerare.fn" which aggressively removes rare trees
  # The third program is "reducecommon.fn" which aggressively removes common trees
  # A fourth program "reduce2.fn" which is aggregated removal, is not used in this application
  # 
  # ctfs.dat - BCI stem mapping plot data, e.g., "bci82.dat" which is the 1982 census data after removing all the NAs from "bci.full1".
  # size   - lattice size in meters, that is scale
  # cc - reduction of population size
  # this is main program for reading species x, y coordinates for each species
  
  abund=numeric()
  sp=ctfs.dat$sp
  splist=unique(sp)
  nsp=length(splist)		# no of species
  
  sptable = sort(table(ctfs.dat$sp))
  Nmin = sptable[[1]]                 # number of individuals of the least abundant species
  Nmax = sptable[[length(sptable)]]   # number of individuals of the most abundant species
  
  x=ctfs.dat$gx
  y=ctfs.dat$gy
  
  xmax=1000 # xy coordinates from bci trees
  ymax=500
  
  nxcell=xmax/size		# no of cells along x-axis
  nycell=ymax/size		# no of cells along y-axis
  
  # Build up lists to hold the dataframes for abundance and vector for occupancy
  ntree.rand=list("ntree"=data.frame(abund=rep(-99, nxcell*nycell)), "noccup"=numeric(), "abund"=numeric())
  ntree.rare=list("ntree"=data.frame(abund=rep(-99, nxcell*nycell)), "noccup"=numeric(), "abund"=numeric())
  ntree.common=list("ntree"=data.frame(abund=rep(-99, nxcell*nycell)), "noccup"=numeric(), "abund"=numeric())
  
  # Loop through all the species
  for (i in 1:nsp) {
    
    #identify xy coordinates for a given species
    xx=x[sp==splist[i]]
    yy=y[sp==splist[i]]
    
    xy.dat=data.frame(x=xx,y=yy)

    # Remove trees 1) randomly, 2) biased against rare species, and 3) biased against common species
    xy0.rand=reduce.fn(xy.dat, cc)		# random removal
    xy0.rar=reducerare.fn(xy.dat, cc, Nmin, theta=0.5)    # rare removal
    xy0.com=reducecommon.fn(xy.dat, cc, Nmax, theta=0.5)  # common removal
    
    ntree.rand=reduced.ntree(i, xy0.rand, ntree.rand, abund, xmax, ymax, nxcell, nycell)
    ntree.rare=reduced.ntree(i, xy0.rar, ntree.rare, abund, xmax, ymax, nxcell, nycell)
    ntree.common=reduced.ntree(i, xy0.com, ntree.common, abund, xmax, ymax, nxcell, nycell)
    
    print(i)
  }
  
    
  return(list("rand"=ntree.rand, "rare"=ntree.rare, "common"=ntree.common))
  # return(data.frame(abund=abund,occup=noccup*size*size))
}

