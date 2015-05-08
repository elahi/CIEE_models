# This script is for exploring Fangliang's analytical model, based on BCI tree data

# Load packages, source functions, and data

library(ggplot2)

source('BCI-model-fxns.R')
source('gridplot.main.R')
source('SAR.occup.R')
source('divmetrics.main.R')

load('bci82.dat')


# Model BCI plots with different levels of disturbance, where disturbance is random removal of individual trees
# set parameters
size = c(10, 25, 50, 100)
cc = c(0, 0.3, 0.5, 0.8)

for (s in size){
  #TODO:
  # this will need to be fine-tuned to calculate abundance reduction in 3 different ways.
  # I think the current reduce2 works on spatial aggregation, 
  # when what we want is reduction on common vs. rare species biased proportionally
  ntree.dat1 = gridplot.main(bci82.dat, size[s], cc[1], plotsize=c(1000,500))
    ntree.dat2 = gridplot.main(bci82.dat, size[s], cc[2], plotsize=c(1000,500))
    ntree.dat3 = gridplot.main(bci82.dat, size[s], cc[3], plotsize=c(1000,500))
    ntree.dat4 = gridplot.main(bci82.dat, size[s], cc[4], plotsize=c(1000,500))
    
    scale = divmetrics.main(ntree.dat1, ntree.dat2, ntree.dat3, ntree.dat4)
  
  if(s = 1){
    scale$scale = size[s]
    results = scale
  }
  if (s>1){
    scale$scale = size[s]
    results = rbind(results, scale)
  }
}

#ggplot needs work!
absdiff = subset(scale, metric %in% c("rich.abs", "Hshannon.abs", "Hsimpson.abs"))
ggplot(absdiff, aes(scale,value, group=stress)) + geom_boxplot() + facet_wrap(~metric)







