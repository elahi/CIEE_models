# This script is for exploring Fangliang He's analytical model, based on BCI tree data
# Code written by Fangliang He, based on He and Legendre 2002 and He 2012
# Code modified by Sarah Supp and Robin Elahi
# Model simulation for UBC working group on Biodiversity Change, 4-7 May 2015

# The folder contains several files for replicating the results of UBC manuscript
#   1. model-BCI.R       : Runs the code and outputs the main results and figures 
#                         (e.g. Figure 5 for UBC manuscript)
#   2. gridplot.main.R   : Reduces BCI trees based on a random or biased reduction algorithm, 
#                         using IUCN population reduction levels (0, 0.3, 0.5, 0.8)
#                          outputs a site by species matrix for abundance
#   3. divmetrics.main.R : Calculates diversity metrics and effect sizes from gridplot.main.R output (S, H, D, Bray-Curtis) 
#   4. BCI-model-fxns.R  : Contains additional functions for running the main scripts for analysis
#   5. SAR.occup.r.R     : Plots the effect of abundance reduction on the number of species sampled. Outputs species area curves (e.g. Figure for Box 1 for UBC manuscript).
#   We can discuss if there is a better way to organize these scripts later (e.g. consolidate into fewer separate .R files)

#Load FH RData, which contains the dataset to model on (bci82.dat) and all the functions (which I've since put into the R scripts so we can see and modify them as needed)
# Please do not edit the RData file, it contains FH code and data, but we should work on a script to replicate the analyses
# TODO: Save just bci82.dat out as a small datafile to load here, instead of the whole RData file (which isn't needed)
load('/Users/sarah/Documents/GitHub/CIEE_models/bci-model/DiversityMetrics.RData')

# Load packages and source functions
library(ggplot2)
library(ggtern)
library(reshape)
library(reshape2)
library(gridExtra)

source('BCI-model-fxns.R')
source('gridplot.main.R')
source('SAR.occup.r.R')
source('divmetrics.main.R')


# Model BCI plots with different levels of disturbance, where disturbance is random removal of individual trees
# set parameters
size = c(10, 25, 50, 100) #main plots to focus on 25 and 100
cc = c(0, 0.3, 0.5, 0.8) #main plots to focus on 0.5 and 0.8

for (s in 1:length(size)){
  #TODO:
  # this will need to be fine-tuned to calculate abundance reduction in 3 different ways.
  # right now gridplot.main just does random reduction (reduce.fn)
  # we should reparameterize it so that it can do the reduction in 3 different ways (random, common, rare)
  # I think the current reduce2.fn works on spatial aggregation, but what we want is reduction on common vs. rare species biased proportionally
   # TODO: make this section flexible so it can handle different lenghts of cc (to do so, divmetrics.main would also need to be modified)
    ntree.dat1 = gridplot.main(bci82.dat, size[s], cc[1], plotsize=c(1000,500))
    ntree.dat2 = gridplot.main(bci82.dat, size[s], cc[2], plotsize=c(1000,500))
    ntree.dat3 = gridplot.main(bci82.dat, size[s], cc[3], plotsize=c(1000,500))
    ntree.dat4 = gridplot.main(bci82.dat, size[s], cc[4], plotsize=c(1000,500))
    
    scale = divmetrics.main(ntree.dat1, ntree.dat2, ntree.dat3, ntree.dat4, size[s])
  
  if(s == 1){
    raw = scale[[1]]
    effect = scale[[2]]
  }
  if (s > 1){
    raw = rbind(raw, scale[[1]])
    effect = rbind(effect, scale[[2]])
  }
}


raw$scale = as.factor(raw$scale)
raw$stress = as.factor(raw$stress)
effect$scale = as.factor(effect$scale)
effect$stress = as.factor(effect$stress)

# make a melted version of the dataframe for plotting
raw.melt = melt(raw, id.vars=c("stress", "scale"))
  names(raw.melt) = c("stress", "scale", "metric", "value")
effect.melt = melt(effect, id.vars=c("stress", "scale"))
  names(effect.melt) = c("stress", "scale", "metric", "value")

# Plot results for Figure 5 in UBC manuscript
#  the dataframe scale should hold all the necessary data, but will need to be subsetted to plot the values you want for the results
#  We want several panels for Figure 5 and for supplementary figures in the appendix:
#            There are three plots in a row (for Richness, Shannon Diversity, and Simpson's)
#                  and two rows (for absolute change in metric and LRR effect size)
#            x-axis is scale (we discussed just using size = c(25, 100) for the main fig)
#            y-axis is either abolute difference or effect size (use 0.5 and/or 0.8 for the main fig)
#            data in the figure is grouped by the three reduction scenarios (random, common-biased, & rare-biased removal)

absdiff = subset(effect.melt, metric %in% c("rich.abs","Hshannon.abs", "Hsimpson.abs"))
p1 = ggplot(absdiff, aes(scale,value, group=interaction(scale, stress))) + geom_boxplot(aes(fill=stress)) + facet_wrap(~metric, scales="free") + theme_bw() +
  ylab("absolute difference")
effectsize = subset(effect.melt, metric %in% c("rich.es", "Hshannon.es", "Hsimpson.es"))
p2 = ggplot(effectsize, aes(scale,value, group=interaction(scale,stress))) + geom_boxplot(aes(fill=stress)) + facet_wrap(~metric) + theme_bw() + 
  ylab("LRR Effect Size")

grid.arrange(p1, p2)

# Plot the raw measures for Richness, Shannon, Simpson, and Bray-Curtis
bcmetric = subset(effect.melt, metric %in% c("bc"))
raw.melt2 = rbind(raw.melt, bcmetric)
ggplot(raw.melt2, aes(scale,value, group=interaction(scale,stress))) + geom_boxplot(aes(fill=stress)) + facet_wrap(~metric, scales="free") + theme_bw()

#
# # Code to make the triangle plots. 
#The data frame for making the Triangle plot has the three measures of diversity, route (stress), and time (scale). 
ggtern(raw, aes(x=Richness, y=Hshannon, z=Hsimpson, group=stress, color=scale)) +
  geom_point() +
  tern_limits(T=.6,L=0.8,R=0.5) +
  theme_bw(base_size = 16) + 
  scale_color_brewer(type = "qual",palette = 6)


