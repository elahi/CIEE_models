# This script is for exploring Fangliang He's analytical model, based on BCI tree data
# Code written by Fangliang He, based on He and Legendre 2002 and He 2012
# Code modified by Sarah Supp
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
size = c(25,100) #main plots to focus on 25 and 100, could also include 10 and 50
cc = c(0, 0.3, 0.5, 0.8) #main plots to focus on 0.5 and 0.8

for (s in 1:length(size)){
  # Calculate abundance reduction in 3 different ways:
  # 1) random reduction
  # 2) rare-biased reduction
  # 3) common-biased reduction
   # TODO: make this section flexible so it can handle different lenghts of cc (to do so, divmetrics.main would also need to be modified)
    #returns a list of lists holding an ntree dataframe for each scenario and an occupancy vector
    ntree.dat1 = gridplot.main(bci82.dat, size[s], cc[1], plotsize=c(1000,500))
    ntree.dat2 = gridplot.main(bci82.dat, size[s], cc[2], plotsize=c(1000,500))
    ntree.dat3 = gridplot.main(bci82.dat, size[s], cc[3], plotsize=c(1000,500))
    ntree.dat4 = gridplot.main(bci82.dat, size[s], cc[4], plotsize=c(1000,500))
    
    #plot effect of abundance reduction on species
    #FIXME: Ideally, this should be an apply statement that inputs from the scenarios, and labels each panel
    SAR.occup.r(ntree.dat1$rand$abund, ntree.dat1$rand$noccup) #TODO
    
    # calculate diversity metrics for the three scenarios
    scale.rand = divmetrics.main(ntree.dat1[[1]][1], ntree.dat2[[1]][1], ntree.dat3[[1]][1], ntree.dat4[[1]][1], size[s])
    scale.rare = divmetrics.main(ntree.dat1[[2]][1], ntree.dat2[[2]][1], ntree.dat3[[2]][1], ntree.dat4[[2]][1], size[s])
    scale.comm = divmetrics.main(ntree.dat1[[3]][1], ntree.dat2[[3]][1], ntree.dat3[[3]][1], ntree.dat4[[3]][1], size[s])
  
    # add tags for removal scenario
    scale.rand[[1]]$removal = rep("random", nrow(scale.rand[[1]]))
    scale.rand[[2]]$removal = rep("random", nrow(scale.rand[[2]]))
    
    scale.rare[[1]]$removal = rep("rare", nrow(scale.rare[[1]]))
    scale.rare[[2]]$removal = rep("rare", nrow(scale.rare[[2]]))
    
    scale.comm[[1]]$removal = rep("common", nrow(scale.comm[[1]]))
    scale.comm[[2]]$removal = rep("common", nrow(scale.comm[[2]]))
    
  if(s == 1){
    raw = rbind(scale.rand[[1]], scale.rare[[1]], scale.comm[[1]])
    effect = rbind(scale.rand[[2]], scale.rare[[2]], scale.comm[[2]])
  }
  if (s > 1){
    raw = rbind(raw, scale.rand[[1]], scale.rare[[1]], scale.comm[[1]])
    effect = rbind(effect, scale.rand[[2]], scale.rare[[2]], scale.comm[[2]])
  }
}

# convert values to factors
raw$scale = as.factor(raw$scale)
raw$stress = as.factor(raw$stress)
raw$removal = factor(raw$removal, levels=c("random", "common", "rare"))
effect$scale = as.factor(effect$scale)
effect$stress = as.factor(effect$stress)
effect$removal = factor(effect$removal, levels=c("random", "common", "rare"))

# make a melted version of the dataframe for plotting
raw.melt = melt(raw, id.vars=c("stress", "scale", "removal"))
  names(raw.melt) = c("stress", "scale", "removal", "metric", "value")
effect.melt = melt(effect, id.vars=c("stress", "scale", "removal"))
  names(effect.melt) = c("stress", "scale", "removal", "metric", "value")

# Plot results for Figure 5 in UBC manuscript
#  the dataframe scale should hold all the necessary data, but will need to be subsetted to plot the values you want for the results
#  We want several panels for Figure 5 and for supplementary figures in the appendix:
#            There are three plots in a row (for Richness, Shannon Diversity, and Simpson's)
#                  and two rows (for absolute change in metric and LRR effect size)
#            x-axis is scale (we discussed just using size = c(25, 100) for the main fig)
#            y-axis is either abolute difference or effect size (use 0.5 and/or 0.8 for the main fig)
#            data in the figure is grouped by the three reduction scenarios (random, common-biased, & rare-biased removal)

absdiff = subset(effect.melt, metric %in% c("rich.abs","Hshannon.abs", "Hsimpson.abs") & stress %in% c(0.5))
p1 = ggplot(absdiff, aes(scale,value, group=interaction(scale, removal, stress))) + geom_boxplot(aes(fill=removal)) + facet_wrap(~metric, scales="free") + theme_bw() +
  ylab("absolute difference")
effectsize = subset(effect.melt, metric %in% c("rich.es", "Hshannon.es", "Hsimpson.es") & stress %in% c(0.5))
p2 = ggplot(effectsize, aes(scale,value, group=interaction(scale, removal, stress))) + geom_boxplot(aes(fill=removal)) + facet_wrap(~metric) + theme_bw() + 
  ylab("LRR Effect Size")
grid.arrange(p1, p2)


# Plot the raw measures for Richness, Shannon, Simpson, and Bray-Curtis
bcmetric = subset(effect.melt, metric %in% c("bc") & stress %in% c(0, 0.5))
raw.melt2 = rbind(raw.melt, bcmetric)
raw.melt2 = subset(raw.melt, stress %in% c(0, 0.5))
ggplot(raw.melt2, aes(scale,value, group=interaction(scale, stress, removal))) + geom_boxplot(aes(fill=removal)) + facet_wrap(~metric, scales="free") + theme_bw()

#
# # Code to make the triangle plots. 
#The data frame for making the Triangle plot has the three measures of diversity, route (stress), and time (scale). 
raw2 = subset(raw, stress==c(0,0.5) & removal=="random" & scale==25)
ggtern(raw2, aes(x=Richness, y=Hshannon, z=Hsimpson, group=scale, color=stress)) +
#  geom_point(aes(shape=stress), size=4, alpha=0.5) +
  geom_point(alpha=0.5) +
  tern_limits(T=.7,L=0.95,R=0.7) +
  theme_bw(base_size = 16) + 
  scale_color_brewer(type = "qual",palette = 6)

#The data frame for making the Triangle plot has the three measures of diversity, route (stress), and time (scale). 
# effectLRR = effect[effect$stress==0.8,c(1,2,10,7,8,9)]
# ggtern(effectLRR, aes(x=rich.es, y=Hshannon.es, z=Hsimpson.es, group=scale, color=removal)) +
#   geom_point() +
#   tern_limits(T=.8,L=0.9,R=0.8) +
#   theme_bw(base_size = 16) + 
#   scale_color_brewer(type = "qual",palette = 6)


