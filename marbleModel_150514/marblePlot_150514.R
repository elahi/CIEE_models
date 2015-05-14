# RE 150514
# Figure 1
# plots

# remove everything from R's memory
rm(list=ls(all=TRUE)) 

# load packages
library(ggplot2)
library(plyr)
library(vegan)
source("./summarizeData_150204.R")
source("./multiplotF.R")

# set disturbance size (# of individuals removed)
distSize <- 200

# run 10 spp simulation
source("./marble10spp_150514.R")
spp10dat

# run 20 spp simulation
source("./marble20spp_150514.R")
spp20dat

# combine datasets
dat <- rbind(spp10dat, spp20dat)

summary(dat)
dat$scale <- factor(dat$scale)
dat$pool <- factor(dat$pool)
dat$time <- factor(dat$time)

names(dat)[1:3] <- c("hill1", "hill2", "hill3")

head(dat)

summaryDat <- summarySE(data = dat, measurevar = "hill1", 
	groupvars = c("scale", "pool", "time"))
summaryDat

# rename factor levels for 'scale' and 'time'
summaryDat$Scale <- factor(c(rep("Regional", 4), rep("Local", 4)))
timeList <- unique(summaryDat$time)
TimeList <- c("Post-disturbance", "Pre-disturbance")
summaryDat$Census <- mapvalues(summaryDat$time, from = timeList, to = TimeList)
# relevel Scale
summaryDat$Scale <- factor(summaryDat$Scale, levels = c("Regional", "Local"))
# relevel Scale
summaryDat$Census <- factor(summaryDat$Census, levels = c("Pre-disturbance", "Post-disturbance"))

### plot this
ULClabel <- theme(plot.title = element_text(hjust = -0.025, vjust = 0.01, size = rel(1)))

summaryPlot <- ggplot(summaryDat, 
	aes(x = pool, y = hill1, color = Scale, shape = Census)) + 
	theme_classic(base_size = 12) + xlab("Regional pool") +
	ylab("Species richness") + 
	geom_point(size = 3) + 
	geom_errorbar(aes(ymin = hill1 - ci, ymax = hill1 + ci), 
		width = 0.02) +
	theme(legend.justification = c(1,0), legend.position = c(0.6, 0.5)) +
	labs(title = "A") + ULClabel 
	
summaryPlot	
##############################################
##############################################
# Calculating change
##############################################
##############################################
summary(dat)
scaleList <- unique(dat$scale)
scaleList
ScaleList <- c("Regional", "Local")
dat$Scale <- mapvalues(dat$scale, from = scaleList, to = ScaleList)
unique(dat$Scale)

scalePoolTimeF <- function(x) return(as.factor(paste(x$scale, 
	substr(x$pool, 1, 2), x$time, sep = "_")))

dat$cat <- scalePoolTimeF(dat)
head(dat)

datPre <- dat[dat$time == "pre", ]
datPost <- dat[dat$time == "post", ]

absChangeF <- function(x, y) {
	crap <- ifelse(x$cat == "global_10_post", 
		x$hill1 - y[y$cat == "global_10_pre", ]$hill1, 
				ifelse(x$cat == "global_20_post", 
		x$hill1 - y[y$cat == "global_20_pre", ]$hill1, 
				ifelse(x$cat == "local_10_post", 
		x$hill1 - y[y$cat == "local_10_pre", ]$hill1, 
		
		x$hill1 - y[y$cat == "local_20_pre", ]$hill1)))
		
	return(crap)
}

logChangeF <- function(x, y) {
	crap <- ifelse(x$cat == "global_10_post", 
		log(x$hill1/y[y$cat == "global_10_pre", ]$hill1), 
				ifelse(x$cat == "global_20_post", 
		log(x$hill1/y[y$cat == "global_20_pre", ]$hill1), 
				ifelse(x$cat == "local_10_post", 
		log(x$hill1/y[y$cat == "local_10_pre", ]$hill1), 
		
		log(x$hill1/y[y$cat == "local_20_pre", ]$hill1) 
		)))
	return(crap)
}

logChangeF(datPost, datPre)

datPost$absChange <- absChangeF(datPost, datPre)
head(datPost)

datPost$logChange <- logChangeF(datPost, datPre)
head(datPost)
summary(datPost)

absChangeDF <- summarySE(data = datPost, 
	measurevar = "absChange", 
	groupvars = c("Scale", "pool"))

logChangeDF <- summarySE(data = datPost, 
	measurevar = "logChange", 
	groupvars = c("Scale", "pool"))

ULClabel <- theme(plot.title = element_text(hjust = -0.05, vjust = 0.01, size = rel(1)))

absChangePlot <- ggplot(absChangeDF, 
	aes(x = pool, y = absChange, color = Scale)) + 
	theme_classic(base_size = 12) + xlab("Regional pool") +
	ylab("Absolute change") + 
	geom_point(size = 3) + 
	geom_errorbar(aes(ymin = absChange - ci, ymax = absChange + ci), 
		width = 0.02) +
	theme(legend.justification = c(1,0), legend.position = c(0.9, 0.75)) +
	labs(title = "B") + ULClabel +
	geom_hline(yintercept = 0, linetype = 'dashed') +
	theme(legend.position = "none")

logChangePlot <- ggplot(logChangeDF, 
	aes(x = pool, y = logChange, color = Scale)) + 
	theme_classic(base_size = 12) + xlab("Regional pool") +
	ylab("Log change") + 
	geom_point(size = 3) + 
	geom_errorbar(aes(ymin = logChange - ci, ymax = logChange + ci), 
		width = 0.02) +
	theme(legend.justification = c(1,0), legend.position = c(0.9, 0.75)) +
	labs(title = "C") + ULClabel +
	geom_hline(yintercept = 0, linetype = 'dashed') +
	theme(legend.position = "none")

###
multiplot(panelA, panelB, summaryPlot, 
	absChangePlot, logChangePlot, 
	layout = matrix(c(1, 2, 3, 3, 4, 5), nrow = 2, byrow = FALSE))
