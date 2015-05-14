# RE 150514
# Figure 1 simulation
# 20 species in regional pool

###############
set.seed(126)

#20 total species
# species pool
commSpp <- c(1:4) # 
medSpp <- c(5:8) # 
medLoSpp <- c(9:12)
rareSpp <- c(13:20) # 

# random sampling of species, 400 total individuals
commonSample <- sample(commSpp, replace = TRUE, size = 150)
medSample <- sample(medSpp, replace = TRUE, size = 125)
medLoSample <- sample(medLoSpp, replace = TRUE, size = 75)
rareSample <- sample(rareSpp, replace = TRUE, size = 50)

# combine into a vector
commSample <- c(commonSample, medSample, medLoSample, rareSample)
commIdentity <- c(rep("common", 175), rep("med", 125), 
	rep("medLo", 75), rep("rare", 25))

# shuffle rows
commSample2 <- data.frame(commSample, commIdentity)
commSample3 <- commSample2[sample(nrow(commSample2)), ]

# get x and y coordinates for the grid map
xVector <- rep(1:20, n = 20)
yVector <- rep(1:20, n = 20)

# expand this into a grid
grid1 <- expand.grid(xVector, yVector)

# preData
preDat <- cbind(grid1, commSample3)
head(preDat)
names(preDat) <- c("xVal", "yVal", "preVec", "abundLevel")
preDat$preVecF <- factor(preDat$preVec)
head(preDat)
unique(preDat$preVecF)
with(preDat, table(preVecF))

plotSpecs <- theme_bw() + 
	theme(legend.position = "none") + 
	theme(axis.title = element_blank()) + 
	theme(axis.ticks = element_blank()) +
	theme(axis.text = element_blank()) + 
	theme(panel.grid = element_blank())

prePlot <- ggplot(data = preDat, 
	aes(xVal, yVal, color = preVecF)) + 
	geom_point() + plotSpecs

prePlot + ggtitle("Pre-disturbance") +
	theme(plot.title = element_text(lineheight = 3, face = "bold"))

# create a matrix
preMatrix <- matrix(preDat$preVec, 20, 20)
preMatrix
sppPool <- length(unique(preDat$preVecF))
sppPool

##############################################
##############################################
##############################################
# Apply a random disturbance that removes individuals from cells
head(preDat)
head(preMatrix)

preVector <- preDat$preVecF

sample.i <- sample(1:length(preVector), 100, replace = FALSE)
sample.i
preVector.i <- preVector

preVector.i[sample.i] <- NA
preVector.i

preDat2 <- cbind(preDat, preVector.i)
dim(preDat2)
head(preDat2)
preDat2

removalDat <- preDat2[is.na(preDat2$preVector.i), ]
removalDat
dim(removalDat)

# plot original, overlay with removed points
randOverlay <- geom_point(aes(xVal, yVal), data = removalDat, shape = 15, color = "white", size = 3) 

prePlot
prePlot + randOverlay + ggtitle("Post-disturbance")

##############################################
##############################################
# Calculating metrics
##############################################
##############################################
# now calculate metrics from preDat
head(preDat)
str(preDat)
preDat$cell <- 1:400
preMatrix

##############
### Scale - All individuals
##############

# calculate richness
head(preDat)
str(preDat)

# create a species x abund matrix (one row)
sppSiteBlank <- matrix(data = NA, nrow = 1, ncol = sppPool)
sppSiteDF <- data.frame(sppSiteBlank)
names(sppSiteDF) <- c(1:sppPool)
blankSppSite <- sppSiteDF

with(preDat, table(preVecF))
preSppCounts <- ddply(preDat, .(preVecF), summarise, freq = length(preVec), .drop = FALSE)
preSppCounts
preDF <- blankSppSite
preDF[1, ] <- preSppCounts$freq
preDF

preRenyi <- renyi(preDF, scales = c(0, 1, 2), hill = TRUE)
preRenyi
##############
### Scale - sample 25/400 individuals (top left corner, 1st 100 individuals)
##############
subDim <- 1:10
preMatSub <- preMatrix[subDim, subDim]
preMatSub
c(preMatSub)
preVectorSub <- c(preMatSub)
length(unique(preVectorSub))

xVector <- rep(subDim, n = length(subDim))
yVector <- rep(subDim, n = length(subDim))

# expand this into a grid
grid1 <- expand.grid(xVector, yVector)
dim(grid1)
grid1

# preData
preDatSub <- cbind(grid1, preVectorSub, factor(preVectorSub))

names(preDatSub) <- c("xVal", "yVal", "preVec", "preVecF")
preDatSub

preSppCounts <- ddply(preDatSub, .(preVecF), summarise, freq = length(preVec), .drop = FALSE)
preSppCounts
dim(preSppCounts)[1]

# create a species x abund matrix (one row)
sppSiteBlank <- matrix(data = NA, nrow = 1, ncol = dim(preSppCounts)[1])
sppSiteBlank
sppSiteDF <- data.frame(sppSiteBlank)
names(sppSiteDF) <- preSppCounts$preVecF
blankSppSite <- sppSiteDF
blankSppSite

preDF <- blankSppSite
preDF[1, ] <- preSppCounts$freq
preDF

rarecurve(preDF)
radfit(preDF)

preRenyiSub <- renyi(x = preDF, scales = c(0, 1, 2), hill = TRUE)
preRenyiSub
preRenyi

##############################################
##############################################
# Running the random removal simulation for the global population
##############################################
##############################################

head(preDat)
dim(preDat)
head(preMatrix)
sppPool <- length(unique(preDat$preVecF))
sppPool

preVector <- preDat$preVecF
length(preVector)
# write a function for random removal
preVector
preDat$preVecF

removeRandomF <- function(x, n) {
	sample.i <- sample(1:length(x), n, replace = FALSE)
	preVector.i <- preVector
	preVector.i[sample.i] <- NA
	preDat2 <- cbind(preDat, preVector.i)
	sppCount.i <- ddply(preDat2, .(preVector.i), summarise, freq = length(preVector.i), .drop = FALSE)[1:sppPool, ]
	# create a species x abund matrix (one row)
	sppSiteBlank <- matrix(data = NA, nrow = 1, ncol = dim(sppCount.i)[1])
	sppSiteDF <- data.frame(sppSiteBlank)
	names(sppSiteDF) <- sppCount.i$preVecF
	blankSppSite <- sppSiteDF
	sppAbund <- sppSiteBlank
	sppAbund[1, ] <- sppCount.i$freq	
	return(renyi(x = sppAbund, scales = c(0, 1, 2), hill = T))
}

removeRandomF(preVector, 100)

# Simulate 100 times
randomDistDF <- replicate(n = 100, removeRandomF(preVector, distSize))
randomResults <- data.frame(unname(t(randomDistDF)))
randomResults$scale <- rep("global", dim(randomResults)[1])
randomResults$pool <- rep(paste(sppPool, "species"), dim(randomResults)[1])
randomResults
str(randomResults)

##############################################
##############################################
# Running the random removal simulation for the local population
##############################################
##############################################

head(preDat)
dim(preDat)
head(preMatrix)

preVector <- preDat$preVecF
length(preVector)
# write a function for random removal
preVector
preDat$preVecF

removeRandomFsub <- function(x, n) {
	sample.i <- sample(1:length(x), n, replace = FALSE)
	postVector.i <- preVector
	postVector.i[sample.i] <- NA
	preDat2 <- cbind(preDat, postVector.i)
	
	# now create matrix for subsampling
	postMatrix.i <- matrix(preDat2$postVector.i, 20, 20)
	postMatrix.i
	postMatrixSub <- postMatrix.i[subDim, subDim]
	postVector <- c(postMatrixSub)
	xVector <- rep(subDim, n = length(subDim))
	yVector <- rep(subDim, n = length(subDim))

	# expand this into a grid
	grid1 <- expand.grid(xVector, yVector)
	
	# postDatSub
	postDatSub <- cbind(grid1, postVector)
	names(postDatSub) <- c("xVal", "yVal", "postVector")

	
	sppCount.i <- ddply(postDatSub, .(postVector), summarise, freq = length(postVector), .drop = FALSE)
	sppCount.i <- sppCount.i[-nrow(sppCount.i), ]
	
	# create a species x abund matrix (one row)
	sppSiteBlank <- matrix(data = NA, nrow = 1, ncol = dim(sppCount.i)[1])
	sppSiteDF <- data.frame(sppSiteBlank)
	names(sppSiteDF) <- sppCount.i$preVecF
	blankSppSite <- sppSiteDF
	sppAbund <- sppSiteBlank
	sppAbund[1, ] <- sppCount.i$freq	
	return(renyi(x = sppAbund, scales = c(0, 1, 2), hill = T))
}

removeRandomFsub(preVector, 100)

# Simulate 100 times, and calculate the hill numbers
randomDistDFsub <- replicate(n = 100, removeRandomFsub(preVector, distSize))
randomDistDFsub
randomResultsSub <- data.frame(unname(t(randomDistDFsub)))
randomResultsSub
str(randomResultsSub)
randomResultsSub$scale <- rep("local", dim(randomResults)[1])
randomResultsSub$pool <- rep(paste(sppPool, "species"), dim(randomResults)[1])
str(randomResultsSub)

####### complete data
randomDat <- rbind(randomResults, randomResultsSub)
dim(randomDat)
summary(randomDat)
randomDat$time <- rep("post", dim(randomDat)[1])
head(randomDat)

#
crap <- data.frame(preRenyi[1:3])
crap
preRenyiDF <- data.frame(t(crap))
preRenyiDF$scale <- "global"
preRenyiDF$pool <- paste(sppPool, "species")
preRenyiDF$time <- "pre"

crap <- data.frame(preRenyiSub[1:3])
preRenyiDFsub <- data.frame(t(crap))
preRenyiDFsub
preRenyiDFsub$scale <- "local"
preRenyiDFsub$pool <- paste(sppPool, "species")
preRenyiDFsub$time <- "pre"
#
preDFcomplete <- rbind(preRenyiDF, preRenyiDFsub)
preDFcomplete
names(preDFcomplete)[1:3] <- c("X1", "X2", "X3")

#########
spp20dat <- rbind(randomDat, preDFcomplete)
