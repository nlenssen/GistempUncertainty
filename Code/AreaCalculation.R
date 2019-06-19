# A helper script used to calculate the land area proportion for a grid given
# a land mask.

# GISTEMP Uncertainty Analysis
# Version 1.0 (May 1, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# load in the necessary data
load(sprintf('Data/%s/anomalyData_%s.Rda',reanalysis,reanalysis))
load(sprintf('Data/%s/landMasks_%s.Rda',reanalysis,reanalysis))

# unpack the list for readability
lon          <- anomalyData$lon
lat          <- anomalyData$lat
anomalyField <- anomalyData$anomalyField
zoneMask     <- anomalyData$zoneMask

cosMat <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180))

# calculate the land area
landArea <- sum(maximalLandMask * cosMat,na.rm=T)
totalArea <- sum(cosMat)

AL <- landArea/totalArea
AS <- 1-AL

# get the hemispheric proportions as well
nLats <- which(lat>=0)
sLats <- which(lat<=0)

ALn <- sum(maximalLandMask[,nLats] * cosMat[,nLats],na.rm=T)/sum(cosMat[,nLats])
ALs <- sum(maximalLandMask[,sLats] * cosMat[,sLats],na.rm=T)/sum(cosMat[,sLats])

ASn <- 1 - ALn
ASs <- 1 - ALs

save(AL,AS,ALn,ALs,ASn,ASs,file=sprintf('Data/%s/area_%s.Rda',reanalysis,reanalysis))