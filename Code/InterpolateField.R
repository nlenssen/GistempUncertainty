# Interpolates the SAT anomaly field on the reanalysis grid given a coverage
# mask. This step is both memory and CPU intensive and may require tuning of
# the number of cores used in the parallel loop.

# GISTEMP Uncertainty Analysis
# Version 1.2.1 (December 12, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# load in the necessary data
load(sprintf('Data/%s/landMasks_%s.Rda',reanalysis,reanalysis))
load(sprintf('Data/%s/decadalMasks_%s.Rda',reanalysis,reanalysis))

# perform the analysis decade by decade
for(i in 1:dim(decadalMasks)[3]){
	cat(paste('Decade',i,'\n'))

	load(sprintf('Data/%s/anomalyData_%s.Rda',reanalysis,reanalysis))

	# unpack the list for more readable code
	lon <- anomalyData$lon
	nlon <- length(lon)
	lat <- anomalyData$lat
	nlat <- length(lat)
	nt <- dim(anomalyData$anomalyField)[3]

	# make the anomaly data a matrix with non-land removed
	tempAnomaly <- apply(anomalyData$anomalyField,3L,c)

	rm(anomalyData)

	goodInds    <- which(!is.na(maximalLandMask) | decadalMasks[,,i]==1)
	goodIndsMat <- which(!is.na(maximalLandMask) | decadalMasks[,,i]==1,arr.ind=T)

	anomalyMat <- tempAnomaly[goodInds,]

	# get the matrix locations of the tables
	stationIndsMat <- which(decadalMasks[,,i]==1,arr.ind=T)
	partialMask <- c(decadalMasks[,,i])[goodInds]
	stationInds <- which(partialMask==1)

	rm(tempAnomaly)
	

	# interpolate the field to the maximal land mask
	interpolatedField <- interpolateFieldPointwiseLowMem(
							radius=radius,
							lon=lon,
							lat=lat,
							goodIndsMat=goodIndsMat,
							stationInds=stationInds,
							stationIndsMat=stationIndsMat,
							anomalyMat=anomalyMat)


	# save the output each loop to allow restarts
	ofname <- sprintf('Data/%s/InterpolatedFields/interpolatedField%02d_%s.Rda',reanalysis,i,reanalysis)
	save(interpolatedField,file=ofname)
	
	# clean up to prevent memory issues
	rm(interpolatedField)
	gc()
}