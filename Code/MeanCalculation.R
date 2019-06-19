# Calculated the global mean temperature anomaly using the GISTEMP equal area
# averaging method. This step is both memory and CPU intensive and may
# require tuning of the number of cores used in the parallel loop.

# GISTEMP Uncertainty Analysis
# Version 1.0 (May 1, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# load the data
files <- system(sprintf('ls Data/%s/InterpolatedFields/',reanalysis),intern=T)

load(sprintf('Data/%s/anomalyData_%s.Rda',reanalysis,reanalysis))
load(sprintf('Data/%s/landMasks_%s.Rda',reanalysis,reanalysis))


lon <- anomalyData$lon
lat <- anomalyData$lat
anomalyField <- anomalyData$anomalyField
zoneMask <- anomalyData$zoneMask

rm(anomalyData)

tMonth <- rep(1:12,dim(anomalyField)[3]/12)

# properly mask the 'true' field
anomalyFieldAdj <- array(NA, dim=dim(anomalyField))

for(j in 1:length(tMonth)){
	anomalyFieldAdj[,,j] <- anomalyField[,,j] * seasonalLandMask[,,tMonth[j]]
}


rm(anomalyField)
gc()

# True Mean
trueMean <- globalMean(anomalyField = anomalyFieldAdj,
					   mask = NULL,
					   lat = lat,
					   nCores=nCores,
					   zoneMask = zoneMask)


rm(anomalyFieldAdj)
gc()

maskList <- list()

# loop over the decades for the mask means
for(i in 1:length(files)){
	cat(paste('Starting Decade', i, '\n'))
	load(sprintf('Data/%s/InterpolatedFields/%s',reanalysis,files[i]))

	interpolatedFieldAdj <- array(NA, dim=dim(interpolatedField))
	
	for(j in 1:length(tMonth)){
		interpolatedFieldAdj[,,j] <- interpolatedField[,,j] * seasonalLandMask[,,tMonth[j]]
	}

	rm(interpolatedField)

	maskMean <- globalMean(anomalyField = interpolatedFieldAdj,
						   mask = NULL,
						   lat = lat,
						   nCores=nCores,
						   zoneMask = zoneMask)

	
	maskList[[i]] <- maskMean

	rm(maskMean)
	gc()
}

resultsList <- list(true = trueMean, mask = maskList)

save(resultsList,file=sprintf('Data/%s/meanTimeSeries_%s.Rda',reanalysis,reanalysis))