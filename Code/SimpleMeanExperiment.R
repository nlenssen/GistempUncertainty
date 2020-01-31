# Performs the uncertainty analysis using a simple area-weighted mean instead
# of the GISTEMP mean.

# GISTEMP Uncertainty Analysis
# Version 1.2.1 (December 12, 2019)
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
					   zoneMask = zoneMask,
					   simpleMean=TRUE,
					   ALn = ALn, ALs = ALs)

rm(anomalyFieldAdj)
gc()

maskList <- list()


for(i in 1:length(files)){
	cat(paste('Starting Decade', i, '\n'))
	load(sprintf('Data/%s/InterpolatedFields/%s',reanalysis,files[i]))


	finalField <- array(NA, dim=dim(interpolatedField))
	
	for(j in 1:length(tMonth)){
		finalField[,,j] <- interpolatedField[,,j] * seasonalLandMask[,,tMonth[j]]
	}

	rm(interpolatedField)

	maskMean <- globalMean(anomalyField = finalField,
						   mask = NULL,
						   lat = lat,
						   nCores=nCores,
						   zoneMask = zoneMask,
						   simpleMean = TRUE,
						   ALn = ALn, ALs = ALs)

	gc()

	maskList[[i]] <- maskMean

}

resultsListSimple <- list(true = trueMean, mask = maskList)

save(resultsListSimple,file=sprintf('Data/%s/simpleMean_%s.Rda',reanalysis,reanalysis))