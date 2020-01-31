# Create coverage masks on the reanalyses grid and save to disk.

# GISTEMP Uncertainty Analysis
# Version 1.2.1 (December 12, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# load station data
load('Data/Shared/GHCN/stationData.Rda')

# Load reanalysis grid metadata
load(sprintf('Data/%s/anomalyData_%s.Rda',reanalysis,reanalysis))
lon <- anomalyData$lon
lat <- anomalyData$lat

rm(anomalyData)

startYears <- seq(startYear,2010,by=10)

decadalMasks <- array(NA, dim=c(length(lon),
								length(lat),
								length(startYears)))
coverageList <- list()

# loop over decades to determine coverage and map to the reanalysis grid
for(i in 1:length(startYears)){
	decade <- (startYears[i]-1):(startYears[i]+9)

	# start by subsetting the data that appears in this data
	partialInds <- which(stationData[,1] %in% decade)

	partial <- stationData[partialInds,]

	# replace all -9999 with NA
	partial[partial==-9999] <- NA

	# pull the unique station ids and figure out what years they 
	ids <- unique(dataID[partialInds])


	# Compute which stations have coverage in the decade
	coverageMatFull <- data.frame()

	for(j in 1:length(ids)){
		station <- subset(partial, dataID[partialInds] == ids[j])
		coverageMatFull <- rbind(coverageMatFull,
			coverageCheck(station,decade,ids[j]))
	}

	# Now take the stations with coverage and map to lon/lat
	coveredStations <- subset(stationID, stationID$id %in% 
						 coverageMatFull[coverageMatFull[,2]==1,][,1])

	coverageList[[i]] <- coveredStations

	# Now we need to fill out the merra grid based on these stations
	decadalMasks[,,i] <- merraGridMask(coveredStations,lon,lat)
}

# save the output
save(decadalMasks,file=sprintf('Data/%s/decadalMasks_%s.Rda',reanalysis,reanalysis))