# Performs the limiting uncertainty experiment. Warning: very long running.

# GISTEMP Uncertainty Analysis
# Version 1.1 (August 15, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# load in the necessary data
load(sprintf('Data/%s/anomalyData_%s.Rda',reanalysis,reanalysis))
load(sprintf('Data/%s/landMasks_%s.Rda',reanalysis,reanalysis))


# run the interpolation
interpolatedFieldLimit <- interpolateFieldPointwise(
	maximalLandMask,radius,nCores,anomalyData,maximalLandMask)ß

# save the output
save(interpolatedFieldLimit,file=sprintf("%s/interpolatedFieldMaximal.Rda",ddir))

gc()


# calculate the means (with option to load in the interpolated field)
load(sprintf("%s/era5Data.Rda",ddir))
load(sprintf("%s/seasonalLandMask.Rda",ddir))

# a helper load if the interpolation is skipped
if(!exists('interpolatedFieldLimit')) load(sprintf("%s/interpolatedFieldMaximal.Rda",ddir))

lon <- anomalyData$lon
lat <- anomalyData$lat
anomalyField <- anomalyData$anomalyField
zoneMask <- anomalyData$zoneMask
rm(anomalyData)

tMonth <- rep(1:12,dim(anomalyField)[3]/12)

# properly mask the "true" field
anomalyFieldAdj      <- array(NA, dim=dim(anomalyField))
interpolatedFieldAdj <- array(NA, dim=dim(anomalyField))

for(j in 1:length(tMonth)){
	anomalyFieldAdj[,,j] <- anomalyField[,,j] * seasonalLandMask[,,tMonth[j]]
	interpolatedFieldAdj[,,j] <- interpolatedFieldLimit[,,j] * seasonalLandMask[,,tMonth[j]]
}

rm(anomalyField,interpolatedFieldLimit)


trueMean <- globalMean(anomalyField = anomalyFieldAdj,
						   mask = NULL,
						   lat = lat,
						   nCores=nCores,
						   zoneMask = zoneMask)


interpMean <- globalMean(anomalyField = interpolatedFieldAdj,
						   mask = NULL,
						   lat = lat,
						   nCores=nCores,
						   zoneMask = zoneMask)



globalVar <- var(trueMean$global - interpMean$global)
globalLimitCI <- sqrt(globalVar)*1.96

regFit     <- lm(trueMean$global ~ interpMean$global)
globalBias <- regFit$coefficients[2]
regCI      <- summary(regFit)$coefficients[2,2]

regFitCorrected <- lm(interpMean$global ~ trueMean$global)
globalBiasCorrected <- regFitCorrected$coefficients[2]
regCICorrected      <- summary(regFitCorrected)$coefficients[2,2]

c(globalBiasCorrected - 1.96 * regCICorrected, globalBiasCorrected + 1.96 * regCICorrected)

save(trueMean,interpMean,globalLimitCI,globalBias,regCI,
	file=sprintf('Data/%s/limitingGlobalCI_%s.Rda',reanalysis,reanalysis))
