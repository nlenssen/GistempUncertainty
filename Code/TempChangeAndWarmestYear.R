# Performs the Pr(warmest year) Monte Carlo simulation.

# GISTEMP Uncertainty Analysis
# Version 1.2.1 (December 12, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# load in the time series and CI object
load(sprintf('Data/%s/totalCI_%s.Rda',reanalysis,reanalysis))

# grab a meanSeries and a ciSeries to write the function that allows us
# to do a MC analysis of possible slopes

# do the warmest year calc for 2014 and 2015
yearInds <- which(gisGlobal2018[,1] < 2015)
results2014 <- resultsMC(gisGlobal2018[yearInds,2],lstsstFinal[yearInds],gisGlobal2018[yearInds,1],alpha=0.5,M=10000)

yearInds <- which(gisGlobal2018[,1] < 2016)
results2015 <- resultsMC(gisGlobal2018[yearInds,2],lstsstFinal[yearInds],gisGlobal2018[yearInds,1],alpha=0.5,M=10000)

resultsPresent <- resultsMC(gisGlobal2018[,2],lstsstFinal[],gisGlobal2018[,1],alpha=0.5,M=10000)


# calculate the probability of warmest year
round(results2014$ind*100 ,3)
round(results2015$ind *100,3)
round(resultsPresent$ind *100,3)

round(results2014$ar *100,3)
round(results2015$ar*100,3)
round(resultsPresent$ar *100,3)