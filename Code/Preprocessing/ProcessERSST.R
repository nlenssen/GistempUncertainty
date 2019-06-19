# loop through the 1000 ensemble members of the ERSST ensemble
for(ensembleNumber in 1:1000){

	# download the data to a local directory (temporarily, removed in loop)
	ftpLoc <- "ftp://ftp.ncdc.noaa.gov/pub/data/cmb/ersst/v4/ensemble/"
	fileName <- sprintf("sst2d.ano.1854.2014.ensemble.%04d.dat",ensembleNumber)
	system(sprintf("wget -t 10 -P %s/tempERSST %s/%s.gz",ddirRaw,ftpLoc,fileName),ignore.stdout = TRUE)

	system(sprintf("gunzip %s/tempERSST/%s.gz",ddirRaw,fileName))

	# load in the data

	dat <- readGridded(sprintf("%s/tempERSST",ddirRaw),sprintf("sst2d.ano.1854.2014.ensemble.%04d.dat",ensembleNumber))

	lon <- dat$lon
	lat <- dat$lat
	sst <- dat$sst

	# work with the lat lon grid to make a zone mask
	zoneMask <- matrix(NA, nrow=length(lon), ncol=length(lat))

	for(i in 1:length(lon)){
		for(j in 1:length(lat)){
			zoneMask[i,j] <- zoneAssign(lon[i], lat[j])
		}
	}


	# take the global mean with just the ocean
	oceanMean <- globalMean(anomalyField = dat$sst,
						   mask = NULL,
						   lat = lat,
						   nCores=nCores,
						   zoneMask = zoneMask)

	# save output
	ofname <- sprintf("sstMean%04d.Rda",ensembleNumber)
	save(oceanMean,file=sprintf("Data/Shared/ERRSST/%s",ofname))

	# remove the data
	system(sprintf("rm %s/tempERSST/%s",ddirRaw,fileName))
}

# remove the temp directory used to hold the raw files
system(sprintf("rm -r %s/tempERSST/",ddirRaw))
