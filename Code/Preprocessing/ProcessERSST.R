# GISTEMP Uncertainty Analysis
# Version 1.2.1 (December 12, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

library(ncdf4)

source('Code/Functions.R')

ddirRaw <- '/Users/lenssen/Documents/gistempRaw'


nCores <- 4

ftpLoc <- "ftp://ftp.ncdc.noaa.gov/pub/data/cmb/ersst/v4/ensemble/"

# make the temp dir to hold the downloaded files
system(sprintf("mkdir %s/fullERSST/",ddirRaw), ignore.stderr = TRUE)

# loop through the 1000 ensemble members of the ERSST ensemble
for(ensembleNumber in 967:1000){
	cat(paste('Starting Ensemble Member',ensembleNumber,'\n'))

	fileName <- sprintf("sst2d.ano.1854.2014.ensemble.%04d.dat",ensembleNumber)

	# download the data to a local directory (temporarily, removed in loop)
	if(!file.exists(sprintf("%s/fullERSST/%s.gz",ddirRaw,fileName)) &
	   !file.exists(sprintf("%s/fullERSST/%s",   ddirRaw,fileName)) ) {		
		system(sprintf("wget -t 10 -P %s/fullERSST %s/%s.gz",ddirRaw,ftpLoc,fileName), ignore.stdout = TRUE)
	}
	
	try(system(sprintf("gunzip %s/fullERSST/%s.gz",ddirRaw,fileName), ignore.stderr = TRUE))
	
	# load in the data
	dat <- readGridded(sprintf("%s/fullERSST",ddirRaw),sprintf("sst2d.ano.1854.2014.ensemble.%04d.dat",ensembleNumber))

	lon <- dat$lon
	lat <- dat$lat
	sst <- dat$sst

	
	# Run some calculations that only need to be taken once for the whole loop
	if(ensembleNumber == 1){
		# work with the lat lon grid to make a zone mask
		zoneMask <- matrix(NA, nrow=length(lon), ncol=length(lat))
	
		for(i in 1:length(lon)){
			for(j in 1:length(lat)){
				zoneMask[i,j] <- zoneAssign(lon[i], lat[j])
			}
		}

		# Calculate the sst areas in the NH and SH
		load(sprintf('Data/%s/landMasks_%s.Rda','MERRA','MERRA'))
		seasonalLandMaskRef <- seasonalLandMask
		
		# get the lon/lat for merra
		handle <- nc_open(sprintf("%s/MERRA/constant/MERRA2_101.const_2d_asm_Nx.00000000.nc4",ddirRaw))

		lonRef  <- ncvar_get(handle,'lon')
		latRef  <- ncvar_get(handle,'lat')

		nc_close(handle)

		# make the new mask on the 2 degree grid
		newMask <- array(NA, dim=c(length(lon),length(lat),12))

		for(sInd in 1:12){
			for(i in 1:length(lon)){
				for(j in 1:length(lat)){
					lonInd <- which.min(abs(lon[i] - lonRef))	
					latInd <- which.min(abs(lat[j] - latRef))

					newMask[i,j,sInd] <- seasonalLandMaskRef[lonInd,latInd,sInd]
				}
			}
		}
		maximalLandMask <- ifelse(apply(newMask,c(1,2),function(x) any(!is.na(x))),NA,1)


		# calculate the area
		cosMat <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180))

		# calculate the sea area 
		seaArea <- sum(maximalLandMask * cosMat,na.rm=T)
		totalArea <- sum(cosMat)

		AS2 <- seaArea/totalArea
		AL2 <- 1-AS2

		# get the hemispheric proportions as well
		nLats <- which(lat>=0)
		sLats <- which(lat<=0)

		ASn2 <- sum(maximalLandMask[,nLats] * cosMat[,nLats],na.rm=T)/sum(cosMat[,nLats])
		ASs2 <- sum(maximalLandMask[,sLats] * cosMat[,sLats],na.rm=T)/sum(cosMat[,sLats])

		ALn2 <- 1 - ASn2
		ALs2 <- 1 - ASs2
	}
	

	# take the global mean with just the ocean
	oceanMean <- globalMean(anomalyField = dat$sst,
						   mask = NULL,
						   lat = lat,
						   nCores=nCores,
						   zoneMask = zoneMask,
						   ALn = ASn2, ALs= ASs2)

	# save output
	ofname <- sprintf("sstMean%04d.Rda",ensembleNumber)
	save(oceanMean,file=sprintf("%s/ERSST/%s",ddirRaw,ofname))
}

