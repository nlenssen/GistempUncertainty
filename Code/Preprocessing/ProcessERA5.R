# CDO command used to map the ERA data to the MERRA grid

# cdo remapcon,r576x361 era5_t2m.nc era5_t2m_merraGrid.nc

###############################################################################
# STEP (1) Load in R objects (R data and functions)
###############################################################################
tYear <- 1979:endYear
nt    <- length(tYear)*12

###############################################################################
# STEP (2) Load all the data
###############################################################################
handle <- nc_open(sprintf('%s/ERA/era5_t2m_merraGrid.nc',ddirRaw))

rawLon <- ncvar_get(handle,'lon')
nlon <- length(rawLon)
lon <- c(rawLon[(nlon/2+1):nlon]-360,rawLon[1:(nlon/2)])

rawLat <- ncvar_get(handle,'lat')
nlat <- length(rawLat)
lat <- rawLat


rawTemp <- ncvar_get(handle,'t2m',start=c(1,1,1),count=c(-1,-1,nt))

nc_close(handle)

tempArray <- rawTemp[c((nlon/2+1):nlon,1:(nlon/2)),,]
rm(rawTemp)
###############################################################################
# STEP (3) Compute Anomalies on the grid
###############################################################################
anomalyField <- array(NA, dim=dim(tempArray))

monthVector <- rep(1:12,length(tYear))

monthNorm <- array(NA,dim=c(dim(tempArray)[1:2],12))

pb   <- txtProgressBar(1, 12, style=3)

for(month in 1:12){
	setTxtProgressBar(pb, month)

	inds <- which(monthVector==month)
	monthNorm[,,month] <- apply(tempArray[,,inds],c(1,2),mean)

	for(j in 1:length(inds)){	
		anomalyField[,,inds[j]] <- tempArray[,,inds[j]] - monthNorm[,,month]
	}
}

###############################################################################
# STEP (4) Create Zone Mask
###############################################################################

zoneMask <- matrix(NA, nrow=length(lon), ncol=length(lat))

for(i in 1:length(lon)){
	for(j in 1:length(lat)){
		zoneMask[i,j] <- zoneAssign(lon[i], lat[j])
	}
}

anomalyData <- list(anomalyField = anomalyField,
				 lon = lon,
				 lat = lat,
				 zoneMask = zoneMask)


save(anomalyData,file=sprintf('Data/%s/anomalyData_%s.Rda','ERA','ERA'))

