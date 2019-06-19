files <- system(sprintf('ls %s/JRA/t2m/',ddirRaw), intern=TRUE)
nt <- length(files)

###############################################################################
# STEP (1) Get grid info and create land mask
###############################################################################

handle <- nc_open(sprintf('%s/JRA/t2m/%s',ddirRaw,files[1]))

rawLon <- ncvar_get(handle,'g4_lon_2')
nlon <- length(rawLon)
lon <- c(rawLon[(nlon/2+1):nlon]-360,rawLon[1:(nlon/2)])

rawLat <- ncvar_get(handle,'g4_lat_1')
nlat <- length(rawLat)
lat <- rev(rawLat)

nc_close(handle)


###############################################################################
# STEP (2b) Populate R array from JRA ncdfs
###############################################################################

tempArray <- array(NA, dim=c(nlon,nlat,12*nt))

for(i in 1:nt){
	handle <- nc_open(sprintf('%s/JRA/t2m/%s',ddirRaw,files[i]))

	rawTemp <- ncvar_get(handle,'TMP_GDS4_HTGL_S123')
	tempArray[,,(1+(i-1)*12):(i*12)] <- rawTemp[c((nlon/2+1):nlon,1:(nlon/2)),rev(1:nlat),]


	nc_close(handle)
}

###############################################################################
# STEP (3) Compute Anomalies
###############################################################################
anomalyField <- array(NA, dim=dim(tempArray))
monthVector <- rep(1:12,nt)
monthNorm <- array(NA,dim=c(dim(tempArray)[1:2],12))

for(month in 1:12){
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

save(anomalyData,file=sprintf('Data/%s/anomalyData_%s.Rda','JRA','JRA'))


