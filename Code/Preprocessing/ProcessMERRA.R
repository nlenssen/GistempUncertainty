# A script that reads in the raw MERRA Data and creates an R Data object
# of monthly SAT anomalies on the MERRA grid

# GISTEMP Uncertainty Analysis
# Version 1.2.1 (December 12, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

###############################################################################
# STEP (1) Populate R array from merra ncdf
###############################################################################
files <- system(sprintf('ls %s/MERRA/t2m',ddirRaw),intern=TRUE)

# Load one netcdf field to get the metadata for the creation of our
# data structures and loops and build out landmask

handle <- nc_open(sprintf("%s/MERRA/constant/MERRA2_101.const_2d_asm_Nx.00000000.nc4",ddirRaw))

lon  <- as.numeric(ncvar_get(handle,"lon"))
lat  <- as.numeric(ncvar_get(handle,"lat"))

nx <- length(lon)
ny <- length(lat)
nt <- length(files)

nc_close(handle)

# pull each monthly temperature field into our array
merraArray <- array(NA,dim=c(nx,ny,nt))

for(i in 1:nt){
	handle <- nc_open(sprintf('%s/MERRA/t2m/%s',ddirRaw,files[i]))
	merraArray[,,i] <- ncvar_get(handle,'T2M')
	nc_close(handle)
}


###############################################################################
# STEP (2) Compute Anomalies on the merra grid
###############################################################################
anomalyField <- array(NA, dim=dim(merraArray))
monthVector <- rep(1:12,nt/12)
monthNorm <- array(NA,dim=c(dim(merraArray)[1:2],12))


for(month in 1:12){
	inds <- which(monthVector==month)
	monthNorm[,,month] <- apply(merraArray[,,inds],c(1,2),mean)

	for(j in 1:length(inds)){	
		anomalyField[,,inds[j]] <- merraArray[,,inds[j]] - monthNorm[,,month]
	}
}

###############################################################################
# STEP (3) Create Zone Mask
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

save(anomalyData,file=sprintf('Data/%s/anomalyData_%s.Rda','MERRA','MERRA'))