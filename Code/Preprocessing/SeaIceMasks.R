# GISTEMP Uncertainty Analysis
# Version 1.2.1 (December 12, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# get all of the MERRA sea ice netcdf files
files <- system(sprintf('ls %s/MERRA/seaIce',ddirRaw),intern=TRUE)

# open the constant field to pull in grid info and make the land mask
handle <- nc_open(sprintf("%s/MERRA/constant/MERRA2_101.const_2d_asm_Nx.00000000.nc4",ddirRaw))

lon  <- ncvar_get(handle,'lon')
nlon <- length(lon)
lat  <- ncvar_get(handle,'lat')
nlat <- length(lat)
nt   <- length(files)

landField  <- ncvar_get(handle,"FRLAND")
iceField   <- ncvar_get(handle,"FRLANDICE")
lakeField  <- ncvar_get(handle,"FRLAKE")

propCutoff <- 0
landMask <- ifelse(landField > propCutoff | iceField > propCutoff | lakeField > 0, 1, NA)

nc_close(handle)

# load the sea ice fraction data into an R array
seaIceFrac <- array(NA,dim=c(nlon,nlat,nt))

for(i in 1:nt){
	handle <- nc_open(sprintf('%s/MERRA/seaIce/%s',ddirRaw,files[i]))

	seaIceFrac[,,i] <- ncvar_get(handle,'FRSEAICE')	

	nc_close(handle)
}

# build some time inds so we don't have to keep this data at the same state at the reanalyses
tYear <- 1980:2018
tMonth <- rep(1:12,length(tYear))

# update the old land mask with a land mask that evolves over time and includes the sea ice
# as 'land' rather than ocean
seasonalLandMask <- array(NA, dim=c(nlon,nlat,12))

for(i in 1:12){
	subInds <- which(tMonth==i)
	tempSeaIceMat <- apply(seaIceFrac[,,subInds]>0.2,c(1,2),any)
	seasonalLandMask[,,i] <- landMask | tempSeaIceMat
}


# make a maximal mask for more efficent interpolation with legacy code
maximalLandMask <- ifelse(apply(seasonalLandMask,c(1,2),function(x) any(!is.na(x))),1,NA)

# save the output
save(seasonalLandMask,maximalLandMask,file=sprintf('Data/%s/landMasks_%s.Rda','MERRA','MERRA'))


# create the sea ice land masks on the other reanalysis grids (function saves the output)
landMaskNewGrid('ERA',lon,lat,seasonalLandMask)
landMaskNewGrid('JRA',lon,lat,seasonalLandMask)
