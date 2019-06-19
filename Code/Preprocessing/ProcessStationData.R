########################################################################
# Import the GHCN station metadata into an easy to work with format
########################################################################
# Read in the (currently) relevant data from the station data

rawInfo <- readLines(sprintf('Data/Shared/GHCN/%s',ddirRaw,infoFile))

# Vectorized to run over all rows of the station file

id  <- substr(rawInfo,1,11) # as string this time...
lat <- as.numeric(substr(rawInfo,13,20))
lon <- as.numeric(substr(rawInfo,22,30))


# Place in a data frame to mirror the raw, human-readable organization
stationID <- data.frame(id,lon,lat, stringsAsFactors=FALSE)


########################################################################
# Import the GISTEMP station data into an easy to work with format
########################################################################
raw <- readLines(sprintf('%s/GHCN/%s',ddirRaw,dataFile))

dataID <- rep(NA,length(raw))
stationData <- matrix(NA,nrow=length(raw),ncol=13)
colnames(stationData) <- c('year', 1:12)

# helper function used to parse the data portion of the input string
# 11 or 12 digit station?

for(i in 1:length(raw)){
	temp <- raw[i]
	station <- as.character(substr(temp,1,12))
	year    <- as.numeric(substr(temp,13,16))

	data <- parseDataString(substr(temp,21,nchar(temp)))

	dataID[i] <- station
	stationData[i,] <- c(year,data)
}


# save the processed data
save(stationID, stationData, dataID,  file = 'Data/Shared/GHCN/stationData.Rda')