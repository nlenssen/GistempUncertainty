wdOld <- '/Users/lenssen/Dropbox/GistempUncertainty/'
wdNew <- '/Users/lenssen/Dropbox/GistempUncertaintyFIXEDMEAN/'

plotdir <- '/Users/lenssen/Dropbox/GistempUncertaintyFIXEDMEAN/Figures/Corrections'



################################################################################
# Process the Old data
################################################################################

# write output to the paper figures directory
setwd(wdOld)
cexScale <- 1.5
gridlwd <- 1.5

# source the area calculation
load(sprintf('Data/%s/area_%s.Rda',reanalysis,reanalysis))

# load in the lst and sst raw result objects
load(sprintf('Data/%s/meanTimeSeries_%s.Rda',reanalysis,reanalysis))

# load in the limiting case for lsat sampling
load(sprintf('Data/%s/landMasks_%s.Rda',reanalysis,reanalysis))


# Load in external datasets
totalHomog <- read.csv('Data/Shared/GHCN/total-homog-uncertainty.csv',header=FALSE)

paramHomog <- read.csv('Data/Shared/GHCN/parametric-uncertainty.csv',header=FALSE)

gisLandOcean <- read.csv('Data/Shared/GISTEMP/gisLandOcean.csv',skip=1,header=TRUE)[-138,]
gisGlobal <- read.csv('Data/Shared/GISTEMP/gisGlobal.csv',skip=2,header=TRUE)[-138,]

gisGlobal2018 <- read.csv('Data/Shared/GISTEMP/gistemp2018.csv',skip=0,header=FALSE)[,]

# all of the time indicies used in plotting
tYear <- 1854:2014
tDec  <- seq(1850, 2010, by=10)
tYearFull <- 1850:2017
tYearFinal <- 1880:2016
tYear2018 <- 1880:2018

# decades to be used when plotting a limited range
subDec <- 4:17


# Final processing of our results from the analysis

# get the SST results (need to put into a function maybe?)
ensemblesS <- system('ls Data/Shared/ERSST', intern=TRUE)

# First, we compute sample mean and sample median annual estimates
annMatS  <- extratTSRdat(ensemblesS,'Data/Shared/ERSST',tYear)

# get the two different land uncertainties 
globalLandCI <- totalLandUncertainty(resultsList,tDec)

# Global Uncertainties:

# land calcs (empirical estimates) [take only through 2014]
sigma2L   <- decadeToYear(globalLandCI$diffVar,tYear,tDec)[1:length(tYear)]

# land variance calcs running over all 1880-2016
sigma2LFinal   <- decadeToYear(globalLandCI$diffVar,tYearFull,tDec)[tYearFull %in% tYearFinal]

# total land CI calculation (adding in homogonization uncertainty)
homogVar        <- ((totalHomog[,3] - totalHomog[,2])/(2*1.96))^2
totalLandVar    <- sigma2LFinal + homogVar

# ocean variance calc 
sigma2S  <- apply(annMatS,1,var)

# ocean variance calcs running over all 1880-2016
sigma2SFinal <- c(sigma2S[tYear %in% tYearFinal], rep(sigma2S[length(sigma2S)],2))


################################################################################
# Process the new data
################################################################################

setwd(wdNew)

# write output to the paper figures directory
cexScale <- 1.5
gridlwd <- 1.5

# source the area calculation
load(sprintf('Data/%s/area_%s.Rda',reanalysis,reanalysis))

# load in the lst and sst raw result objects
load(sprintf('Data/%s/meanTimeSeries_%s.Rda',reanalysis,reanalysis))

# load in the limiting case for lsat sampling
load(sprintf('Data/%s/landMasks_%s.Rda',reanalysis,reanalysis))


# decades to be used when plotting a limited range
subDec <- 4:17

# Final processing of our results from the analysis


# get the two different land uncertainties 
globalLandCINEW  <- totalLandUncertainty(resultsList,tDec)

# Global Uncertainties:

# land calcs (empirical estimates) [take only through 2014]
sigma2LNEW   <- decadeToYear(globalLandCINEW$diffVar,tYear,tDec)[1:length(tYear)]

# land variance calcs running over all 1880-2016
sigma2LFinalNEW    <- decadeToYear(globalLandCINEW$diffVar,tYearFull,tDec)[tYearFull %in% tYearFinal]

# total land CI calculation (adding in homogonization uncertainty)
homogVar        <- ((totalHomog[,3] - totalHomog[,2])/(2*1.96))^2
totalLandVarNEW    <- sigma2LFinalNEW + homogVar


################################################################################
# show change in the global sampling uncertainty
################################################################################
inds <- which(tYearFull %in% tYearFinal)

setwd(wdOld)
load(sprintf('Data/%s/meanTimeSeries_%s.Rda','ERA','ERA'))
globalLandCIERA <- totalLandUncertainty(resultsList,tDec)
sigmaLERA     <- decadeToYear(globalLandCIERA$diffVar,tYearFull,tDec)
totalLERA     <- sigmaLERA[inds] + homogVar

sciERA <- sqrt(sigmaLERA)* 1.96
tciERA <- sqrt(totalLERA) * 1.96

setwd(wdNew)
load(sprintf('Data/%s/meanTimeSeries_%s.Rda','ERA','ERA'))
globalLandCIERANEW <- totalLandUncertainty(resultsList,tDec)
sigmaLERANEW    <- decadeToYear(globalLandCIERANEW$diffVar,tYearFull,tDec)
totalLERANEW     <- sigmaLERANEW[inds] + homogVar

sciERANEW <- sqrt(sigmaLERANEW)* 1.96
tciERANEW <- sqrt(totalLERANEW) * 1.96


# make the plot
yr <- c(0,max(sciERANEW,sciERA))

pdf(sprintf('%s/globalSamplingCorrection.pdf',plotdir),10,7)
par(mfrow=c(1,1), mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFull, sciERANEW, ylim=yr,type='l', lwd=2,
	xlab='Year', ylab='95% LSAT Confidence Interval (°C)', main = 'Correction of Global LSAT Uncertainty', 
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,col='black')

points(tYearFull,sciERA,type='l', lwd=2,col='red')
grid(lwd=gridlwd)

legend('bottomleft',c('Correction', 'Original'),
	col=rep(c('black','red'),1),lwd=2,bg='white',cex=1.25)
dev.off()

################################################################################
# show change in the SH sampling uncertainty
################################################################################

# old 
setwd(wdOld)
load(sprintf('Data/%s/meanTimeSeries_%s.Rda','ERA','ERA'))
tempResultsList <- resultsList

nhTrue <- monthToYearMeans(tempResultsList$true$nh)
shTrue <- monthToYearMeans(tempResultsList$true$sh)

nDec <- length(tDec)

betaHemi <- matrix(NA,nrow=nDec,ncol=2)
sigmaHemi <- matrix(NA,nrow=nDec,ncol=2)
rawSDHemi <- matrix(NA,nrow=nDec,ncol=2)


for(i in 1:nDec){

	nhTemp <- monthToYearMeans(tempResultsList$mask[[i]]$nh)
	shTemp <- monthToYearMeans(tempResultsList$mask[[i]]$sh)

	#NH calc
	tempLM <- lm(nhTrue ~ nhTemp)
	betaHemi[i,1] <- tempLM$coefficients[2]
	sigmaHemi[i,1] <- summary(tempLM)$sigma
	rawSDHemi[i,1] <- sd(nhTemp - nhTrue)

	#SH calc
	tempLM <- lm(shTrue ~ shTemp)
	betaHemi[i,2] <- tempLM$coefficients[2]
	sigmaHemi[i,2] <- summary(tempLM)$sigma
	rawSDHemi[i,2] <- sd(shTemp - shTrue)
}


sdnL <- decadeToYear(rawSDHemi[,1],tYearFull,tDec)
sdsL <- decadeToYear(rawSDHemi[,2],tYearFull,tDec)

# new
setwd(wdNew)
load(sprintf('Data/%s/meanTimeSeries_%s.Rda','ERA','ERA'))

tempResultsList <- resultsList

nhTrue <- monthToYearMeans(tempResultsList$true$nh)
shTrue <- monthToYearMeans(tempResultsList$true$sh)

nDec <- length(tDec)

betaHemi <- matrix(NA,nrow=nDec,ncol=2)
sigmaHemi <- matrix(NA,nrow=nDec,ncol=2)
rawSDHemi <- matrix(NA,nrow=nDec,ncol=2)


for(i in 1:nDec){

	nhTemp <- monthToYearMeans(tempResultsList$mask[[i]]$nh)
	shTemp <- monthToYearMeans(tempResultsList$mask[[i]]$sh)

	#NH calc
	tempLM <- lm(nhTrue ~ nhTemp)
	betaHemi[i,1] <- tempLM$coefficients[2]
	sigmaHemi[i,1] <- summary(tempLM)$sigma
	rawSDHemi[i,1] <- sd(nhTemp - nhTrue)

	#SH calc
	tempLM <- lm(shTrue ~ shTemp)
	betaHemi[i,2] <- tempLM$coefficients[2]
	sigmaHemi[i,2] <- summary(tempLM)$sigma
	rawSDHemi[i,2] <- sd(shTemp - shTrue)
}


sdnLNEW <- decadeToYear(rawSDHemi[,1],tYearFull,tDec)
sdsLNEW <- decadeToYear(rawSDHemi[,2],tYearFull,tDec)


# check to make sure NH doesn't change
all(sdnL == sdnLNEW)

# make the plot
pdf(sprintf('%s/shSamplingCorrection.pdf',plotdir),10,7)

par(mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFull,sdsL*1.96,type='l',lwd=2,col = 'red',
	ylim = c(0,max(sdsL*1.96,sdsLNEW*1.96)),
	xlab = 'Year',
	ylab = 'S. Hemisphere LSAT 95% Sampling CI (°C)', 
	main= 'Corrected SH Sampling Uncertainty',
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
points(tYearFull,sdsLNEW*1.96,type='l',lwd=2,col = 'black')
grid(lwd=gridlwd)

legend('bottomleft',c('Corrected','Original'),
	col=rep(c('black','red','blue'),1),lwd=2,bg='white',cex=1.25)

dev.off()

################################################################################
# show change in the mean comparison
################################################################################
# old
setwd(wdOld)
load(sprintf('Data/%s/simpleMean_%s.Rda',reanalysis,reanalysis))
load(sprintf('Data/%s/limitingGlobalCI_%s.Rda',reanalysis,reanalysis))

limitingCI <- globalLimitCI

globalLandCISimple <- totalLandUncertainty(resultsListSimple,tDec)
x <- tYearFull
y1 <- decadeToYear(sqrt(globalLandCI$diffVar)*1.96,tYearFull,tDec)
y2 <- decadeToYear(sqrt(globalLandCISimple$diffVar)*1.96,tYearFull,tDec)


# new
setwd(wdNew)
load(sprintf('Data/%s/simpleMean_%s.Rda',reanalysis,reanalysis))
load(sprintf('Data/%s/limitingGlobalCI_%s.Rda',reanalysis,reanalysis))

limitingCINEW <- globalLimitCI


globalLandCISimple <- totalLandUncertainty(resultsListSimple,tDec)
x <- tYearFull
y1NEW <- decadeToYear(sqrt(globalLandCINEW$diffVar)*1.96,tYearFull,tDec)


# make the plot comparing the change in gistemp sampling uncertainty

pdf(sprintf('%s/meanCompCorrection.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)

plot(x,y1,ylim=c(0,max(y1,y2,y1NEW)),
	type='l',col='red',lwd=2, 
	xlab='Year',ylab='95% LSAT Sampling Confidence Interval',
	cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=1.5)
points(x,y1NEW,type='l',col='black',lwd=2)

points(x,y2,type='l',col='green2',lwd=2)


abline(v=1880,lty=3,lwd=2)

abline(h=limitingCI,col='blue',lwd=2,lty=3)
abline(h=limitingCINEW,col='blue',lwd=2,lty=2)
points(x,y1,type='l',col='red',lwd=2)

legend('topright', c('GISTEMP Correction', 'GISTEMP Original','Simple Mean','Corrected Limiting Uncertainty', 'Orig. Limiting Uncertainty'), 
	col=c('black', 'red','green2', 'blue', 'blue'),lty=c(1,1,1,2,3),lwd=2,bg='white',cex=cexScale)

dev.off()

################################################################################
# show change in the total global uncertainty
################################################################################

# old
lstsstFinal <- sqrt(AL^2 * totalLandVar + AS^2 * sigma2SFinal)*1.96


# new
lstsstFinalNEW <- sqrt(AL^2 * totalLandVarNEW + AS^2 * sigma2SFinal)*1.96

pdf(sprintf('%s/globalCICorrection.pdf',plotdir),10,7)

par(mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFinal,lstsstFinal*1.96,type='l',lwd=2,col = 'red',
	ylim = c(0,max(lstsstFinal*1.96,lstsstFinalNEW*1.96)),
	xlab = 'Year',
	ylab = 'Global Total 95% CI (°C)', 
	main= 'Corrected Total Global Uncertainty',
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
points(tYearFinal,lstsstFinalNEW*1.96,type='l',lwd=2,col = 'black')
grid(lwd=gridlwd)

legend('bottomleft',c('Correction','Original'),
	col=rep(c('black','red'),1),lwd=2,bg='white',cex=1.25)

dev.off()
