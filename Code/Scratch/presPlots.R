source('Namelists/Namelist_ERA.R')
# Performs final analysis steps and replicates all Figures in Lenssen et al
# uncertainty paper.

library(RColorBrewer)

plotdir <- '/Users/lenssen/Dropbox/DEES/IMSC19/gistempTalk/Figures/codeFigures'
# GISTEMP Uncertainty Analysis
# Version 1.0 (May 1, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# write output to the paper figures directory
cexScale <- 1.5
gridlwd <- 1.5

# source the area calculation
load(sprintf('Data/%s/area_%s.Rda',reanalysis,reanalysis))

# load in the lst and sst raw result objects
load(sprintf('Data/%s/meanTimeSeries_%s.Rda',reanalysis,reanalysis))

# load in the limiting case for lsat sampling
load(sprintf('Data/%s/landMasks_%s.Rda',reanalysis,reanalysis))

# load in the mask
load(sprintf('Data/%s/decadalMasks_%s.Rda',reanalysis,reanalysis))

# get the list of interploated fields
files <- system(sprintf('ls Data/%s/InterpolatedFields/',reanalysis),intern=T)

# load in the anomaly field and land mask
load(sprintf('Data/%s/anomalyData_%s.Rda',reanalysis,reanalysis))
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


################################################################################
# Series of plots demonstrating the interpolation process 1880s
################################################################################
dec <- 4
load(sprintf('Data/%s/InterpolatedFields/%s',reanalysis, files[dec]))


eraTimePoint <- 137

# build a time map to know which time point we are at
tYearTemp <- 1979:endYear
timeMap <- cbind(rep(tYearTemp,each=12),rep(1:12, length(tYearTemp)))

season <- eraTimePoint %% 12

trueField <- anomalyData$anomalyField[,,eraTimePoint]

landMask <- seasonalLandMask[,,season]

range(trueField)
# smoosh the anomaly colorbar down
zMax <- 6

zr <- c(-zMax, zMax)

anomClip <- trueField
anomClip[anomClip>zMax] <- zMax 
anomClip[anomClip< -zMax] <- -zMax

interpClip <- interpolatedField[,,eraTimePoint]
interpClip[interpClip>zMax] <- zMax
interpClip[interpClip< -zMax] <- -zMax


pal <- rev(designer.colors(256,brewer.pal(11,'RdBu')))

lon <- anomalyData$lon
lat <- anomalyData$lat

z <- decadalMasks[,,dec]
z[z==0] <- NA

gridList <- expand.grid(list(lon=lon,lat=lat))
stationLoc <- which(!is.na(z))

pointLoc <- gridList[stationLoc,]
anomLoc  <- c(trueField)[stationLoc]

pointPal <- pal[as.numeric(cut(anomLoc,breaks=256))]


binField     <- ifelse(!is.na(interpClip),1,NA)
noCoverfield <- ifelse(is.na(interpClip),1,NA)

# true field
pdf(sprintf('%s/trueField.pdf',plotdir),10,7)
image.plot(lon,lat,anomClip*landMask,col=pal, zlim=zr,
	xlab='', ylab='',
	main='ERA5 True Field')
world(add=T)

dev.off()

# approx station obs
pdf(sprintf('%s/stationField.pdf',plotdir),10,7)
image.plot(lon,lat,anomClip*z,col=pal, zlim=zr,
	xlab='', ylab='',
	main='1880s Station Distribution')
world(add=T)
points(pointLoc,pch=20,col=pointPal)
dev.off()

# interpolated field
pdf(sprintf('%s/interpField.pdf',plotdir),10,7)
image.plot(lon,lat,interpClip,col=pal, zlim=zr,
	xlab='', ylab='',
	main='1880s Interploated Field')
world(add=T)

# points(pointLoc,pch=20,col=pointPal)
dev.off()

# interpolated field
pdf(sprintf('%s/interpFieldGrey.pdf',plotdir),10,7)
image.plot(lon,lat,interpClip,col=pal, zlim=zr,
	xlab='', ylab='',
	main='1880s Interploated Field')
world(add=T)
image(lon,lat,noCoverfield*landMask,col=adjustcolor('black',alpha=0.2),add=T)
# points(pointLoc,pch=20,col=pointPal)
dev.off()


pdf(sprintf('%s/binaryField.pdf',plotdir),10,7)

image.plot(lon,lat,interpClip,col=pal, zlim=zr,
	xlab='', ylab='',
	main='1880s Interploated Coverage',
	legend=FALSE)
image()
world(add=T)
# plot the results of the global mean

dif <- resultsList$mask[[dec]]$global- resultsList$true$global
x <- seq(-1.5,1.5, length=300)
y <- dnorm(x,mean(dif),sd(dif))

pdf(sprintf('%s/difHist.pdf',plotdir),10,7)
hist(dif, breaks=25, freq=FALSE,xlab='', main='Mask Global Mean Error (1880s)',xlim=c(-1,1))
points(x,y,type='l',lwd=2,col='blue')

legend('topright', c(paste('Standard Dev:',round(sd(dif),3))),bty='n', text.col='blue')
dev.off()




################################################################################
# Series of plots demonstrating the interpolation process 1980s
################################################################################
dec <- 14
load(sprintf('Data/%s/InterpolatedFields/%s',reanalysis, files[dec]))


# smoosh the anomaly colorbar down
interpClip <- interpolatedField[,,eraTimePoint]
interpClip[interpClip>zMax] <- zMax
interpClip[interpClip< -zMax] <- -zMax


z <- decadalMasks[,,dec]
z[z==0] <- NA

gridList <- expand.grid(list(lon=lon,lat=lat))
stationLoc <- which(!is.na(z))

pointLoc <- gridList[stationLoc,]
anomLoc  <- c(trueField)[stationLoc]

pointPal <- pal[as.numeric(cut(anomLoc,breaks=256))]


# true field
pdf(sprintf('%s/trueField2.pdf',plotdir),10,7)
image.plot(lon,lat,anomClip*landMask,col=pal, zlim=zr,
	xlab='', ylab='',
	main='ERA5 True Field')
world(add=T)

dev.off()

# approx station obs
pdf(sprintf('%s/stationField2.pdf',plotdir),10,7)
image.plot(lon,lat,anomClip*z,col=pal, zlim=zr,
	xlab='', ylab='',
	main='1980s Station Distribution')
world(add=T)
points(pointLoc,pch=20,col=pointPal)
dev.off()

# interpolated field
pdf(sprintf('%s/interpField2.pdf',plotdir),10,7)
image.plot(lon,lat,interpClip,col=pal, zlim=zr,
	xlab='', ylab='',
	main='1980s Interpolated Field')
world(add=T)
# points(pointLoc,pch=20,col=pointPal)
dev.off()

dif <- resultsList$mask[[dec]]$global- resultsList$true$global
x <- seq(-1.5,1.5, length=300)
y <- dnorm(x,mean(dif),sd(dif))

pdf(sprintf('%s/difHist2.pdf',plotdir),10,7)
hist(dif, breaks=25, freq=FALSE,xlab='', main='Mask Global Mean Error (1980s)',xlim=c(-1,1))
points(x,y,type='l',lwd=2,col='blue')

legend('topright', c(paste('Standard Dev:',round(sd(dif),3))),bty='n', text.col='blue')
dev.off()

################################################################################
# bias visualization
################################################################################

pdf(sprintf('%s/lmFit.pdf',plotdir),14,7)
par(mfrow=c(1,2), mar=c(5, 5, 4, 3) + 0.1)


plot(resultsList$mask[[4]]$global, resultsList$true$global,
	xlab='Mask Global Mean', ylab='True Global Mean',
	main='Global Mean Bias (1880s)', 
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,col='black',
    ylim=c(-1, 1), xlim= c(-1,1))

grid(lwd=1.5)
abline(lm(resultsList$true$global ~ resultsList$mask[[4]]$global), col='black',lwd=2)
abline(0,1,col='darkgreen',lwd=1.5)


plot(resultsList$mask[[14]]$global, resultsList$true$global,
	xlab='Mask Global Mean', ylab='True Global Mean',
	main='Global Mean Bias (1980s)', 
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,col='black',
    ylim=c(-1, 1), xlim= c(-1,1))

grid(lwd=1.5)
abline(lm(resultsList$true$global ~ resultsList$mask[[14]]$global), col='black',lwd=2)
abline(0,1,col='darkgreen',lwd=1.5)


dev.off()


################################################################################
# Final processing of our results from the analysis
################################################################################

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
# SST visualization
################################################################################

inds <- which(tYear==1880):length(tYear)
yr <- range(annMatS[inds,])

pdf(sprintf('%s/sstDist.pdf',plotdir),10,7)
par(mar=c(5,5,4,2)+0.1)
plot(tYear[inds], annMatS[inds,1],ylim=yr,type='l', col=adjustcolor('black', alpha=0.2),
		xlab='Year', ylab='Ocean Temp. Anomaly (°C)', 
		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,
		main='GISTEMP Means of ERSSTv4 1000 Member Ensemble')
grid(lwd=1.5)

for(i in 2:1000){
	points(tYear[inds], annMatS[inds,i],type='l', col=adjustcolor('black', alpha=0.2))
}

dev.off()


pdf(sprintf('%s/sstDist2.pdf',plotdir),10,7)
par(mar=c(5,5,4,2)+0.1)
plot(tYear[inds], annMatS[inds,1],ylim=yr,type='l', col=adjustcolor('black', alpha=0.2),
		xlab='Year', ylab='Ocean Temp. Anomaly (°C)', 
		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,
		main='GISTEMP Means of ERSSTv4 1000 Member Ensemble')
grid(lwd=1.5)
for(i in 2:1000){
	points(tYear[inds], annMatS[inds,i],type='l', col=adjustcolor('black', alpha=0.2))
}

points(tYear[inds], apply(annMatS[inds,],1,mean),type='l',col='red', lwd=2)
legend('bottomright', 'Ensemble Mean', col='red', lwd=2, bty='n',cex=1.5)
dev.off()
################################################################################
# Figure 2
################################################################################
load(file='Data/Shared/GHCN/ghcnPlottingObjects.Rda')

inds <- 1:17
yr <- c(0,1)

pdf(sprintf('%s/02merraStationCoverage.pdf',plotdir),height=7,width=10)
par(mar=c(5,5,4,2)+0.1)
plot(tDec[inds],coverage3[inds],type='l',col='black',lwd=2,
	ylim=yr,xlab='Year',ylab='Land Coverage Proportion with 1200km Smoothing', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)

grid(lwd=gridlwd)

points(tDec[inds],coverage4[inds],type='l',col='black',lwd=2,lty=2)

points(tDec[inds],nhCoverage3[inds],type='l',col='blue',lwd=2)
points(tDec[inds],nhCoverage4[inds],type='l',col='blue',lwd=2,lty=2)

points(tDec[inds],shCoverage3[inds],type='l',col='red',lwd=2)
points(tDec[inds],shCoverage4[inds],type='l',col='red',lwd=2,lty=2)

abline(v=1880,lty=3,lwd=2)

legend('bottomright', 
	c('Global GHCNv3','Global GHCNv4', 'NH GHCNv3', 'NH GHCNv4', 'SH GHCNv3', 'SH GHCNv4'),
	col=c('black','black', 'blue', 'blue', 'red', 'red'),
		lwd=2,
		lty=c(1,3,1,3,1,3),cex=cexScale)

dev.off()

################################################################################
# Figure 3a: Comparison of LST error sources (all reanalyses)
################################################################################
inds <- which(tYearFull %in% tYearFinal)
load(sprintf('Data/%s/meanTimeSeries_%s.Rda','JRA','JRA'))
globalLandCIJRA <- totalLandUncertainty(resultsList,tDec)
sigmaLJRA     <- decadeToYear(globalLandCIJRA$diffVar,tYearFull,tDec)
totalLJRA     <- sigmaLJRA[inds] + homogVar

sciJRA <- sqrt(sigmaLJRA)* 1.96
tciJRA <- sqrt(totalLJRA) * 1.96

load(sprintf('Data/%s/meanTimeSeries_%s.Rda','MERRA','MERRA'))
globalLandCIMERRA <- totalLandUncertainty(resultsList,tDec)
sigmaLMERRA    <- decadeToYear(globalLandCIMERRA$diffVar,tYearFull,tDec)
totalLMERRA     <- sigmaLMERRA[inds] + homogVar

sciMERRA <- sqrt(sigmaLMERRA)* 1.96
tciMERRA <- sqrt(totalLMERRA) * 1.96

load(sprintf('Data/%s/meanTimeSeries_%s.Rda','ERA','ERA'))
globalLandCIERA <- totalLandUncertainty(resultsList,tDec)
sigmaLERA     <- decadeToYear(globalLandCIERA$diffVar,tYearFull,tDec)
totalLERA     <- sigmaLERA[inds] + homogVar

sciERA <- sqrt(sigmaLERA)* 1.96
tciERA <- sqrt(totalLERA) * 1.96


homogCI <- sqrt(homogVar) * 1.96
gis <- c(rep(0.2,20),NA,rep(0.15,50),rep(NA,9),rep(0.08,57))

yr <- c(0,max(sciERA,sciMERRA,sciJRA))

pdf(sprintf('%s/03lstSourcesComp.pdf',plotdir),20,7)
par(mfrow=c(1,2), mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFull, sciJRA, ylim=yr,type='l', lwd=2,
	xlab='Year', ylab='95% LSAT Confidence Interval (°C)', main = 'Comparison of LSAT Uncertainty', 
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,col='red')

points(tYearFull, sciMERRA, col='black', type='l', lwd=2)
points(tYearFull, sciERA, col='blue', type='l', lwd=2)

points(tYearFinal, tciJRA, col='red', type='l', lwd=2,lty=3)
points(tYearFinal, tciMERRA, col='black', type='l', lwd=2,lty=3)
points(tYearFinal, tciERA, col='blue', type='l', lwd=2,lty=3)

grid(lwd=gridlwd)

abline(v=1880,lty=3,lwd=2)
legend('topright', '(a)', cex = cexScale, bty='n')
legend('bottomleft',c('JRA55 Total','MERRA2 Total','ERA5 Total', 'JRA55 Sampling','MERRA2 Sampling','ERA5 Sampling'),
	col=rep(c('red','black','blue'),2),lwd=2,bg='white',lty=c(3,3,3,1,1,1),
	cex=1.25)

plot(tYearFull, sciERA, ylim=yr,type='l', lwd=2,
	xlab='Year', ylab='95% LSAT Confidence Interval (°C)', main = 'LSAT Uncertainty with ERA5',
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,col='blue')
grid(lwd=gridlwd)

points(tYearFinal,tciERA,type='l',col='black',lwd=2,lty=1)
points(tYearFinal,homogCI,type='l',col='red',lwd=2,lty=1)
points(tYearFinal,gis,type='l',col='blue',lwd=2,lty=2)

abline(v=1880,lty=3,lwd=2)
legend('topright', '(b)', cex = cexScale, bty='n')
legend('bottomleft', 
	c('Total LSAT Uncertainty', 'Homgenization Uncertainty', 'Sampling Uncertainty', 'Hansen et al. 2010'),
	col=c('black', 'red', 'blue','blue'), lwd=2,lty=c(1,1,1,2), cex=1.25,bg='white')

dev.off()


# clean plot of the sampling uncertainty
yr <- c(0,max(sciERA,sciMERRA,sciJRA))

pdf(sprintf('%s/samplingUncertainty.pdf',plotdir),10,7)
par(mfrow=c(1,1), mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFull, sciJRA, ylim=yr,type='l', lwd=2,
	xlab='Year', ylab='95% LSAT Confidence Interval (°C)', main = 'Comparison of LSAT Uncertainty', 
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,col='red')

points(tYearFull, sciMERRA, col='black', type='l', lwd=2)
points(tYearFull, sciERA, col='blue', type='l', lwd=2)


grid(lwd=gridlwd)

abline(v=1880,lty=3,lwd=2)
legend('topright', '(a)', cex = cexScale, bty='n')
legend('topright',c('JRA55 Sampling','MERRA2 Sampling','ERA5 Sampling'),
	col=rep(c('red','black','blue'),1),lwd=2,bg='white',lty=c(1,1,1),
	cex=1.25)

dev.off()

# plot of the ERA5 breakdown
pdf(sprintf('%s/eraUncertainty.pdf',plotdir),10,7)
plot(tYearFull, sciERA, ylim=yr,type='l', lwd=2,
	xlab='Year', ylab='95% LSAT Confidence Interval (°C)', main = 'LSAT Uncertainty with ERA5',
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,col='blue')
grid(lwd=gridlwd)

points(tYearFinal,tciERA,type='l',col='black',lwd=2,lty=1)
points(tYearFinal,homogCI,type='l',col='red',lwd=2,lty=1)
points(tYearFinal,gis,type='l',col='blue',lwd=2,lty=2)

abline(v=1880,lty=3,lwd=2)
legend('topright', 
	c('Total LSAT Uncertainty', 'Homgenization Uncertainty', 'Sampling Uncertainty', 'Hansen et al. 2010'),
	col=c('black', 'red', 'blue','blue'), lwd=2,lty=c(1,1,1,2), cex=1.25,bg='white')
dev.off()

# clean plot of the total uncertainty
pdf(sprintf('%s/totalUncertainty.pdf',plotdir),10,7)
par(mfrow=c(1,1), mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFull, sciJRA, ylim=yr,type='l', lwd=2,
	xlab='Year', ylab='95% LSAT Confidence Interval (°C)', main = 'Comparison of LSAT Uncertainty', 
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,col='red')

points(tYearFull, sciMERRA, col='black', type='l', lwd=2)
points(tYearFull, sciERA, col='blue', type='l', lwd=2)

points(tYearFinal, tciJRA, col='red', type='l', lwd=2,lty=3)
points(tYearFinal, tciMERRA, col='black', type='l', lwd=2,lty=3)
points(tYearFinal, tciERA, col='blue', type='l', lwd=2,lty=3)

grid(lwd=gridlwd)

abline(v=1880,lty=3,lwd=2)
legend('topright', '(a)', cex = cexScale, bty='n')
legend('topright',c('JRA55 Total','MERRA2 Total','ERA5 Total', 'JRA55 Sampling','MERRA2 Sampling','ERA5 Sampling'),
	col=rep(c('red','black','blue'),2),lwd=2,bg='white',lty=c(3,3,3,1,1,1),
	cex=1.25)

dev.off()


################################################################################
# Figure 4: Comparison of averaging methods
################################################################################
load(sprintf('Data/%s/simpleMean_%s.Rda',reanalysis,reanalysis))
load(sprintf('Data/%s/limitingGlobalCI_%s.Rda',reanalysis,reanalysis))

globalLandCISimple <- totalLandUncertainty(resultsListSimple,tDec)
x <- tYearFull
y1 <- decadeToYear(sqrt(globalLandCI$diffVar)*1.96,tYearFull,tDec)
y2 <- decadeToYear(sqrt(globalLandCISimple$diffVar)*1.96,tYearFull,tDec)

pdf(sprintf('%s/04meanComp.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)

plot(x,y1,ylim=c(0,max(y1,y2)),
	type='l',col='black',lwd=2, 
	xlab='Year',ylab='95% LSAT Sampling Confidence Interval',
	cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=1.5)
points(x,y2,type='l',col='red',lwd=2)

abline(v=1880,lty=3,lwd=2)

abline(h=globalLimitCI,col='blue',lwd=2,lty=2)
points(x,y1,type='l',col='black',lwd=2)

legend('topright', c('Simple Mean', 'GISTEMP Mean','Limiting Uncertainty'), 
	col=c('red', 'black','blue'),lty=c(1,1,2),lwd=2,bg='white',cex=cexScale)

dev.off()


################################################################################
# Figure 5: Hemispheric LST Uncertainties
################################################################################
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

# homog errors
nhHomog <- read.csv('Data/Shared/GHCN/NH-total-homogenization-error.csv',
				header = FALSE)
shHomog <- read.csv('Data/Shared/GHCN/SH-total-homogenization-error.csv',
				header = FALSE)
#assume symmetric errors
nhCI    <- (nhHomog[,3] - nhHomog[,2])/2
shCI    <- (shHomog[,3] - shHomog[,2])/2

# get the total land hemi var

landNHVar <- (sdnL)^2 + (c(rep(NA,31),nhCI)/1.96)^2
landSHVar <- (sdsL)^2 + (c(rep(NA,31),shCI)/1.96)^2

# Output plot (Maybe add weights?)


pdf(sprintf('%s/05hemiLSTUncertainty.pdf',plotdir),10,7)

par(mar=c(5, 5, 4, 3) + 0.1)
plot(tYearFull,sdnL*1.96,type='l',lwd=2,col = 'red',
	ylim = c(0,max(sdnL*1.96,sdsL*1.96,shCI,0.5)),
	xlab = 'Year',
	ylab = 'Hemispheric LSAT 95% Confidence Interval (°C)', 
    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
points(tYearFull,sdsL*1.96,type='l',lwd=2,col = 'blue')
grid(lwd=gridlwd)

legend('topright', 
	c('NH Sampling','SH Sampling','NH Homogenization','SH Homogenization'),
	lwd=2,lty=c(1,1,3,3),col=c('red','blue','red','blue'),cex=cexScale,bg='white')

points(tYearFinal,nhCI,type='l',col='red',lty=3,lwd=2)
points(tYearFinal,shCI,type='l',col='blue',lty=3,lwd=2)
abline(v=1880,lty=3,lwd=2)

dev.off()

################################################################################
# Figure 6: Band LST Uncertainties
################################################################################
tempResultsList <- resultsList

trueBand <- tempResultsList$true$bands

nDec <- length(tempResultsList$mask)

betaBand <- matrix(NA,nrow=nDec,ncol=8)
sigmaBand <- matrix(NA,nrow=nDec,ncol=8)
sdBand    <- matrix(NA, nrow=nDec,ncol=8)

for(i in 1:nDec){
	maskBand <- tempResultsList$mask[[i]]$bands


	for(j in 1:ncol(trueBand)){
		trueYear <- monthToYearMeans(trueBand[,j])
		maskYear <- monthToYearMeans(maskBand[,j])
		
		if(!all(is.nan(maskYear))){
			tempLM <- lm(trueYear ~ maskYear)
			betaBand[i,j] <- tempLM$coefficients[2]
			sigmaBand[i,j] <- summary(tempLM)$sigma
			sdBand[i,j] <- sd(trueYear - maskYear)
		}
	}
}

tYearLong <- 1850:2016

x <- tYearLong
y <- sdBand * 1.96

pdf(sprintf('%s/06bandLSTUncertainty.pdf',plotdir),10,7)

par(mfrow=c(4,2),mar=c(2,5,1,3)+0.1)

yr <- c(0,max(y,na.rm=T))

# right now band 1 is antartica running North to band 8 arctic. we want to
# rearrange them as R plots accross the rows of mfrow

order <- c(5,4,6,3,7,2,8,1)
ylabs <- c(rep('Tropics 95% CI',2),
		   rep('Sub-Tropics 95% CI',2),
		   rep('Mid-Latitude 95% CI',2),
		   rep('Polar 95% CI',2))

subDec <- 1:17
for(i in 1:ncol(sigmaBand)){
	if(i ==1 ){
		yTemp <- decadeToYear(y[subDec,order[i]],tYearLong,tDec[subDec])
		plot(x,yTemp,type='l', lwd=1.5,
			ylim=yr, main='Northern Hemishpere',
			xlab='Decade',ylab=ylabs[i], 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale+0.2)
		grid(lwd=gridlwd)
		points(x,yTemp,type='l', lwd=1.5)
		text(1855,max(yr)-0.05,paste('(', letters[i],')'), bty = 'n',cex=cexScale)
		abline(v=1880,col='black',lty=3,lwd=1.5)

	} else if (i ==2){
		yTemp <- decadeToYear(y[subDec,order[i]],tYearLong,tDec[subDec])
		plot(x,yTemp,type='l',ylim=yr, main='Southern Hemishpere',
			xlab='Decade',ylab=ylabs[i], 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale+0.2)
		grid(lwd=gridlwd)
		points(x,yTemp,type='l', lwd=1.5)
		text(1855,max(yr)-0.05,paste('(', letters[i],')'), bty = 'n',cex=cexScale)
		abline(v=1880,col='black',lty=3,lwd=1.5)
	} else{
		yTemp <- decadeToYear(y[subDec,order[i]],tYearLong,tDec[subDec])
		plot(x,yTemp,type='l',ylim=yr,
			xlab='Decade',ylab=ylabs[i], 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale+0.2)
		grid(lwd=gridlwd)
		points(x,yTemp,type='l', lwd=1.5)
		text(1855,max(yr)-0.05,paste('(', letters[i],')'), bty = 'n',cex=cexScale)
		abline(v=1880,col='black',lty=3,lwd=1.5)
	}
}


dev.off()

subDec <- 4:17
################################################################################
# Figure 7: Land Mean Series with CI (1880-2016)
################################################################################
totalLandCI    <- sqrt(totalLandVar) * 1.96
landSamplingCI <- sqrt(sigma2LFinal) * 1.96
homogCI        <- sqrt(homogVar) * 1.96

pdf(sprintf('%s/07lstCI.pdf',plotdir),10,14)
par(mfrow=c(2,1),mar=c(5, 5, 4, 3) + 0.1)
x <- tYearFinal
y <- gisLandOcean[,2]

ciEnvelope1 <- landSamplingCI
ciEnvelope2 <- totalLandCI - ciEnvelope1

yr <- c(min(y-ciEnvelope1),max(y+ciEnvelope2))
plot(x,y,lwd=1,type='l',col='black',ylim=yr,
		main='Global Annual Land Surface Temperature Anomaly',
	xlab='Year', ylab='Land Temp. Anomaly (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,y,lwd=1,type='l',col='black')

polygon(x=c(x,rev(x)), y=c((y-totalLandCI), (rev(y) + rev(totalLandCI))) ,col='red',border=NA)
polygon(x=c(x,rev(x)), y=c((y-ciEnvelope1), (rev(y) + rev(ciEnvelope1))) ,col='blue',border=NA)
points(x,y,lwd=2,type='l',col='black')

legend('topleft', c('Global Mean LSAT','Sampling Uncertainty', 'Sampling + Homogenization Uncertainty'),
	col=c('black','blue','red'),lwd=c(2,6,6),cex=cexScale,bg='white')

legend('bottomright', '(a)', cex = cexScale, bty='n')

# plot the smoothed time series with the same CIs
x <- tYearFinal
y <- gisLandOcean[,3]

ciEnvelope1 <- landSamplingCI
ciEnvelope2 <- totalLandCI - ciEnvelope1

plot(x,y,lwd=1,type='l',col='black',ylim=yr,
	main='Global Annual Land Surface Temperature Anomaly',
	xlab='Year', ylab='Land Temp. Anomaly (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,y,lwd=1,type='l',col='black')

polygon(x=c(x,rev(x)), y=c((y-totalLandCI), (rev(y) + rev(totalLandCI))) ,col='red',border=NA)
polygon(x=c(x,rev(x)), y=c((y-ciEnvelope1), (rev(y) + rev(ciEnvelope1))) ,col='blue',border=NA)
points(x,y,lwd=2,type='l',col='black')
legend('bottomright', '(b)', cex = cexScale, bty='n')


dev.off()

################################################################################
# Figure 8: Slope plots with SE (1850 - Present)
################################################################################
x <- tDec[subDec]+5

y <- globalLandCI$beta[subDec]
ci <- globalLandCI$se[subDec]
upper <- y + ci*1.96
lower <- y - ci*1.96

yr <- range(0.7, 1.05)

pdf(sprintf('%s/08lstSlopes.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)
plot(x,1/y,pch=19,ylim=yr,xlab='Decadal Distribution of Stations', ylab='Estimated Bias (95% CI)',
	main='Estimated Bias for 1980-2016 Land-Only Global Mean', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)

abline(h=1,lty=2, lwd=2, col='blue')
abline(v=1980, lwd=2, lty=3, col='red')
for(i in 1:length(x)) lines(x=c(x[i],x[i]),y=c(1/lower[i],1/upper[i]), lwd=2)
dev.off()


################################################################################
# Figure 8b: Slope plots with SE (All Reanalyses)
################################################################################
x1 <- tDec[subDec]+3.5
y1 <- globalLandCIJRA$beta[subDec]
ci1 <- globalLandCIJRA$se[subDec]
upper1 <- y1 + ci1*1.96
lower1 <- y1 - ci1*1.96

x2 <- tDec[subDec]+5
y2 <- globalLandCIMERRA$beta[subDec]
ci2 <- globalLandCIMERRA$se[subDec]
upper2 <- y2 + ci2*1.96
lower2 <- y2 - ci2*1.96

x3 <- tDec[subDec]+6.5
y3 <- globalLandCIERA$beta[subDec]
ci3 <- globalLandCIERA$se[subDec]
upper3 <- y3 + ci3*1.96
lower3 <- y3 - ci3*1.96

yr <- range(1/upper1,1/upper2,1/upper3,1/lower1,1/lower2,1/lower3)

pdf(sprintf('%s/08lstSlopesComp.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)
plot(x1,1/y1,pch=19,ylim=yr,xlab='Decadal Distribution of Stations', ylab='Estimated Bias (95% CI)',
	main='Estimated Bias for Land-Only Global Mean', col='red',
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)

points(x2,1/y2,pch=19,col='black')
points(x3,1/y3,pch=19,col='blue')

abline(h=1,lty=2, lwd=2, col='black')
abline(v=1980, lwd=2, lty=3, col='red')

for(i in 1:length(x)) lines(x=c(x1[i],x1[i]),y=c(1/lower1[i],1/upper1[i]), lwd=2, col='red')
for(i in 1:length(x)) lines(x=c(x2[i],x2[i]),y=c(1/lower2[i],1/upper2[i]), lwd=2, col='black')
for(i in 1:length(x)) lines(x=c(x3[i],x3[i]),y=c(1/lower3[i],1/upper3[i]), lwd=2, col='blue')

legend('bottomright',c('JRA55','MERRA2','ERA5'),col=c('red','black','blue'),lwd=2,pch=19,bg='white',
	cex=cexScale)
dev.off()


################################################################################
# Figure 8b CORRECTED: Slope plots with SE (All Reanalyses) CORRECTED
################################################################################
x1 <- tDec[subDec]+3.5
y1 <- globalLandCIJRA$beta2[subDec]
ci1 <- globalLandCIJRA$se2[subDec]
upper1 <- y1 + ci1*1.96
lower1 <- y1 - ci1*1.96

x2 <- tDec[subDec]+5
y2 <- globalLandCIMERRA$beta2[subDec]
ci2 <- globalLandCIMERRA$se2[subDec]
upper2 <- y2 + ci2*1.96
lower2 <- y2 - ci2*1.96

x3 <- tDec[subDec]+6.5
y3 <- globalLandCIERA$beta2[subDec]
ci3 <- globalLandCIERA$se2[subDec]
upper3 <- y3 + ci3*1.96
lower3 <- y3 - ci3*1.96

yr <- range(upper1,upper2,upper3,lower1,lower2,lower3)

pdf(sprintf('%s/08lstSlopesCompCORRECTED.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)
plot(x1,y1,pch=19,ylim=yr,xlab='Decadal Distribution of Stations', ylab='Estimated Bias (95% CI)',
	main='Estimated Bias for Land-Only Global Mean (Corrected)', col='red',
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)

points(x2,y2,pch=19,col='black')
points(x3,y3,pch=19,col='blue')

abline(h=1,lty=2, lwd=2, col='black')
abline(v=1980, lwd=2, lty=3, col='red')

for(i in 1:length(x)) lines(x=c(x1[i],x1[i]),y=c(1/lower1[i],1/upper1[i]), lwd=2, col='red')
for(i in 1:length(x)) lines(x=c(x2[i],x2[i]),y=c(1/lower2[i],1/upper2[i]), lwd=2, col='black')
for(i in 1:length(x)) lines(x=c(x3[i],x3[i]),y=c(1/lower3[i],1/upper3[i]), lwd=2, col='blue')

legend('bottomright',c('JRA55','MERRA2','ERA5'),col=c('red','black','blue'),lwd=2,pch=19,bg='white',
	cex=cexScale)
dev.off()


################################################################################
# Figure 9: Ocean Mean Series with CI from sigma2SG (1880-2014)
################################################################################
pdf(sprintf('%s/09sstCI.pdf',plotdir),10,14)
par(mfrow=c(2,1), mar=c(5, 5, 4, 3) + 0.1)

ciCol <- adjustcolor('black',alpha=0.3)


x <-tYearFinal
y <- gisLandOcean[,4]
ci <- 1.96 * sqrt(sigma2SFinal)
yr <- c(min(y-ci),max(y+ci))


plot(x,y,type='l',lwd=1,ylim=yr,
		xlab='Year', ylab='Ocean Temp. Anomaly (°C)', 
		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,
		main='Global Annual Ocean Surface Temperature Anomaly')
grid(lwd=gridlwd)
points(x,y,type='l',lwd=1)

polygon(x=c(x,rev(x)), y=c((y-ci), (rev(y) + rev(ci))) ,col=ciCol,border=NA)
legend('topleft', '(a)', cex = cexScale, bty='n')


x <-tYearFinal
y <- gisLandOcean[,5]


plot(x,y,type='l',lwd=1,ylim=yr,
		xlab='Year', ylab='Ocean Temp. Anomaly (°C)', 
		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,
		main='Global Annual Ocean Surface Temperature Anomaly')
grid(lwd=gridlwd)
points(x,y,type='l',lwd=1)

polygon(x=c(x,rev(x)), y=c((y-ci), (rev(y) + rev(ci))) ,col=ciCol,border=NA)

legend('bottomright', c('Sea Surface Temperature', '95% Confidence Interval'),
		lwd=c(1,6), col=c('black',ciCol) ,cex=cexScale,bg='white')
legend('topleft', '(b)', cex = cexScale, bty='n')

dev.off()


pdf(sprintf('%s/sstCI.pdf',plotdir),10,7)
par(mfrow=c(1,1), mar=c(5, 5, 4, 3) + 0.1)

ciCol <- adjustcolor('black',alpha=0.3)


x <-tYearFinal
y <- gisLandOcean[,4]
ci <- 1.96 * sqrt(sigma2SFinal)
yr <- c(min(y-ci),max(y+ci))


plot(x,y,type='l',lwd=1,ylim=yr,
		xlab='Year', ylab='Ocean Temp. Anomaly (°C)', 
		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale,
		main='Global Annual Ocean Surface Temperature Anomaly')
grid(lwd=gridlwd)
points(x,y,type='l',lwd=1)

polygon(x=c(x,rev(x)), y=c((y-ci), (rev(y) + rev(ci))) ,col=ciCol,border=NA)
dev.off()
################################################################################
# Figure 10: Hemispheric Ocean Uncertanties (1854-2014)
################################################################################
annMatSNH <- extratTSRdatNH(ensemblesS,'Data/Shared/ERSST',tYear)
sigma2SNH  <- apply(annMatSNH,1,var)

annMatSSH <- extratTSRdatSH(ensemblesS,'Data/Shared/ERSST',tYear)
sigma2SSH  <- apply(annMatSSH,1,var)



# deal with date stuff
sdn <- c(sqrt(sigma2SNH)[tYear %in% tYearFinal], 
			rep(sqrt(sigma2SNH[length(sigma2SNH)]),2))
sds <- c(sqrt(sigma2SSH)[tYear %in% tYearFinal], 
			rep(sqrt(sigma2SSH[length(sigma2SSH)]),2))


pdf(sprintf('%s/10sstHemiCI.pdf',plotdir),10,7)
par(mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFinal,sdn*1.96,type='l',lwd=2,col = 'red',
	ylim = c(0,1.96*max(sdn,sds)),
	xlab = 'Year',
	ylab = 'SST Hemispheric 95% Confidence Interval (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(tYearFinal,sdn*1.96,type='l',lwd=2,col = 'red')

points(tYearFinal,sds*1.96,type='l',lwd=2,col = 'blue')
legend('topright', c('NH','SH'),lwd=2,lty=1,col=c('red','blue'),
	cex=cexScale,bg='white')

dev.off()


################################################################################
# Figure 11: Final Global Annual time series (No correlation)
################################################################################

lstFinal <- sqrt(AL^2 * totalLandVar)*1.96
lstsstFinal <- sqrt(AL^2 * totalLandVar + AS^2 * sigma2SFinal)*1.96

ciEnvelope2 <- c(lstsstFinal,rep(lstsstFinal[length(lstsstFinal)],2))
ciEnvelope3 <- c(lstFinal,rep(lstFinal[length(lstFinal)],2))


pdf(sprintf('%s/11totalTimeSeriesNoCorr.pdf',plotdir),10,14)
par(mfrow=c(2,1),mar=c(5, 5, 4, 3) + 0.1)

x <- tYear2018
y <- gisGlobal2018[,2]

plot(x,y,lwd=1,type='l',col='black',ylim=c(min(y-ciEnvelope2),max(y+ciEnvelope2)),
	xlab='Year', ylab='Total Global Annual Mean Surface Temperature (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,y,lwd=1,type='l',col='black')

polygon(x=c(x,rev(x)), y=c((y-ciEnvelope2), (rev(y) + rev(ciEnvelope2))) ,col='blue',border=NA)
polygon(x=c(x,rev(x)), y=c((y-ciEnvelope3), (rev(y) + rev(ciEnvelope3))) ,col='green',border=NA)

points(x,y,lwd=2,type='l',col='black')
legend('topleft', '(a)', cex = cexScale, bty='n')

# plot the smoothed time series with the same CIs
goodInds <- which(!is.na(gisGlobal2018[,3]))

x <- tYear2018[goodInds]
y <- gisGlobal2018[goodInds,3]

plot(x,y,lwd=1,type='l',col='black',ylim=c(min(y-ciEnvelope2[goodInds],na.rm=T),max(y+ciEnvelope2[goodInds],na.rm=T)),xlim=range(tYear2018),
	xlab='Year', ylab='Total Global Annual Mean Surface Temperature (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,y,lwd=1,type='l',col='black')

polygon(x=c(x,rev(x)), y=c((y-ciEnvelope2[goodInds]), (rev(y) + rev(ciEnvelope2[goodInds]))) ,col='blue',border=NA)
polygon(x=c(x,rev(x)), y=c((y-ciEnvelope3[goodInds]), (rev(y) + rev(ciEnvelope3[goodInds]))) ,col='green',border=NA)

points(x,y,lwd=2,type='l',col='black')

legend('bottomright', c('Global Mean Temperature','LSAT Uncertainty', 'LSAT + SST Uncertainty'),
	col=c('black','green','blue'),lwd=c(2,6,6),
	cex=cexScale,horiz=F,bg='white')
legend('topleft', '(b)', cex = cexScale, bty='n')
dev.off()

save(gisGlobal2018,ciEnvelope2,file=c(sprintf('Data/%s/totalCI_%s.Rda',reanalysis,reanalysis)))



pdf(sprintf('%s/totalSeries.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)

x <- tYear2018
y <- gisGlobal2018[,2]

plot(x,y,lwd=1,type='l',col='black',ylim=c(min(y-ciEnvelope2),max(y+ciEnvelope2)),
	xlab='Year', ylab='Global Annual Anomaly (°C)', 
	main = 'Global Annual Mean Surface Temperature Anomaly',
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,y,lwd=1,type='l',col='black')

polygon(x=c(x,rev(x)), y=c((y-ciEnvelope2), (rev(y) + rev(ciEnvelope2))) ,col='blue',border=NA)
polygon(x=c(x,rev(x)), y=c((y-ciEnvelope3), (rev(y) + rev(ciEnvelope3))) ,col='green',border=NA)

points(x,y,lwd=2,type='l',col='black')

legend('topleft', c('Global Anomaly','LSAT Uncertainty', 'LSAT + SST Uncertainty'),
	col=c('black','green','blue'),lwd=c(2,6,6),
	cex=cexScale,horiz=F,bty='n')

dev.off()
################################################################################
# Figure 11b: Final Global Annual time series (Flipped Order)
################################################################################

sstFinal <- sqrt(AS^2 * sigma2SFinal)*1.96
lstsstFinal <- sqrt(AL^2 * totalLandVar + AS^2 * sigma2SFinal)*1.96

ciEnvelope2 <- c(lstsstFinal,rep(lstsstFinal[length(lstsstFinal)],2))
ciEnvelope3 <- c(sstFinal,rep(sstFinal[length(sstFinal)],2))


pdf(sprintf('%s/11totalTimeSeriesNoCorrFlip.pdf',plotdir),10,14)
par(mfrow=c(2,1),mar=c(5, 5, 4, 3) + 0.1)

x <- tYear2018
y <- gisGlobal2018[,2]

plot(x,y,lwd=1,type='l',col='black',ylim=c(min(y-ciEnvelope2),max(y+ciEnvelope2)),
	xlab='Year', ylab='Total Global Annual Mean Surface Temperature (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,y,lwd=1,type='l',col='black')

polygon(x=c(x,rev(x)), y=c((y-ciEnvelope2), (rev(y) + rev(ciEnvelope2))) ,col='green',border=NA)
polygon(x=c(x,rev(x)), y=c((y-ciEnvelope3), (rev(y) + rev(ciEnvelope3))) ,col='blue',border=NA)

points(x,y,lwd=2,type='l',col='black')
legend('topleft', '(a)', cex = cexScale, bty='n')

# plot the smoothed time series with the same CIs
goodInds <- which(!is.na(gisGlobal2018[,3]))

x <- tYear2018[goodInds]
y <- gisGlobal2018[goodInds,3]

plot(x,y,lwd=1,type='l',col='black',ylim=c(min(y-ciEnvelope2[goodInds],na.rm=T),max(y+ciEnvelope2[goodInds],na.rm=T)),xlim=range(tYear2018),
	xlab='Year', ylab='Total Global Annual Mean Surface Temperature (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,y,lwd=1,type='l',col='black')

polygon(x=c(x,rev(x)), y=c((y-ciEnvelope2[goodInds]), (rev(y) + rev(ciEnvelope2[goodInds]))) ,col='green',border=NA)
polygon(x=c(x,rev(x)), y=c((y-ciEnvelope3[goodInds]), (rev(y) + rev(ciEnvelope3[goodInds]))) ,col='blue',border=NA)

points(x,y,lwd=2,type='l',col='black')

legend('bottomright', c('Global Mean Temperature','SST Uncertainty', 'LSAT + SST Uncertainty'),
	col=c('black','blue','green'),lwd=c(2,6,6),
	cex=cexScale,horiz=F,bg='white')
legend('topleft', '(b)', cex = cexScale, bty='n')
dev.off()

save(gisGlobal2018,lstsstFinal,file=sprintf('Data/%s/totalCI_%s.Rda',reanalysis,reanalysis))

csvOut <- data.frame(year=tYear2018,gistemp=gisGlobal2018[,2],ci95=ciEnvelope2)

write.csv(csvOut,file=sprintf('Data/%s/totalCI_%s.csv',reanalysis,reanalysis),row.names=FALSE)


# figure for website

pdf(sprintf('%s/websiteSeries.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)

x <- tYear2018[goodInds]
y <- gisGlobal2018[goodInds,3]

plot(x,y,lwd=1,type='l',col='black',ylim=c(min(y-ciEnvelope2[goodInds],na.rm=T),max(y+ciEnvelope2[goodInds],na.rm=T)),xlim=range(tYear2018),
	xlab='Year', ylab='Total Global Annual Mean Surface Temperature (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,y,lwd=1,type='l',col='black')

polygon(x=c(x,rev(x)), y=c((y-ciEnvelope2[goodInds]), (rev(y) + rev(ciEnvelope2[goodInds]))) ,col='green',border=NA)
polygon(x=c(x,rev(x)), y=c((y-ciEnvelope3[goodInds]), (rev(y) + rev(ciEnvelope3[goodInds]))) ,col='blue',border=NA)

points(x,y,lwd=2,type='l',col='black')

legend('bottomright', c('Global Mean Temperature','SST Uncertainty', 'LSAT + SST Uncertainty'),
	col=c('black','blue','green'),lwd=c(2,6,6),
	cex=cexScale,horiz=F,bg='white')
dev.off()


################################################################################
# Figure 12: Final Hemi Annual uncertainty series 
################################################################################
x <- tYearFinal


# get the three standard deviations (2hemi 1global) in easy to play around with
# variables to figure out what they are describing
sdnLFinal <- decadeToYear(sigmaHemi[subDec,1],tYearFinal,tDec[subDec])
sdsLFinal <- decadeToYear(sigmaHemi[subDec,2],tYearFinal,tDec[subDec])

# get the land hemi var
landNHVar <- (sdnLFinal)^2 + (nhCI/1.96)^2
landSHVar <- (sdsLFinal)^2 + (shCI/1.96)^2

# get the total CIs
totalNHCI <- sqrt((1-ASn)^2 * landNHVar + ASn^2 * sdn^2) * 1.96
totalSHCI <- sqrt((1-ASs)^2 * landSHVar + ASs^2 * sds^2) * 1.96

pdf(sprintf('%s/12totalHemiCI.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)

plot(x,totalNHCI,type='l',lwd=2,col = 'red',
	ylim = c(0,max(totalNHCI,totalSHCI)),
	xlab = 'Year',
	ylab = 'Total Hemispheric 95% Confidence Interval (°C)', 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
points(x,totalNHCI,type='l',lwd=2,col='red')

points(x,totalSHCI,type='l',lwd=2,col = 'blue')
legend('topright', c('NH','SH'),lwd=2,lty=1,col=c('red','blue'),
	cex=cexScale,bg='white')

dev.off()


################################################################################
# Figure 13: Final CI comparison
################################################################################

giss <- read.csv("Data/ERA/totalCI_ERA.csv",header=FALSE,skip=1)
best <- read.table("Data/Shared/otherProducts/Best_uncertainty.txt",skip=48,header=FALSE)
had  <- read.table("Data/Shared/otherProducts/HadCRUT.4.6.0.0.annual_ns_avg.txt",header=FALSE)

pdf(sprintf('%s/13compCI.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)

plot(had[,1],(had[,12]-had[,11])/2,type="l",ylim=c(0,0.25),
	yaxs="i",lwd=2,xlim=c(1850,2020),xaxs="i",
	ylab="95% Total Uncertainty (ºC)",xlab="Year",main="Comparison of Uncertainty Estimates", 
    		cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)
grid(lwd=gridlwd)
lines(best[,1],best[,3],col=4,lwd=2)
lines(giss[,1],giss[,3],col=2,lwd=2)
legend(1950,0.23,col=c(1,2,4),lwd=c(3,3,3),box.lwd=NA,legend=c("HadCRUT4","GISTEMP","Berkeley Earth"),cex=1.5)


dev.off()


################################################################################
# Visualize the warmest year
################################################################################
resultsPresent <- resultsMC(gisGlobal2018[,2],lstsstFinal[],gisGlobal2018[,1],alpha=0.5,M=10000)

inds <- which(giss[,1] > 2009)

meanVec <- giss[inds,2]
sdVec   <- giss[inds,3]/1.96

sortIndex <- sort(meanVec,index.return=TRUE)$ix
pal <- brewer.pal(9, 'Paired')[sortIndex]

xMin <- 0.5
xMax <- 1.1
yMax <- 20

pdf(sprintf('%s/warmestYear.pdf',plotdir),10,7)
x <- seq(xMin, xMax, length=256)
y <- dnorm(x,meanVec[1], sdVec[1])
plot(x,y,ylim=c(0,yMax),type='l',col=pal[1], lwd=1.5,xlab='Mean Temperature Anomaly (°C)',ylab='',yaxt='n')
text(meanVec[1],max(y)+0.5, paste(giss[inds[1],1]),col=pal[1])

for(i in 2:length(inds)){
	y <- dnorm(x,meanVec[i], sdVec[i])
	points(x,y,ylim=c(0,yMax),type='l',col=pal[i], lwd=1.5)
	text(meanVec[i],max(y)+0.5, paste(giss[inds[i],1]),col=pal[i])

	if(giss[inds[i],1]==2016){
		text(meanVec[i], max(y)-8, '86%', col=pal[i])
		text(meanVec[i], max(y)-9, '(87%)', col=pal[i])
	}
}
dev.off()