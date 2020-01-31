plotdir <- '~/Documents/gistempWebsite/websiteErrataMaterials'

# all of the time indicies used in plotting
tYear <- 1854:2014
tDec  <- seq(1850, 2010, by=10)
tYearFull <- 1850:2017
tYearFinal <- 1880:2016
tYear2018 <- 1880:2018

# decades to be used when plotting a limited range
subDec <- 4:17

# Old SST
annMatSWrong  <- extratTSRdat(ensemblesS,'~/Documents/gistempRaw/ERSSTMeanWrong',tYear)
sigma2SWrong  <- apply(annMatSWrong,1,var)
sigma2SFinalWrong <- c(sigma2SWrong[tYear %in% tYearFinal], rep(sigma2SWrong[length(sigma2SWrong)],2))

# New SST
annMatS  <- extratTSRdat(ensemblesS,'Data/Shared/ERSST',tYear)
sigma2S  <- apply(annMatS,1,var)
sigma2SFinal <- c(sigma2S[tYear %in% tYearFinal], rep(sigma2S[length(sigma2S)],2))

# Show changes to the SST uncertainty
pdf(sprintf('%s/sstCorrection.pdf',plotdir),10,7)

par(mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFinal,sqrt(sigma2SFinalWrong)*1.96,type='l',lwd=1.5,col='red',ylim=c(0,0.2),xlab='Year',
	ylab='Annual Mean SST Uncertainty 95% Confidence Interval', main='Corrections to the GISTEMP Global SST Uncertainty')
points(tYearFinal,sqrt(sigma2SFinal)*1.96,type='l',lwd=1.5)
grid(lwd=1.5)
legend('bottomleft',c('Original','Correction'),
	col=rep(c('red','black'),1),lwd=2,bg='white',cex=1.25)
dev.off()





# Show changes to the resulting total uncertainty

# load things
reanalysis <- 'ERA'
load(sprintf('Data/%s/meanTimeSeries_%s.Rda',reanalysis,reanalysis))
load(sprintf('Data/%s/area_%s.Rda',reanalysis,reanalysis))

totalHomog <- read.csv('Data/Shared/GHCN/total-homog-uncertainty.csv',header=FALSE)

paramHomog <- read.csv('Data/Shared/GHCN/parametric-uncertainty.csv',header=FALSE)

# calculate the LSAT
globalLandCI <- totalLandUncertainty(resultsList,tDec)

# Global Uncertainties:

# land calcs (empirical estimates) [take only through 2014]
sigma2L   <- decadeToYear(globalLandCI$diffVar,tYear,tDec)[1:length(tYear)]

# land variance calcs running over all 1880-2016
sigma2LFinal   <- decadeToYear(globalLandCI$diffVar,tYearFull,tDec)[tYearFull %in% tYearFinal]

# total land CI calculation (adding in homogonization uncertainty)
totalCIWrong <- sqrt(AL^2 * totalLandVar + AS^2 * sigma2SFinalWrong)*1.96
totalCI      <- sqrt(AL^2 * totalLandVar + AS^2 * sigma2SFinal)*1.96


# make the plot!
pdf(sprintf('%s/totalCorrection.pdf',plotdir),10,7)

par(mar=c(5, 5, 4, 3) + 0.1)

plot(tYearFinal,totalCIWrong,type='l',lwd=1.5,col='red',ylim=c(0,0.2),xlab='Year',
	ylab='Annual Mean SST Uncertainty 95% Confidence Interval', main='Corrections to the GISTEMP Global SST Uncertainty')
points(tYearFinal,totalCI,type='l',lwd=1.5)
grid(lwd=1.5)
legend('bottomleft',c('Original','Correction'),
	col=rep(c('red','black'),1),lwd=2,bg='white',cex=1.25)
dev.off()
