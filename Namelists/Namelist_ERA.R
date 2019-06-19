# The 'namelist' for the ERA5 reanalysis. Modifications to the code should
# be made through the namelist as much as possible to ensure major calculations
# are kept consistent.

# GISTEMP Uncertainty Analysis
# Version 1.0 (May 1, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

# reanalysis to use ('ERA', 'MERRA', 'JRA')
reanalysis <- 'ERA'

# control the parallelization of the interpolation
nCores <- 3

# which steps to run 
steps <- c(1:5)

plotdir <- 'Figures'

# Packages
library(ncdf4)
library(fields)
library(doParallel)
library(foreach)

# user funtions
source('Code/Functions.R')

# start year (1850, 1852, or 1880 are most common)
startYear <- 1850

# End year for the reanalysis
endYear   <- 2018

# smoothing radius
radius <- 1200


# option to set things differently due to memory issues
interpolateCores <- nCores
meanCores <- nCores


# Preprocessing stuff (step 0)
ddirRaw <- '/Users/lenssen/Documents/gistempRaw'

infoFile <- 'inv.ghcn4'
dataFile <- 'gistemp.step2.ghcn4'

