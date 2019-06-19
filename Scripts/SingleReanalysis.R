# A script to run the entire uncertainty analysis for one reanalysis product.
# Proper use is setting the steps to run using a namelist and then sourcing
# the script. In general, step 0 should not be run unless the analysis is
# being performed from the raw data as it is time intensive.

# GISTEMP Uncertainty Analysis
# Version 1.0 (May 1, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)

###############################################################################
# (STEP 0) These steps are generally skipped as they rely on the full
# datasets which are very large and will not be distributed
###############################################################################
if(0 %in% steps){
	# process all of the raw reanalysis data
	source('Code/Preprocessing/ProcessMERRA.R')
	source('Code/Preprocessing/ProcessERA5.R')
	source('Code/Preprocessing/ProcessJRA.R')

	# process the ERSST ensemble (Long to Run!)
	source('Code/Preprocessing/ProcessERSST.R')				

	# create the seaIceMask for the datasets (must run MERRA first)
	source('Code/Preprocessing/SeaIceMasks.R') 

	# process the station data from the GISTEMP output
	source('Code/Preprocessing/ProcessStationData.R')
}

###############################################################################
# (STEP 1) Create the coverage masks on the on the reanalysis grid and 
# calculate areas
###############################################################################
if(1 %in% steps){
	print('Step 1')

	source('Code/CreateMask.R')
	source('Code/AreaCalculation.R')
}

###############################################################################
# (STEP 2) Make the interpolated field
###############################################################################
if(2 %in% steps){
	print('Step 2')
	source('Code/InterpolateField.R')
}

###############################################################################
# (STEP 3) Take the means
###############################################################################
if(3 %in% steps){
	print('Step 3')
	source('Code/MeanCalculation.R')
}

###############################################################################
# (STEP 4) Run extended analyses (not reqired for operational 
# uncertainty estimates)
###############################################################################
if(4 %in% steps){
	print('Step 4')
	source('Code/LandInterpolationLimit.R')
	source('Code/SimpleMeanExperiment.R')
}

###############################################################################
# (STEP 5) Process Results and draw plots
###############################################################################
if(5 %in% steps){
	print('Step 5')
	source('Code/PaperPlots.R')
	source('Code/TempChangeAndWarmestYear.R')
}

