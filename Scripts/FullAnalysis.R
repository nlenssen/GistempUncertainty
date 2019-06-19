# A script to run the entire uncertainty analysis for all three reanalyses
# resulting in a full replication of the Lenssen et al paper.

# GISTEMP Uncertainty Analysis
# Version 1.0 (May 1, 2019)
# Nathan Lenssen (lenssen@ldeo.columbia.edu)


# run the MERRA analysis first as it generates some shared results
source('Namelists/Namelist_MERRA.R')
source('Scripts/SingleReanalysis.R')

rm(list=ls())
gc()

# then JRA
source('Namelists/Namelist_MERRA.R')
source('Scripts/SingleReanalysis.R')

rm(list=ls())
gc()

# Finally, ERA (which will draw figures)
source('Namelists/Namelist_MERRA.R')
source('Scripts/SingleReanalysis.R')

rm(list=ls())
gc()