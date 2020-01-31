The code and data repository for the GISTEMP uncertainty analysis. More details can be found in Lenssen et al. (2019)

# Quick Start Guide
0. Make sure the required packages are installed for your R configuration by running
```install.packages(c('fields', 'doParallel', 'foreach', 'ncdf4'))```
1. Check the Namelists to make sure that the number of cores used in the calculation (variable nCores) is compatible for your machine. We suggest starting with 2 cores for 8 GB of RAM and 3 cores for 16 GB.
2. Run the full analysis with `source('Scripts/FullAnalysis.R')` Note: takes around 8 hours with 3 cores on 2018 Macbook Pro

# Download Options
* `gistempUncertaintyFull.tar.gz` (14.3 GB) Contains all raw, processed, intermediate, and final data and results.
* `gistempUncertaintyRaw.tar.gz` (1.3 GB) Contains all raw and processed data. Runs the analysis without any additional downloaded data.
* `gistempUncertaintyCode.tar.gz` (96 kB) Contains only the code. Raw data needs to be downloaded from the respective repositories.

We recommend downloading `gistempUncertaintyRaw.tar` (1.3 GB) to avoid downloading and processing the raw reanalysis and ERSST data which is not provided as part of our repository. 

# Code Organization

The codebase is to be run with `GistempUncertainty/` as the working directory. It should have the following subdirectories

* `Code/` The analysis scrips and user functions for the analysis.
..* `Preprocessing` Scripts used in the preprocessing of raw data products.
* `Data/` All raw and intermediate data used in the analysis.
* `Figures/` Plots saved from `Code/paperPlots.R` .
* `Namelists/` Scripts that set run parameters. These should be altered before diving deeper into the code.
* `Scripts/` Full analysis scripts.


# Version Changes
**1.2.1** Changes to ERSST preprocessing for correction to SST ensemble mean calculation.

**1.2.0** Changes correcting an error in the mean calculation. Results of sampling uncertainty for Global and SH means affected

**1.1.0** Changes correcting an error in the bias uncertainty calculation

**1.0.0** Original code release
