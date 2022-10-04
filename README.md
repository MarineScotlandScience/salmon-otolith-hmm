#   Inferring marine migration from otolith biogeochemical signatures of a wide-ranging fish

## Requirements

### Otolith data

TODO: link to data portal

Assumed path to data:

    data/ototlith_sst.csv

### OISST data

NOAA High Resolution SST data provided by the NOAA/OAR/ESRL PSL, Boulder, Colorado, USA, 
from their Web site at 
https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html

Need the files

    sst.day.mean.2007.nc
    sst.day.mean.2008.nc
    sst.day.mean.2009.nc
    sst.day.mean.2010.nc
    lsmask.oisst.v2.nc

For the model need to combine all sst data into one file with NCO software
available here: http://nco.sourceforge.net/#Binaries

Ubuntu / debian use command

    ncrcat sst.day.mean.20??.nc sst.day.mean.nc

Also need NetCDF installation for RNetCDF library

Assumed paths to data:

    data/sst.day.mean.nc
    data/lsmask.oisst.v2.nc

### R LIBRARIES

    Rcpp
    RNetCDF 
    progress
    tidyverse
    lubridate


## Usage

Run scripts in order. 

+ Run `00_simulate_migration.R` with both the `north` and `west` configuration.
+ `01_wrangle_data.R` expects simulation files from these simulations to
piece together with the real data. 
+ `02_hmm.R` Runs the HMM on the specified time series
+ `03_summarise.R` Generate summary plots etc. on specified HMM output.
