# Constants

# Body lengths per second
BL_PS <- 1.25

# Standard deviation in likelihood calculation
SD <- 1

# Time steps
DAYS_PER_TS = 7
SECONDS_TO_TS = DAYS_PER_TS*24*60*60
SWIM_SPEED_PBL = BL_PS*SECONDS_TO_TS

# Units 
# DISTANCE = KM
MAX_LAT = 80
DEG_NS_TO_KM = 110.574
DEG_EW_TO_KM = 111.320 # This should be multiplied by cos(lat)
R_km <- 6371

# Region for simulation (and plot)
# Hi res resolution is 0.25  degrees. 
# With values at .125 .375 .625 .875
SIM_RES = 0.25 # Degrees (width of grid element) 
SIM_LON = seq(-75.125, 25.125, by = SIM_RES)
SIM_LAT = seq(45.125, MAX_LAT + .125, by = SIM_RES)

# Functions

# Adjust lat / lot to a grid point in the model
snap_to_grid = function(x){ return(x - x%%SIM_RES + SIM_RES/2) }

# Functions to transform our coordinate system to lon and lat indices of sst 
# data
# For the high res sims - 0.25 degree resolution with .125 offset
# Lon indices go from 0 - 360 so need modulo 360 to deal with negative coords,
# - 0.125 to go to index value, plus 1 to get R indices (i.e 0 degrees = index 1)
# times 4 to retrieve index
lon_index = function(x){ ((x - .125) %% 360 ) * 4 + 1 }
# Lat indices go from 0 - 180, with 0 at south pole so need 
# to add 90 to get corresponding indices
# - 0.125 to go to index value, *4 and add 1
lat_index = function(x){ ((x - 0.125) + 90) * 4 + 1 }