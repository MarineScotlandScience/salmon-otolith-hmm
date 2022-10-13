# HMM
#
# Code to run the hmm comparing otolith derived (or simulated) sst values
# with observed sst values.
#
# Largely follows:
#
# Thygesen, U. H., Pedersen, M. W., & Madsen, H. (2009). 
# Geolocating Fish Using Hidden Markov Models and Data Storage Tags. 
# In J. L. Nielsen, H. Arrizabalaga, N. Fragoso, A. Hobday, M. Lutcavage, & J. Sibert (Eds.), 
# Tagging and tracking of marine animals with electronic devices (pp. 277–293). 
# Springer Netherlands. https://doi.org/10.1007/978-1-4020-9640-2_17
#
# Requires:
#   data/hmm_data.rds
#   data/lsmask.oisst.v2.nc
#   data/sst.day.mean.nc
#
# Produces:
#   output/results_ID[FISH ID].rds
#
# FISH ID: Is the true or simulated individual.
#
# Author: James Ounsley

rm(list=ls())
gc()

library(Rcpp)
library(RNetCDF)
library(progress)
library(tidyverse)
library(lubridate)

source('./R/utilities.R')

# CPP implementation of convolution, improves runtime
sourceCpp('./cpp/conv.cpp')

#### CONFIGURATION ####


# Specify fish to use
P_ID <- "2009M116" # MSW
# P_ID <- "2009M053" # 1SW
# P_ID <- "sim_north_2009M053" # Simulation north
# P_ID <- "sim_west_2009M053" # Simulation west

# Save all output, set to false for just the likelihoods etc
SAVE_ALL = TRUE

#### FUNCTIONS ####

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

# Convolution scheme - this depends on current lat and current swim speed 
# (i.e body length)
# Requires index of latitude in SIMLAT, not actual latitude
H_unif = function(lat_index,p_length){
  
  # Get the lat
  p_lat = SIM_LAT[lat_index+1]
  
  # We want KM per DAY 
  # so scale BL (m) to KMs
  d_km = SWIM_SPEED_PBL * (p_length/1000)
  
  d_lat = SWIM_SPEED_PBL * (p_length/1000) / (DEG_NS_TO_KM)
  d_lon = SWIM_SPEED_PBL * (p_length/1000) / (DEG_EW_TO_KM * cos(p_lat*pi/180))
  d_lat = max(d_lat - d_lat%%SIM_RES,SIM_RES)
  if(d_lon < SIM_RES){
    d_lon = SIM_RES
  }else{
    d_lon = max(d_lon - d_lon%%SIM_RES,SIM_RES)
    #d_lon
  }
  
  
  # # Construct a uniform convolution scheme
  # # [0] [0] [x] [0] [0]
  # # [F] [x] [x] [x] [0]
  # # [x] [x] [x] [x] [x]
  # # [0] [x] [x] [x] [0]
  # # [0] [0] [x] [0] [0]
  H_out = expand.grid(seq(-d_lon,d_lon, by=SIM_RES),
                      seq(-d_lat,d_lat, by=SIM_RES))
  
  H_out = (H_out$Var1*DEG_EW_TO_KM*cos(p_lat*pi/180))^2 + (H_out$Var2*DEG_NS_TO_KM)^2 < d_km^2
  #H_out = abs(H_out$Var2) <= (d_lon - abs(H_out$Var1)) * d_lat/d_lon
  H_out = matrix(H_out,nrow =  length(seq(-d_lat,d_lat,by=SIM_RES)),byrow = TRUE)
  H_out = H_out/sum(H_out)
  return(H_out)
}


#### DATA ####

# Define the conversion from SST / LSMASK to our coords
domain = expand.grid(lat_index(SIM_LAT),lon_index(SIM_LON))
domain = cbind(domain$Var2,domain$Var1)

# NOAA OISST Data
# NOAA High Resolution SST data provided by the NOAA/OAR/ESRL PSL, Boulder, Colorado, USA, 
# from their Web site at 
# https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html

# Land sea mask for NOAA OISST
lm_con = open.nc('./data/lsmask.oisst.v2.nc')
lsmask = var.get.nc(lm_con,'lsmask')
close.nc(lm_con)

# Restrict the mask to only the domain of interest
mask = lsmask[domain]
mask = matrix(mask,nrow = length(SIM_LAT))

# Open noaa sst data connection, need to read day by day in main loop to save
# on memory
sst_data_con = open.nc('./data/sst.day.mean.nc')
times = var.get.nc(sst_data_con,'time')
sst_lon <- var.get.nc(sst_data_con,'lon')
sst_lat <- var.get.nc(sst_data_con,'lat')

hmm_data <- readRDS('./data/hmm_data.rds')

# Get data for specified fish
hmm_data <- hmm_data %>% filter(id==P_ID)


#### INITIALISE ####

# Number of timesteps
TS = nrow(hmm_data)

# Not every sim_time will be in the otolith temperature values,
# i.e we do not have data for every time step of the mode. Derive
# subset of timesteps with data
data_times = !is.na(hmm_data$sst)

start_data <- hmm_data %>% filter(datetime==min(datetime))
end_data <- hmm_data %>% filter(datetime==max(datetime))

# Create object to store likelihoods and intialise first time step to 
# release point
phi = array(0,dim=c(length(SIM_LAT),length(SIM_LON),TS))
phi[which(SIM_LAT==start_data$lat),which(SIM_LON==start_data$lon),1] = 1

# Define likelihood for final timestep (catch point)
L_final = array(0,dim=c(length(SIM_LAT),length(SIM_LON)))
L_final[which(SIM_LAT==end_data$lat),which(SIM_LON==end_data$lon)] = 1

# Store log likelihood, max likelihood, min temperature difference and temp
# at max likelihood
# NOTE other than log likelihood these are indicative only as are calculated
# in forward pass before smoothing in backward pass
max_likelihood = log_likelihood = array(NA,dim = TS)
max_likelihood[1] = log_likelihood[1] = 0
min_temp_diff = array(NA,dim = TS)
temp_max_ll = array(NA,dim = TS)

# Object to store only the movement, not likelihood update
process_step = phi

#### PROCESS ####

print('Running forward filter...')

pb <- progress_bar$new(total = length(seq(2,TS)),
                       format = "[:bar] :percent in :elapsed",
                       clear = FALSE)

# FORWARD SWEEP
# Slows down as simulation progresses
# Run over timesteps 
for(t in seq(2,TS)){
  pb$tick()
  
  # Current fish length
  p_length = hmm_data$length[t]/100
  
  # Process step
  # Update position with convolution process, given convolution scheme defined
  # by current body length
  
  # phi(t_{k+1} | t_k) = H(t_k | t_{k+1}) * phi(t_{k} | t_{k})
  p = phi[,,t-1]
  p = conv_cpp(H_unif,p,p_length)
  
  # Saved for the backward filter
  # phi(t_{k+1} | t_{k})
  process_step[,,t] = p
  
  # Set position probs to 0 on land
  p[mask==0] = 0
  
  # Observation step
  # Use sst values and observed temperature data to generate posterior on 
  # positions given data
  # Only do this when data is available for this time
  # step (or it is the last timestep).
  if(data_times[t]==TRUE || t==TS || t==1){
    # Read data for this date from NETCDF file
    sst_vals = var.get.nc(sst_data_con,'sst',
                          start = c(NA,NA,which(times==hmm_data$datetime[t])), 
                          count = c(NA,NA,1))
    # Restrict to our domain
    sst_vals = sst_vals[domain]
    sst_vals = matrix(sst_vals,nrow = length(SIM_LAT))
    
    if(t==TS){
      # If we are in the final time step, we know the likelihood
      L = L_final
      log_likelihood[t] = log(sum(L))
    }else{
      
      temp_val <- hmm_data$sst[t]
      
      # Likelihood of each sst value from normal distribution given
      # mean of otolith temperature and SD
      L = dnorm(sst_vals,temp_val, sd = SD)
      
      # Store minimum temperature difference
      min_temp_diff[t] = min(abs(sst_vals-temp_val),na.rm = TRUE)
    }
    
    # Combine likelihood and position
    # phi(t_{k+1} | t_{k+1}) = L(t_{k+1}) x phi(t_{k+1} | t_k) / lambda(t_{k+1})
    p = p*L
    
    # Set probs to 0 on land
    p[mask==0] = 0
    
    # Normalise
    #phi[,,t] = p/sum(p)
    
    # Store values of interest
    # Total likelihood at this timestep, used for max likelihood parameter
    # search by summing all values, do not require smoothing step. 
    # See Thygesen et al. 2009 Likelihood estimate of parameters
    log_likelihood[t] = log(sum(p))
    
    # Non smoothed estimates of interest
    max_likelihood[t] = log(max(p))
    temp_max_ll[t] = sst_vals[which.max(p)]
    
  } # No data for this time step, ignore observation
  
  # Normalise
  phi[,,t] = p/sum(p)
  
  
}



# Smoothing recursion - backward sweep
# similar run time to forward sweep
#
# Derive posterior distribution of states conditional on all
# data (forward and backward)
print('Running backward smoothing filter...')
pb = progress_bar$new(total = length((TS-1):2),
                      format = "[:bar] :percent in :elapsed",
                      clear = FALSE)

# Structure of smoothed values equal to unsmoothed
# Final time step already conditional on all data
rho = phi

# Go backwards through time steps
for(t in (TS-1):2){
  pb$tick()
  p_length = hmm_data$length[t+1]/100
  
  # phi(t_{k+1} | \inf) / phi(t_{k+1} | t_{k})
  p = rho[,,t+1]/process_step[,,t+1]
  p[is.nan(p)] = 0
  
  # K == H as convolution scheme is symmetrical
  # phi(t_k | t_k) x [K(t_k,t_k{+1} * phi(t_{k+1} | inf) / phi(t_{k+1} | t_{k})]
  rho[,,t] = phi[,,t]*conv_cpp(H_unif,p,p_length)
  rho[,,t][is.nan(rho[,,t])] = 0
  
  
  rho[,,t] = rho[,,t]/sum(rho[,,t])
}


#### SAVE ####

# Do we want to save the raw results
save_phi = save_rho = save_process_step = NULL
if(SAVE_ALL){
  save_phi = phi
  save_rho = rho
  save_process_step = process_step
}

# Save meta data and results
results = list(
  # Meta data
  P_ID = P_ID,
  SD = SD,
  BL_PS = BL_PS,
  TS = TS,
  DAYS_PER_TS = DAYS_PER_TS,
  SIM_LAT=SIM_LAT,
  SIM_LON=SIM_LON,
  times=hmm_data$datetime,
  data_times = data_times,
  
  # Aggregate statistics
  log_likelihood = log_likelihood,
  max_likelihood = max_likelihood,
  min_temp_diff = min_temp_diff,
  temp_max_ll = temp_max_ll,
  
  # Raw results
  likelihoods = save_phi,
  likelihoods_sm = save_rho,
  process_step = save_process_step
)

if(!dir.exists('./output/')) dir.create('./output/')
saveRDS(results,file = paste0('./output/results_ID',P_ID,'.rds'))







