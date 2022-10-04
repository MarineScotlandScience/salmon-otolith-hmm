# Code to generate SST values from known migration
# trajectories
#
# Author: James Ounsley

rm(list=ls())
gc()

library(RNetCDF)
library(tidyverse)
library(lubridate)
library(mapdata)

source('./R/utilities.R')

# Seed for simulations in the ms: 111
set.seed(111)


#### CONFIG ####

# Specify fish to use for length data
P_ID <- "2009M053"

# Name for simulation
SIM_TYPE <- "sim_west"
SIM_TYPE <- "sim_north"
DIRECTION <- ifelse(grepl("sim_north",SIM_TYPE),"N","W")

P1 <- 120 # Length of phase 1 in days (outward migration)
P2 <- 230 # Length of phase 2 in days (perpendicular movement)

#### LOAD DATA ####

# Land sea mask for NOAA OISST
lm_con = open.nc('./data/lsmask.oisst.v2.nc')
lsmask = var.get.nc(lm_con,'lsmask')
close.nc(lm_con)

sst_data_con = open.nc('./data/sst.day.mean.nc')
times <- var.get.nc(sst_data_con,'time')
sst_lon <- var.get.nc(sst_data_con,'lon')
sst_lat <- var.get.nc(sst_data_con,'lat')

# Interpolated length data
# Defines body length at time and start and end dates
length_data = read_csv('./data/otolith_sst.csv') %>% filter(id==P_ID)

#### FUNCTIONS ####

# Wrapper to get sst from global env values
get_sst <- function(t){
  mean_sst <- var.get.nc(sst_data_con,'sst',
                         start = c(which(sst_lon == snap_to_grid(lon[t]%%360)),
                                   which(sst_lat == snap_to_grid(lat[t])),
                                   which(times == length_data$datetime[t])), 
                         count = c(1,1,1))
  
  return(mean_sst)
  
}

get_new_lat_lon <- function(x_km,y_km,old_lat,old_lon){
  # Equitorial projection method (low accuracy, fast performance)
  new_lat <- y_km / R_km + old_lat*pi/180
  new_lat <- new_lat*180/pi
  
  new_lon <- x_km / (R_km*cos(old_lat*pi/180)) + old_lon*pi/180
  new_lon <- new_lon*180/pi
  
  return(list(lat=new_lat,lon=new_lon))
}

get_distance_km <- function(d){
  # n2
  #runif(1) * BL_PS * length_data$length[d]/100 * 24 * 60 * 60 / 1000
  
  # Random uniform scaling of BL_PS times length scaled to km per day
  runif(1) * BL_PS * length_data$length[d]/100 * 24 * 60 * 60 / 1000
}

#### INITIALISE ####

# Initial and final locations
# North coast (Melvich) origin
origin_lat <- snap_to_grid(58.5)
origin_lon <- snap_to_grid(-4.24)

N <- nrow(length_data)

lat <- lon <- sst_val <- c()
lat[1] <- origin_lat
lon[1] <- origin_lon
sst_val[1] <- get_sst(1)

lat[N] <- origin_lat
lon[N] <- origin_lon
sst_val[N] <- get_sst(N)


#### PROCESS ####


if(DIRECTION=="N"){
  #### NORTHWARD MIGRATION ####
  
  x_direction <- 1
  
  for(d in 2:(N-1)){
    
    old_lat <- lat[d-1]
    old_lon <- lon[d-1]
    
    
    if(d <= P1){
      # Phase 1
      # If within first P1 days move north
      
      # Move north, daily distance
      # km per day
      y_km <- get_distance_km(d)
      x_km <- 0
      
      
    }else if(d <= (P1 + P2)){
      # Phase 2
      # Move longitudinally for P2 days
      
      # Move along x axis within coords
      y_km <- 0
      x_km <- x_direction * get_distance_km(d)
      
      # Check where we would end up
      new_lat_lon <- get_new_lat_lon(x_km,y_km,old_lat,old_lon)
      
      # Change direction if we exceed longitudinal boundaries
      if(new_lat_lon$lon >  12) x_direction <- -1
      if(new_lat_lon$lon < -12) x_direction <-  1
      
      # Set direction
      x_km <- x_direction * get_distance_km(d)
      
    }else{
      # d > P1 + P2
      # Phase 3
      # Move toward stationary location directly North of start point
      
      if(lon[d-1] < stat_lon) x_direction <- 1
      if(lon[d-1] > stat_lon) x_direction <- -1
      
      y_km <- 0
      x_km <- x_direction * get_distance_km(d)
      
      # Check where we would end up
      new_lat_lon <- get_new_lat_lon(x_km,y_km,old_lat,old_lon)
      
      # Check if we would overshoot
      if(x_direction > 0 & new_lat_lon$lon > stat_lon) x_direction <- 0
      if(x_direction < 0 & new_lat_lon$lon < stat_lon) x_direction <- 0
      
      # Have we reached the stationary point?
      if(x_direction == 0){
        # Populate remaining data from stationary point
        lat[d:(N-1)] <- stat_lat
        lon[d:(N-1)] <- stat_lon
        sst_val[d:(N-1)] <- sapply(d:(N-1),FUN = get_sst)
        
        # End loop
        break
      }
    }
    
    # Update position given movement values
    new_lat_lon <- get_new_lat_lon(x_km,y_km,old_lat,old_lon)
    lat[d] <- new_lat_lon$lat
    lon[d] <- new_lat_lon$lon
    
    # Save sst at this position
    sst_val[d] = get_sst(d)
    
    if(d==P1){
      # Save the stationary point at end of P1 duration
      stat_lat <- lat[d]
      stat_lon <- lon[d] 
    }
    
    
  }
  
  # Phase 4
  # Move backward from end point to stationary point
  for(d in (N-1):2){
    old_lat <- lat[d+1]
    old_lon <- lon[d+1]
    
    # Daily distance
    # km per day
    y_km <- get_distance_km(d)
    x_km <- 0
    
    new_lat_lon <- get_new_lat_lon(x_km,y_km,old_lat,old_lon)
    
    lat[d] <- new_lat_lon$lat
    lon[d] <- new_lat_lon$lon
    
    # Have we exceeded statinoary point?
    if(lat[d] > stat_lat){
      lat[d] = stat_lat
      lon[d] = stat_lon
      sst_val[d] = get_sst(d)
      
      # Exit loop
      break
    }else{
      sst_val[d] = get_sst(d)
    }
  }
}

#### WESTWARD MIGRATION ####

if(DIRECTION == "W"){
  y_direction <- 1
  
  for(d in 2:(N-1)){
    
    old_lat <- lat[d-1]
    old_lon <- lon[d-1]
    
    
    if(d <= P1){
      # Phase 1
      # If within first P1 days move north
      
      # Move west, daily distance
      # km per day
      y_km <- 0
      x_km <- -1*get_distance_km(d)
      
      
    }else if(d <= (P1 + P2)){
      # Phase 2
      # Move longitudinally for P2 days
      
      # Move along x axis within coords
      y_km <- y_direction * get_distance_km(d)
      x_km <- 0
      
      # Check where we would end up
      new_lat_lon <- get_new_lat_lon(x_km,y_km,old_lat,old_lon)
      
      # Change direction if we exceed longitudinal boundaries
      if(new_lat_lon$lat > 62.5) y_direction <- -1
      if(new_lat_lon$lat < 55) y_direction <-  1
      
      # Set direction
      y_km <- y_direction * get_distance_km(d)
      
    }else{
      # Phase 3
      # Move toward stationary location directly North of start point
      
      if(lat[d-1] < stat_lat) y_direction <- 1
      if(lat[d-1] > stat_lat) y_direction <- -1
      
      y_km <- y_direction * get_distance_km(d)
      x_km <- 0
      
      # Check where we would end up
      new_lat_lon <- get_new_lat_lon(x_km,y_km,old_lat,old_lon)
      
      # Check if we would overshoot
      if(y_direction > 0 & new_lat_lon$lat > stat_lat) y_direction <- 0
      if(y_direction < 0 & new_lat_lon$lat < stat_lat) y_direction <- 0
      
      # Have we reached the stationary point?
      if(y_direction == 0){
        # Populate remaining data from stationary point
        lat[d:(N-1)] <- stat_lat
        lon[d:(N-1)] <- stat_lon
        sst_val[d:(N-1)] <- sapply(d:(N-1),FUN = get_sst)
        
        # End loop
        break
      }
    }
    
    # Update position given movement values
    new_lat_lon <- get_new_lat_lon(x_km,y_km,old_lat,old_lon)
    lat[d] <- new_lat_lon$lat
    lon[d] <- new_lat_lon$lon
    
    # Save sst at this position
    sst_val[d] = get_sst(d)
    
    if(d==P1){
      # Save the stationary point at end of P1 duration
      stat_lat <- lat[d]
      stat_lon <- lon[d] 
    }
    
    
  }
  
  # Phase 4
  # Move backward from end point to stationary point
  for(d in (N-1):2){
    old_lat <- lat[d+1]
    old_lon <- lon[d+1]
    
    # Daily distance
    # km per day
    y_km <- 0
    x_km <- -1*get_distance_km(d)
    
    new_lat_lon <- get_new_lat_lon(x_km,y_km,old_lat,old_lon)
    
    lat[d] <- new_lat_lon$lat
    lon[d] <- new_lat_lon$lon
    
    # Have we exceeded statinoary point?
    if(lon[d] < stat_lon){
      lat[d] = stat_lat
      lon[d] = stat_lon
      sst_val[d] = get_sst(d)
      # Exit loop
      break
    }else{
      sst_val[d] = get_sst(d)
    }
  }
}
#### SAVE RESULTS ####
sim_data <- tibble(lat=lat,lon=lon,sst=sst_val,
                   datetime=length_data$datetime,
                   length=length_data$length,
                   id=paste0(SIM_TYPE,"_",P_ID))
sim_data <- sim_data %>% 
  mutate(sst_true = sst,
         sst = rnorm(n(),sst_true,SD))

if(!dir.exists('./sim/')) dir.create('./sim/')
saveRDS(sim_data,paste0('./sim/',sim_data$id[1],'.rds'))


#### OUTPUT VISUALISATION ####

# Map for plotting
w2hr <- map_data("world")

# Plot
tibble(lat=lat,lon=lon,sst=sst_val) %>% 
  ggplot(aes(x=lon,y=lat,col=sst_val)) +
  
  # Add the map
  geom_polygon(data=w2hr,
               aes(x = long, y = lat, group = group),
               fill='grey',
               inherit.aes = FALSE) +
  
  coord_quickmap(xlim = c(-75,25),
                 ylim = c(45,80)) +
  
  # Plot contours for CI specified by prob_cut
  geom_point() +
  
  # Add the theme
  theme_bw() + 
  xlab('Longitude') + ylab('Latitude') 

# # Animate (without map)
tibble(lat=lat,lon=lon,sst=sst_val,frame=1:N) %>%
  ggplot(aes(x=lon,y=lat,col=sst)) +
  coord_quickmap(xlim = c(-75,25),
                 ylim = c(45,80)) +
  
  # Plot contours for CI specified by prob_cut
  geom_point() +
  
  # Add the theme
  theme_bw() +
  xlab('Longitude') + ylab('Latitude') +
  gganimate::transition_manual(frames = frame)


