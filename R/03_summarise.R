# Summarise HMM
#
# Create files needed to calculate summary statistics for the HMM fits, and where
# data is from a simulation use the true values to scrutinise the fit.
# 
# Requires:
#   output/results_ID[FISH ID].rds
#   data/lsmask.oisst.v2.nc  
#   data/sst.day.mean.nc
#   data/hmm_data.rds
#
# Produces:
#   output/plt_data.RData
#
rm(list=ls())
library(tidyverse)
library(lubridate)
library(mapdata)
library(gganimate)
library(patchwork)

source('R/utilities.R')

# Load universal data ####
w2hr <- map_data("world")

# SST data
# Define the conversion from SST / LSMASK to our coords
domain = expand.grid(lat_index(SIM_LAT),lon_index(SIM_LON))
domain = cbind(domain$Var2,domain$Var1)

# Land sea mask for NOAA OISST
lm_con = open.nc('./data/lsmask.oisst.v2.nc')
lsmask = var.get.nc(lm_con,'lsmask')
close.nc(lm_con)

# Restrict the mask to only the domain of interest
mask = lsmask[domain]
mask = matrix(mask,nrow = length(SIM_LAT))

sst_data_con = open.nc('./data/sst.day.mean.nc')
times = var.get.nc(sst_data_con,'time')

# Specify data to evaluate ####
P_ID <- c("2009M116", "2009M053", "sim_west_2009M053", "sim_north_2009M053")

# Start loop ####

big_track_tbl <- list()
big_contour_tbl <- list()
big_sst_tbl <- list()
big_history_tbl <- list()

for(i in 1:length(P_ID)){
  
  SIM <- grepl("sim",P_ID[i])
  
  # Load specific data
  results <- readRDS( paste0('./output/results_ID',P_ID[i],'.rds')) 
  results$date <- as_date(results$times,origin = "1800/1/1")
  
  sim_data <- readRDS('./data/hmm_data.rds') %>% 
    filter(id==P_ID[i]) %>% 
    mutate(date = as_date(datetime,origin = "1800/1/1"))
  
  
  # Init ####
  N <- length(results$data_times)
  data_indices <- 1:N
  master_track_tbl <- list()
  master_contour_tbl <- list()
  master_history_tbl <- list()
  master_sst_tbl <- list()
  
  # Loop
  for(t in data_indices){
    
    # Extract the likelihoods for this timestep
    surface <- results$likelihoods_sm[,,t]
    colnames(surface) <- results$SIM_LON
    
    surface_tbl <- as_tibble(surface) %>%
      # Add latitude
      mutate(lat = results$SIM_LAT) %>%
      # Convert to long form
      gather(-lat, key = 'lon',value = 'lik') %>%
      mutate(lon = as.numeric(lon)) 
    
    # Mean locations ####
    
    # From Pederson PhD
    #
    # Mean of smoothed distribution
    # \bar{x}_k = \sum_x x [\sum_y \phi_{x,y}(t_k,Z_N) ]
    # \bar{y}_k = \sum_y y [\sum_x \phi_{x,y}(t_k,Z_N) ]
    #
    # \phi_{x,y}(t_k,Z_N) is the smoothed filtered posterior
    lon_mean <- surface_tbl %>% group_by(lon) %>% 
      summarise(summed_lik = sum(lik)) %>% 
      mutate(tmp = lon*summed_lik) %>% 
      ungroup() %>% 
      summarise(lon_mean = sum(tmp)) %>% 
      pull(lon_mean)
    
    lat_mean <- surface_tbl %>% group_by(lat) %>% 
      summarise(summed_lik = sum(lik)) %>% 
      mutate(tmp = lat*summed_lik) %>% 
      ungroup() %>% 
      summarise(lat_mean = sum(tmp)) %>% 
      pull(lat_mean)
    
    # SST values ####
    # Read data for time t from SST NETCDF file
    sst_vals = var.get.nc(sst_data_con,'sst',
                          start = c(NA,NA,which(times==results$times[t])),
                          count = c(NA,NA,1))
    
    # Restrict to our domain and create data frame
    sst_vals = sst_vals[domain]
    sst_vals = matrix(sst_vals,nrow = length(SIM_LAT))
    dimnames(sst_vals) <- list(SIM_LAT,SIM_LON)
    sst_vals <- as.data.frame.table(sst_vals)
    colnames(sst_vals) <- c("lat","lon","sst")
    sst_vals$lat <- as.numeric(as.character(sst_vals$lat))
    sst_vals$lon <- as.numeric(as.character(sst_vals$lon))
    
    # get SST value in mean lat / lon location 
    sst_mean <- sst_vals %>%
      filter(lat == snap_to_grid(lat_mean) & lon == snap_to_grid(lon_mean)) %>%
      pull(sst)
    
    # get all SST values within 68% CI of the simulation SST (ca 1SD)
    sst_match <- sst_vals %>%
      filter(sst > qnorm(0.16,as.numeric(sim_data[t,"sst"]),1) &
               sst < qnorm(0.84,as.numeric(sim_data[t,"sst"]),1)) %>% # restrict to SST range (+/- 1 SD)
      mutate(diff = sst - as.numeric(sim_data[t,"sst"]),
             P_ID = results$P_ID,
             data_week = t,
             no_data = !results$data_times[t],
             date = as.Date(paste(results$date[t])),
             datetime = results$time[t])
    
    master_sst_tbl[[t]] <- sst_match
    
    # Summary values including mean track
    track_tbl <- list(lon_mean = lon_mean, lat_mean = lat_mean,
                      P_ID = results$P_ID,
                      data_week = t,
                      no_data = !results$data_times[t],
                      date = as.Date(paste(results$date[t])),
                      datetime = results$time[t],
                      sst_mean = sst_mean)
    
    
    # Distance from true values(only populated for SIM = TRUE) ####
    
    sim_lon = sim_data %>% filter(date == results$date[t]) %>% pull(lon)
    sim_lat = sim_data %>% filter(date == results$date[t]) %>% pull(lat)
    
    # Pythag on equatorial proj.
    # x = delta long cos(lat)
    # y = delta lat
    # d = R sqrt(x2 + y2)
    surface_tbl <- surface_tbl %>% 
      mutate(delta_lat = (sim_lat-lat)*pi/180,
             delta_lon = (sim_lon-lon)*pi/180,
             x = delta_lon * cos(sim_lat*pi/180),
             y = delta_lon,
             d = R_km*sqrt(x^2 + y^2))
    
    d_err <- surface_tbl %>% 
      mutate(w_d = lik*d/sum(lik)) %>% 
      summarise(d_err = sum(w_d)) %>% 
      pull(d_err)
    
    # Likelihood ####
    
    # Note that these values are not the likelihoods over the timeseries,
    # they have been normalised to sum to 1 within each timestep, so 
    # they give the probability of the sim value at this timestep only.
    ll_true <- surface_tbl %>% 
      filter(lat==sim_lat,
             lon==sim_lon) %>% 
      pull(lik)
    
    
    track_tbl$km_err = d_err
    track_tbl$ll_true = ifelse(length(ll_true) > 0, ll_true, NA)
    
    # Cumulative likelihood for heatmap plot ####
    
    contour_tbl <- surface_tbl %>%
      # Calculate cumulative density from high to low prob
      arrange(desc(lik)) %>%
      mutate(tot_lik = cumsum(lik)) %>%
      # For .1 and .01 contours need a bit more
      #filter(tot_lik <= PROB_CUT + .5) %>%
      mutate(P_ID = results$P_ID,
             data_week = t,
             no_data = !results$data_times[t],
             date = as.Date(paste(results$date[t])),
             datetime = results$time[t])
    
    # liklihood value at threshold of 90% total probability 
    # to show if true location is within 90% of inferred locations
    thresh <- tail(contour_tbl %>%
                     filter(tot_lik<=0.9),1) %>%
      pull(lik)
    in_ci <- ifelse(length(thresh) > 0, ll_true > thresh, NA) # assign NA to first and last position
    track_tbl$in_ci = in_ci
    
    master_track_tbl[[t]] <- track_tbl  
    master_contour_tbl[[t]] <- contour_tbl
    
    # For contour animation generate history
    tmp <- tibble(lat_mean = lat_mean, lon_mean = lon_mean)
    # If there is a history, add the rows
    if(t>1){
      tmp <- tmp %>% 
        bind_rows(master_history_tbl[[t-1]] %>% select(lat_mean,lon_mean))
    }
    tmp <-  tmp %>% 
      mutate(date = as.Date(paste(results$date[t])),
             alpha = N:(N-(n()-1))/N)
    
    master_history_tbl[[t]] <- tmp
    
  }
  
  
  # Add sea age and simulation categories
  master_track_tbl <- bind_rows(master_track_tbl) %>% 
    left_join(sim_data %>% select(datetime,lat,lon,sst,sst_true,length,sea_age),by='datetime') %>%
    mutate(sim = grepl("sim",P_ID[i]))
  master_contour_tbl <- bind_rows(master_contour_tbl) %>%
    mutate(sea_age = first(sim_data$sea_age),
           sim = grepl("sim",P_ID[i]))
  master_sst_tbl <- bind_rows(master_sst_tbl) %>%
    mutate(sea_age = first(sim_data$sea_age),
           sim = grepl("sim",P_ID[i]))
  master_history_tbl <- bind_rows(master_history_tbl) %>%
    mutate(sea_age = first(sim_data$sea_age),
           sim = grepl("sim",P_ID[i]))
  
  # gather all data together ####
  big_track_tbl[[i]] <- master_track_tbl
  big_contour_tbl[[i]] <- master_contour_tbl
  big_sst_tbl[[i]] <- master_sst_tbl
  big_history_tbl[[i]] <- master_history_tbl
  
}

big_track_tbl <- bind_rows(big_track_tbl) %>%
  group_by(P_ID) %>%
  mutate(time = reorder(factor(format(date, "%b-%y")), date))
big_contour_tbl <- bind_rows(big_contour_tbl) %>%
  group_by(P_ID) %>%
  mutate(time = reorder(factor(format(date, "%b-%y")), date))
big_sst_tbl <- bind_rows(big_sst_tbl) %>%
  group_by(P_ID) %>%
  mutate(time = reorder(factor(format(date, "%b-%y")), date))

big_history_tbl <- bind_rows(big_history_tbl)

if(!dir.exists('./output/')) dir.create('./output/')
save(big_track_tbl, big_contour_tbl, big_sst_tbl, big_history_tbl,file = './output/plot_data.RData')


