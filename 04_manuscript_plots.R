# Create manuscript plots
#
# Requires:
#   output/results_ID[FISH ID].rds
#
# Produces:
#   Figures 1 - 4, supplementary Figures S1 - S2, Animations S3 - S6 
#
rm(list=ls())
library(tidyverse)
library(lubridate)
library(mapdata)
library(gganimate)

source('R/utilities.R')

# Specify simulated data to evaluate against
# P_ID <- "2009M116" # MSW
# P_ID <- "2009M053" # 1SW
# P_ID = "sim_west_2009M053"
P_ID = "sim_north_2009M053"

SIM <- grepl("sim",P_ID)

# Load data ####
w2hr <- map_data("world")

results <- readRDS( paste0('./output/results_ID',P_ID,'.rds')) 
results$date <- as_date(results$times,origin = "1800/1/1")

sim_data <- readRDS('./data/hmm_data.rds') %>% 
  filter(id==P_ID) %>% 
  mutate(date = as_date(datetime,origin = "1800/1/1"))

# Init ####
N <- length(results$data_times)
data_indices <- 1:N
master_track_tbl <- list()
master_contour_tbl <- list()
master_history_tbl <- list()

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
  
  # Summary values including mean track
  track_tbl <- list(lon_mean = lon_mean, lat_mean = lat_mean,
                    P_ID = results$P_ID,
                    data_week = t,
                    no_data = !results$data_times[t],
                    date = as.Date(paste(results$date[t])),
                    datetime = results$time[t])
  
  if(SIM){
    # Distance from true values ####

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
    track_tbl$ll_true = ll_true
  }
  

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
  in_ci <- ll_true > thresh
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