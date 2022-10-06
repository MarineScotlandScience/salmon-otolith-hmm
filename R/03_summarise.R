# Summarise HMM
#
# Calculate summary statistics for the HMM fits. Where the data is from a
# simulation use the true values to scrutinise the fit.
# 
# Requires:
#   output/results_ID[FISH ID].rds
#
# Produces:
#   summary/hmm_summary_ID[FISH ID].csv
#   plt/mean_path_ID[FISH ID].png
#   anim/contour_ID[FISH ID].gif
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

# Load
w2hr <- map_data("world")

results <- readRDS( paste0('./output/results_ID',P_ID,'.rds')) 
results$date <- as_date(results$times,origin = "1800/1/1")

sim_data <- readRDS('./data/hmm_data.rds') %>% 
  filter(id==P_ID) %>% 
  mutate(date = as_date(datetime,origin = "1800/1/1"))

# Init
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
  
  # Mean locations 
  ################
  
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
    # Distance from true values
    ###########################
    
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
    
    # Likelihood
    ############
    
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
  
  master_track_tbl[[t]] <- track_tbl  
  
  
  # Cumulative likelihood for animation
  # with countours
  #####################################
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


# Summary stats
###############

master_track_tbl <- bind_rows(master_track_tbl) 

if(SIM){
  master_track_tbl <- master_track_tbl %>% 
    left_join(sim_data %>% select(datetime,lat,lon),by='datetime')
  
  # Visualise error over time
  master_track_tbl %>% 
    ggplot(aes(x=data_week,y=km_err)) + 
    geom_line() + 
    geom_point(aes(col=no_data)) + 
    labs(title = paste('Mean km err:',round(mean(master_track_tbl$km_err))))
  
  master_track_tbl %>% 
    summarise(mean_err_km = mean(km_err),
              sum_ll = sum(log(ll_true)),
              median_err_km = median(km_err),
              q_25 = quantile(km_err,0.25),
              q_75 = quantile(km_err,0.75),
              IQR = q_75 - q_25)
}

if(!dir.exists('./summary/')) dir.create('./summary/')
write_csv(master_track_tbl,
          paste0('./summary/hmm_summary_',P_ID,'.csv'))

# Mean track plot
#################

plt_mean <- master_track_tbl %>%
  ggplot(aes(x=lon_mean,y=lat_mean,alpha=date)) +
  # Add the map
  geom_polygon(data=w2hr,
               aes(x = long, y = lat, group = group),
               fill='grey',
               inherit.aes = FALSE) +
  coord_quickmap(xlim = c(-75,25),
                 ylim = c(45,80))

if(SIM){
  # Add true position from simulation
  plt_mean <- plt_mean + 
    geom_point(aes(lon,lat,alpha=date),col='grey',inherit.aes = FALSE)
}

plt_mean <- plt_mean +
  # Plot mean positions
  geom_point(col='darkblue') +
  geom_path(col='darkblue') +
  
  # Add the theme
  theme_bw() + theme(plot.margin = unit(c(0,0.1,0,0.1),"lines")) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(title = paste('Mean path',P_ID))

plt_mean

# Save
if(!dir.exists('./plt/')) dir.create('./plt/')
ggsave(paste0('./plt/mean_path_ID',P_ID,'.png'),plt_mean)

# Contour animation
###################

master_contour_tbl <- bind_rows(master_contour_tbl)
master_history_tbl <- bind_rows(master_history_tbl)

plt <- master_contour_tbl %>%
  ggplot(aes(x=lon,y=lat,z=tot_lik,linetype=no_data)) +
  
  # Add the map
  geom_polygon(data=w2hr,
               aes(x = long, y = lat, group = group),
               fill='grey',
               inherit.aes = FALSE) +
  
  coord_quickmap(xlim = c(-75,25),
                 ylim = c(45,80)) +
  
  # Plot contours for CI specified by prob_cut
  geom_contour(breaks=c(.01,0.1,0.5,0.9),col='black') +
  
  geom_point(data = master_history_tbl,aes(x=lon_mean,y=lat_mean,alpha=alpha),
             inherit.aes = FALSE,colour='lightblue') +
  geom_point(data = master_history_tbl %>% filter(alpha==1),
             aes(x=lon_mean,y=lat_mean),
             inherit.aes = FALSE,colour='darkblue') +
  
  # Add the theme
  theme_bw() + theme(legend.position = 'none',
                     plot.margin = unit(c(0,0.1,0,0.1),"lines")) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(title = "{current_frame}")

if(SIM){
  plt <-  plt +
    geom_point(data=sim_data,aes(x=lon,y=lat),inherit.aes = FALSE,colour='red')
}

plt <- plt + 
  transition_manual(frames = date,cumulative = FALSE)

anim <- animate(plt,nframes = N,fps=2.5,
                height=9,width=9,units='in',res=300,
                renderer = gifski_renderer())

if(!dir.exists('./anim/')) dir.create('./anim/')
anim_save(paste0('./anim/contour_ID',P_ID,'.gif'),
          animation=anim,renderer = gifski_renderer())


