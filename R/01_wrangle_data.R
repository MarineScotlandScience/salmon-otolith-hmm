library(tidyverse)

source('R/utilities.R')

# Read in raw otolith sst values and lengths
otolith_data <- read_csv('data/otolith_sst.csv') %>% 
  # There are instances of multiple data points per day.
  # Take the first
  group_by(id,datetime) %>% 
  summarise_all(first)

# Read in simulation files used in manuscript
sim_data <- list(readRDS('./sim/sim_north_2009M053.rds'),
                 readRDS('./sim/sim_west_2009M053.rds')) %>% 
  bind_rows() %>% 
  # Wrangle for compatibility with with otolith data
  mutate(datetime = as.numeric(datetime),
         date = dmy('01/01/1800')+datetime) %>% 
  ungroup()

# For the simulations we want to thin the sampled sst values to match those
# of 2009M053
thin_times <- otolith_data %>% 
  filter(id == "2009M053") %>% 
  pull(sst) %>% 
  is.na()
sim_data[rep(thin_times,2),]$sst <- NA

# Add the simulations data to the real otolith data
hmm_data <- otolith_data %>% 
  bind_rows(sim_data)
  
hmm_data <- hmm_data %>% 
  # Filter data to weekly values
  mutate(datetime = floor(as.numeric(datetime) / DAYS_PER_TS)*DAYS_PER_TS) %>% 
  group_by(id,datetime) %>% 
  # Summarise weekly values, take first
  summarise(sst = first(sst,order_by = is.na(sst)),
            lat=first(lat,order_by = is.na(lat)),
            lon=first(lon,order_by = is.na(lon)),
            length=first(length)) %>% 
  # Adjust start and end locations to match SST resolution
  mutate(lat = snap_to_grid(lat),
         lon = snap_to_grid(lon)) %>% 
  ungroup()

saveRDS(hmm_data,'data/hmm_data.rds')





