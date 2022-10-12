# Wrangle data
#
# Combine actual otolith sst data with simulated values to create a single 
# data object. Wrangle columns to meet requirements of HMM.
#
# Requires:
#   data/otolith_sst.csv
#   sim/sim_north_2009M053.rds    - from 00_simulate_migration.R
#   sim/sim_west_2009M053.rds     - from 00_simulate_migration.R
#
# Produces:
#   data/hmm_data.rds
#   data/daily_data.rds
#
library(tidyverse)

source('R/utilities.R')

# Read in raw otolith sst values and lengths
otolith_data <- read_csv('data/otolith_sst.csv') 

# Read in simulation files used in manuscript
sim_data <- list(readRDS('./sim/sim_north_2009M053.rds'),
                 readRDS('./sim/sim_west_2009M053.rds')) %>% 
  bind_rows() %>% 
  # Wrangle for compatibility with with otolith data
  mutate(datetime = as.numeric(datetime),
         date = lubridate::dmy('01/01/1800')+datetime) %>% 
  ungroup()

# For the simulations we want to thin the sampled sst values to match those
# of 2009M053
thin_times <- otolith_data %>% 
  filter(id == "2009M053") %>% 
  pull(sst) %>% 
  is.na()
sim_data[rep(thin_times,2),]$sst <- NA
sim_data$sea_age <- "1SW"


# Filter simulation file based on presence of data only and bind otolith data
hmm_data <- sim_data %>% 
  bind_rows(otolith_data)

hmm_data <- hmm_data %>% 
  # Group data over weeks
  mutate(tmp_datetime = datetime,
         datetime = floor(as.numeric(datetime) / DAYS_PER_TS)*DAYS_PER_TS) %>% 
  group_by(id,datetime) %>% 
  # Create a sorting index for summarise
  # If an sst value is not na in the week, it will be first by this index,
  # otherwise the first value in the week will be first by the index.
  mutate(is_na_sst = is.na(sst)) %>% 
  # Summarise taking values as determined by sorting index
  summarise(sst = first(sst,order_by = is_na_sst),
            lat = first(lat,order_by = is_na_sst),
            lon = first(lon,order_by = is_na_sst),
            sst_datetime = first(tmp_datetime,order_by = is_na_sst),
            sst_true = first(sst_true,order_by = is_na_sst),
            # Do we want the length to be the same date as sst
            # Length and sst are already unlinked, so lets take
            # regular values for the length (median is middle of the week)
            length = median(length)) %>% 
  # Adjust start and end locations to match SST resolution
  mutate(lat = snap_to_grid(lat),
         lon = snap_to_grid(lon)) %>% 
  ungroup()

saveRDS(hmm_data,'data/hmm_data.rds')

# Add the unfiltered simulation data to the real otolith data to create a 
# daily file
daily_data <- sim_data %>% 
  bind_rows(otolith_data)

saveRDS(daily_data,'data/daily_data.rds')


