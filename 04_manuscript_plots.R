# Create manuscript plots
#
# Requires:
#   output/plot_data.RData.rds
#
# Produces:
#   plt/
#

rm(list=ls())
library(tidyverse)
library(lubridate)
library(mapdata)
library(gganimate)
library(patchwork)

load('./output/plot_data.RData')

# Summary stats ################################################################  
  # Visualise error over time
  big_track_tbl %>% 
  filter(sim == TRUE) %>%
    ggplot(aes(x=data_week,y=km_err)) + 
    geom_line() + 
    geom_point(aes(col=no_data)) + 
    labs(title = paste('Mean km err:',round(mean(master_track_tbl$km_err)))) +
    facet_wrap(~P_ID)
  
  big_track_tbl %>% 
    filter(sim == TRUE) %>%
    group_by(P_ID) %>%
    summarise(mean_err_km = mean(km_err),
              sum_ll = sum(log(ll_true)),
              median_err_km = median(km_err),
              q_25 = quantile(km_err,0.25),
              q_75 = quantile(km_err,0.75),
              IQR = q_75 - q_25,
              in_ci = sum(in_ci, na.rm = T)/sum(!is.na(master_track_tbl$in_ci)))

  
# Mean track plots #############################################################
big_track_tbl$lab <- with(big_track_tbl, ifelse(grepl("north",P_ID),"Northwards",ifelse(grepl("west",P_ID),"Westward",sea_age)))
  
plt_mean_fun <- function(track_data) {
plt_mean <- track_data %>%
  ggplot(aes(x=lon_mean,y=lat_mean,col=date), size = 2) +
  # Add the map
  geom_polygon(data=w2hr,
               aes(x = long, y = lat, group = group),
               fill='grey',
               inherit.aes = FALSE) +
  coord_quickmap(xlim = c(-75,25),
                 ylim = c(45,80))

if(track_data$sim){
  # Add true position from simulation
  plt_mean <- plt_mean + 
    geom_point(aes(lon,lat,alpha=date),col='grey30',inherit.aes = FALSE, size = 2)
}

plt_mean <- plt_mean +
  # Plot mean positions
  geom_point(size = 2) +

  # Add the theme
  scale_color_gradient(low = 'gold', high = 'darkorange3') +
  theme_bw() + theme(legend.position = 'none',
                     plot.margin = unit(c(0,0.1,0,0.1),"lines"),
                     panel.background = element_rect(fill = 'white'),
                     panel.grid.major = element_line(colour = "white"),
                     panel.grid.minor = element_line(colour = "white"),
                     aspect.ratio = 1/1) +
  
  xlab('Longitude') + ylab('Latitude') +
  facet_wrap(~lab)

plt_mean

}

(fig1b <- plt_mean_fun(filter(big_track_tbl, sim == TRUE) ) )
(fig3c <- plt_mean_fun(filter(big_track_tbl, sim == FALSE) ) )

# SST profiles #################################################################
plt_sst_fun <- function(track_data) {
  
plt_sst <- track_data %>%
  ggplot(aes(x = date, y = sst)) 
  
  if(track_data$sim){
    # Add true sst from simulation
  plt_sst <- plt_sst + geom_line(aes(date, sst_true),col='grey80',inherit.aes = FALSE, size = 2) 
  }

plt_sst <- plt_sst +
  
  # Add sst at hmm inferred location
  geom_point(aes(date, sst_mean), col = 'orange', inherit.aes = FALSE) +
  
  # Add sst from otolith
  geom_point(col = 'darkblue') +

  # Add the theme
  theme_bw() + theme(legend.position = 'none',
                     plot.margin = unit(c(0,0.1,0,0.1),"lines"),
                     panel.background = element_rect(fill = 'white'),
                     panel.grid.major = element_line(colour = "white"),
                     panel.grid.minor = element_line(colour = "white"),
                     aspect.ratio = 1/1,
                     text = element_text(family = 'Times')) +
    
    xlab('') + ylab('SST (C)') +
    # facet_wrap(~lab)+
    scale_x_date(date_labels = "%b")+
  facet_wrap(~lab)

plt_sst

}

(fig1a <- plt_sst_fun(filter(big_track_tbl, sim == TRUE) ) )
(fig3b <- plt_sst_fun(filter(big_track_tbl, sim == FALSE) ) )

# Likelihood distribution maps #################################################
big_contour_tbl$lab <- with(big_contour_tbl, ifelse(grepl("north",P_ID),"Northwards",ifelse(grepl("west",P_ID),"Westward",sea_age)))
big_sst_tbl$lab <- with(big_sst_tbl, ifelse(grepl("north",P_ID),"Northwards",ifelse(grepl("west",P_ID),"Westward",sea_age)))

plt_distrib_fun <- function(contour_data, sst_data, track_data) {
  
plt_distrib <- contour_data %>%
  group_by(datetime) %>%
  filter(tot_lik <= 0.9) %>%
  ggplot(aes(x = lon, y = lat, fill = tot_lik)) +
  
  # Plot distribution of SST within 1C of observed in otolith
  geom_tile(data = sst_data, aes(x = lon, y = lat), fill = 'lightblue', inherit.aes = FALSE) +

  # Plot distribution of 90% CI
  geom_tile() +

  # Add the map
  geom_polygon(data=w2hr,
               aes(x = long, y = lat, group = group),
               fill='grey',
               inherit.aes = FALSE) +
  
  coord_quickmap(xlim = c(-75,25),
                 ylim = c(45,80)) +
  
  # Add the theme
  viridis::scale_fill_viridis(direction = -1) +
  
  theme_bw() +
  theme(legend.position = 'none',
        axis.ticks.length = unit(0, "pt"),
        panel.spacing = unit(0,"pt"),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        plot.margin = margin(0,0,0,0,"pt")) +
  xlab('Longitude') + ylab('Latitude')


if(track_data$sim){
  # Add true location from simulation
  plt_distrib <- plt_distrib + geom_point(data = track_data, aes(x = lon, y = lat), col = 'red', inherit.aes = FALSE) 
}

plt_distrib + facet_wrap(~time)

plt_distrib 

}

# Need to trim to desired facet months
# get index of first week in each month where data exists
tmp <- big_contour_tbl %>%
  filter(no_data == FALSE) %>%
  group_by(P_ID, time) %>%
  mutate(month = month(date)) %>%
  arrange(date) %>%
  filter(row_number() == 1 & data_week > 3) %>%
  ungroup()

# trim to X panels depending on sea age
onesw <- filter(tmp, P_ID == "2009M053")
onesw <- onesw[seq(match(6,onesw$month),match(6,onesw$month)+11,3),]$time
twosw <- filter(tmp, P_ID == "2009M116")
twosw <- twosw[seq(match(6,twosw$month),match(6,twosw$month)+24,3),]$time


plt_distrib_fun(filter(big_contour_tbl, lab == "2SW" & data_week %in% twosw), 
                filter(big_sst_tbl, lab == "2SW" & data_week %in% twosw) )

# # Save
# if(!dir.exists('./plt/')) dir.create('./plt/')
# ggsave(paste0('./plt/mean_path_ID',P_ID,'.png'),plt_mean)

# Contour animation (Supplementary Figs) #######################################

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


