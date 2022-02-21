###########################################################################################
################################# Plot maps for ###########################################
######################### Change in human population density ##############################
###################################### and ################################################
######################## Number of vertebrate extinctions #################################
###########################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Created: Jan 2022
# Last updated: Feb 2022

# clear workspace
rm(list=ls())


##########################################################################################
################################### Set up ###############################################


################################### User inputs ##########################################


# set time name to select dataset for analysis
# must match a time_name in prep_data.R
time_name <- "1700-2000_100y"

# set time period to plot
time <- 1900

# set paths to data and outputs
path_data <- "~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/"


####################################### load packages ####################################


library("sf")
library("rnaturalearth")
library("dplyr")
library("tidyr")
#install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
library("rnaturalearthhires")
library("ggplot2")
theme_set(theme_bw())
theme_update(axis.line=element_blank(),
             axis.text.x=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks=element_blank(),
             axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             panel.background=element_blank(),
             panel.border=element_blank(),
             panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank())
library("gridExtra")
library("SpatialPack") # for modified t test


##################################### Load data ########################################


# import extinction and population density data
popden_country <- readRDS(paste0(path_data, time_name, "/", time_name, "_popden_by_country.rds"))
vertex <- readRDS(paste0(path_data, time_name, "/", time_name, "_vertex.rds"))

# load world data
world <- ne_countries(scale = 'large', returnclass = 'sf')


########################################################################################
#################################### Plots #############################################


# mergewith map data for plotting
popden_world <- merge(popden_country, world, by.x="Entity", by.y="name")

# filter by specified time period
popden_world_time <- popden_world[which(popden_world$Year == time),]
vertex_time <- vertex[which(vertex$Year_Block_Var == time),]

# make into mapping objects
popden_world_time_sf <- st_as_sf(popden_world_time)
vertex_time_sf <- st_as_sf(vertex_time) 

# find centroids to plot extinctions as points
vertex_time_sf_cent <- st_centroid(vertex_time_sf) #cent <- st_centroid(grid_count_sub)

# plot human population density and extinctions on one map
ggplot() + 
  geom_sf(data=popden_world_time_sf, aes(fill=log1p(Prop_change)), size = 0.2) +
  scale_fill_continuous(type = "gradient", low="yellow", high="red") +
  geom_sf(data = world, color='grey', fill=NA, size=0.2) +
  geom_sf(data = vertex_time_sf_cent, size=1) +
  labs(fill = "Population Density Change (log)")


#########################################################################################
################################# Spatial Correlation ###################################

##### Unfinished

# run Pearson's spatial correlation to see if number of extinctions and change in human population
# density is correlated spatially

# find centroids for popden data
cent_popden <- st_centroid(popden_world_time_sf)

# how to compare for different numbers of records?
#correlation <- modified.ttest(popden_world$NoExSpec, popden_world$PopDenChange, popden_world_xy, nclass = 10)


## end of script
