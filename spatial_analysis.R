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
world <- world %>% select(name, geometry)

st_crs(world)


########################################################################################
#################################### Plots #############################################


# merge with map data for plotting
popden_world <- merge(popden_country, world, by.x="Entity", by.y="name")

# filter by specified time period
popden_world_time <- popden_world[which(popden_world$Year == time),]
vertex_time <- vertex[which(vertex$Year_Block_Var == time),]

# make into mapping objects
popden_world_time_sf <- st_as_sf(popden_world_time)
vertex_time_sf <- st_as_sf(vertex_time) 

# set coordinate system
coord_sys <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

popden_world_time_sf <- st_sf(st_transform(popden_world_time_sf, crs=4326))
vertex_time_sf <- st_sf(st_transform(vertex_time_sf, crs=4326))

# find centroids to plot extinctions as points
vertex_time_sf_cent <- st_centroid(vertex_time_sf) #cent <- st_centroid(grid_count_sub)

# plot human population density and extinctions on one map
pdf(file="map.pdf")
ggplot() + 
  geom_sf(data=popden_world_time_sf, aes(fill=log1p(Prop_change)), size = 0.2) +
  scale_fill_continuous(type = "gradient", low="yellow", high="red") +
  geom_sf(data = world, color='grey', fill=NA, size=0.2) +
  geom_sf(data = vertex_time_sf_cent, size=1) +
  labs(fill = "Population Density Change (log)")
dev.off()



world_map <- readRDS("~/Dropbox/luke/documents/academia/phd/maps/World/moll_grid.rds")

# set Mollweide projection coordinates
coord_sys <- '+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs'

# fix error in map by making sure object is valid
world_map <- st_make_valid(world_map)
# fix projection error due to different versions GDAL
st_crs(world_map) <- coord_sys

# project to specified coordinate system
world_map <- st_sf(st_transform(world_map, crs=4326))
#world_map<- st_make_valid(world_map)

sf::sf_use_s2(FALSE)

# join with map
grid <- st_join(x=world_map,
                y=st_cast(vertex_time_sf) %>% dplyr::select(Taxa, Class.x, Class_Num, EX.Last.seen., Cont.2..Island.1.),
                left=TRUE, join = st_intersects)

# sum by number of species
#grid_count_sub <- grid_sub %>%
#  group_by(gridId) %>%
#  summarise(SpeciesRichness = n_distinct(Taxa, na.rm = TRUE),
#            spsList = paste(Taxa, collapse = ';'),
#            .groups = 'drop')

#ggplot() + 
#  geom_sf(data = grid_count_sub %>%  
#            mutate(SpeciesRichness=ifelse(SpeciesRichness==0, NA, SpeciesRichness)),
#          aes(fill=SpeciesRichness), col=NA, size=0.2) +
#  scale_fill_continuous(type = "gradient", low="yellow", high="red") +
#  geom_sf(data = world_plot, color='black', fill=NA, size=0.2) +
#  labs(fill='Species Richness')


#########################################################################################
################################# Spatial Correlation ###################################

##### Unfinished

# run Pearson's spatial correlation to see if number of extinctions and change in human population
# density is correlated spatially

# find centroids for popden data
cent_popden <- st_centroid(popden_world_time_sf)

# how to compare for different numbers of records?
#correlation <- modified.ttest(popden_world$NoExSpec, popden_world$PopDenChange, popden_world_xy, nclass = 10)




# try to set country for each exinction based on centroids
# we can then total number of extinctions per country at the time period
# to then compare with human pop density data

library("rworldmap")

coords <- st_coordinates(vertex_time_sf_cent$geometry)
# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees
coords2country = function(points)
{  
  countriesSP <- getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # convert our list of points to a SpatialPoints object
  
  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  #setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  # return the ADMIN names of each country
  indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}

vertex_time_sf_cent$Country <- coords2country(coords)
vertex_time_sf_cent_df <- as.data.frame(vertex_time_sf_cent)
vertex_time_country_sum <- vertex_time_sf_cent_df %>%
                           group_by(Country) %>%
                           summarise(n_extinct = n())

# drop NAs from summary
vertex_time_country_sum <- vertex_time_country_sum[!is.na(vertex_time_country_sum$Country),]

# find centroids for popden data
cent_popden <- st_centroid(popden_world_time_sf)

cent_popden$NoEx <- 0
for (i in 1:nrow(cent_popden)) {
  if (cent_popden$Entity[i] %in% vertex_time_country_sum$Country) {
    cent_popden$NoEx[i] <- vertex_time_country_sum$n_extinct[vertex_time_country_sum$Country == cent_popden$Entity[i]]
  }
}

# get country centroid coordinates
country_coords <- as.matrix(st_coordinates(cent_popden))

library("spatialEco")
cc <- crossCorrelation(cent_popden$NoEx, coords=country_coords)
summary(cc)

# run on exact coordinates
ex_site_coords <- as.matrix(st_coordinates(vertex_time_sf_cent))
cc_sites <- crossCorrelation(vertex_time_sf_cent$Taxa, coords=ex_site_coords)
summary(cc_sites)

## end of script
