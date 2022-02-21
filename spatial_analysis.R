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


# load packages
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

# set time name to select dataset for analysis
# must match a time_name in prep_data.R
time_name <- "1700-2000_100y"

# set paths to data and outputs
path_data <- "~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/"
#path_loc <- "~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/species_locations/"

# import extinction and population density data
popden_country <- readRDS(paste0(path_data, time_name, "/", time_name, "_popden_by_country.rds"))
vertex <- readRDS(paste0(path_data, time_name, "/", time_name, "_vertex.rds"))

# import location data
# species locations
#species_loc <- read_sf(paste0(path_loc, "database_SPECIES.shp"))
# world grid for combining with species data
#world_grid <- readRDS("~/Dropbox/luke/documents/academia/phd/maps/World/moll_grid.rds")
# load world data
world <- ne_countries(scale = 'large', returnclass = 'sf')
# world map for plotting
#world_map <- readRDS("~/Dropbox/luke/documents/academia/phd/maps/World/moll.rds")

# set Mollweide projection coordinates
#coord_sys <- '+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs'

# fix error in map by making sure object is valid
#world_map <- st_make_valid(world_map)
# fix projection error due to different versions GDAL
#st_crs(world_map) <- coord_sys

# project to specified coordinate system
#area_projected <- st_transform(st_union(world), crs=coord_sys)

# fix error by making sure object is valid
#area_projected_fix <- st_make_valid(area_projected)

# find area of countries
#world$area_km2 <- st_area(world$geometry)/1e6
#world_sub <- as.data.frame(cbind(world$name, world$area_km2))
#names(world_sub) <- c("name", "area_km2")
#world_sub$area_km2 <- as.numeric(as.character(world_sub$area_km2))


########################################################################################
#################################### Plots #############################################


# set time period to plot
time <- 2000

# mergewith map data for plotting
popden_world <- merge(popden_country, world, by.x="Entity", by.y="name")

# filter by specified time period
popden_world_time <- popden_world[which(popden_world$Year == time),]

# make into mapping object
popden_world_time_sf <- st_as_sf(popden_world_time)

# plot human pop density at specified time period
#ggplot(data=popden_world_time_sf) + 
#  geom_sf(aes(fill=log1p(Prop_change)), size = 0.2) +
#  scale_fill_continuous(type = "gradient", low="yellow", high="red") +
#  labs(fill = "log1p(Change in Population Density)")


########################### Combine with extinction data #################################


# combine species exintction and location data
#species_all <- merge(vertex, species_loc, by.x='Taxa', by.y='specsNm')

# check extinction dates match
#species_all$EX.Last.seen. == species_all$Extnc_Y

# set coordinate system
#species_map <- st_transform(species_loc, crs=mollweide_crs)

# fix projection error due to different versions GDAL (if using both local and hpc)
#st_crs(world_map) <- mollweide_crs

vertex_sf <- st_as_sf(vertex) # species_all_sf <- st_as_sf(species_all)
#st_crs(species_all_sf) <- mollweide_crs
#print("Joining data with map objects...")
# join with map
#grid <- st_join(x=world_map,
                #y=st_cast(species_all_sf) %>% dplyr::select(Taxa, Class.x, Class_Num, EX.Last.seen., Cont.2..Island.1.),
                #left=TRUE, join = st_intersects)

# save combined object
#saveRDS(grid, paste0(path_data, "species_moll_grid_all.rds"))

# put through time periods blocking
#time_periods <- c(1700, 1800, 1900, 2000)
# do after st_join so st_join only has to be run once
# group by specified time periods
#grid$Year_Block_Var <- NA
#for (t in 2:length(time_periods)-1) {
#  for (spec in 1:nrow(grid)) {
#    if (!is.na(grid$EX.Last.seen.[spec])) {
#      if (time_periods[t] < grid$EX.Last.seen.[spec] && grid$EX.Last.seen.[spec] <= time_periods[t+1]) {
#        grid$Year_Block_Var[spec] <- time_periods[t+1]  
#      }
#    }
#  }
#}



# subset by time period
vertex_sf <- vertex_sf[which(vertex_sf$Year_Block_Var == time),]

# sum by number of species
#grid_count_sub <- grid_sub %>%
#  group_by(gridId) %>%
#  summarise(SpeciesRichness = n_distinct(Taxa, na.rm = TRUE),
#            spsList = paste(Taxa, collapse = ';'),
#            .groups = 'drop')

# save object
# set name for time periods for output names
#time_block <- "1900"
#saveRDS(grid_count_sub, paste0(path_data, time_block, "_grid_spec_count.rds"))

# plot data
#ggplot() + 
#  geom_sf(data = grid_count_sub %>%  
#            mutate(SpeciesRichness=ifelse(SpeciesRichness==0, NA, SpeciesRichness)),
#          aes(fill=SpeciesRichness), col=NA, size=0.2) +
#  scale_fill_continuous(type = "gradient", low="yellow", high="red") +
#  geom_sf(data = world_plot, color='black', fill=NA, size=0.2) +
#  labs(fill='Species Richness')

#ggplot(data=popden_world_time_sf) + 
#  geom_sf(aes(fill=log1p(Prop_change)), size = 0.2) +
#  scale_fill_continuous(type = "gradient", low="yellow", high="red")

vertex_cent <- st_centroid(vertex_sf) #cent <- st_centroid(grid_count_sub)
#cent_sub <- cent[which(cent$SpeciesRichness != 0),]

#vertex_world <- as.data.frame(st_coordinates(st_centroid(grid_richness)))
#names(vertex_world) <- c("Longitude", "Latitude")

#ggplot() + 
#geom_sf(data = cent_sub, aes(col=SpeciesRichness), size=1) +
#scale_color_gradient(low="green", high="black") +
#geom_sf(data = world_plot, color='black', fill=NA, size=0.2) +
#labs(col ='Species Richness')


ggplot() + 
  geom_sf(data=popden_world_time_sf, aes(fill=log1p(Prop_change)), size = 0.2) +
  scale_fill_continuous(type = "gradient", low="yellow", high="red") +
  geom_sf(data = world, color='grey', fill=NA, size=0.2) +
  geom_sf(data = vertex_cent, size=1) +
  #scale_color_gradient(low="green", high="black") +
  labs(fill = "Population Density Change (log)")

#ggplot() + 
#  geom_sf(data=popden_world_time_sf, aes(fill=log1p(Prop_change)), size = 0.2) +
#  scale_fill_continuous(type = "gradient", low="yellow", high="red") +
#  geom_sf(data = cent_sub, aes(col=SpeciesRichness), size=1) +
#  scale_color_gradient(low="green", high="black") +
#  geom_sf(data = world_map, color='black', fill=NA, size=0.2) +
#  labs(col ='Species Richness', fill = "Population Density Change")









# run Pearson's spatial correlation to see if number of extinctions and change in human population
# density is correlated spatially

cent_popden <- st_centroid(popden_world_time_sf)

popden_world_xy <- st_coordinates(st_centroid(popden_world_time_sf))

library("SpatialPack")
correlation <- modified.ttest(popden_world$NoExSpec, popden_world$PopDenChange, popden_world_xy, nclass = 10)


## end of script
