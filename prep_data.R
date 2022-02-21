###########################################################################################
################################# Prepare data  ###########################################
###########################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Created: Jan 2022
# Last updated: Feb 2022

# clear workspace
rm(list=ls())


########################################################################################
################################### Set up #############################################


# load required packages
library("dplyr")
library("ggplot2")
library("gridExtra")
library("sf")
library("rnaturalearth")
#install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
library("rnaturalearthhires")

# set time periods
# year will be final year of interval
# e.g. c(1700, 1750, 1800) will be 1701-1750 and 1751-1800
time_periods <- c(1700, 1800, 1900, 2000)
#time_periods <- seq(1950, 1995, 5)
#time_periods <- seq(1810, 1995, 10)

# set name for time periods for output names
time_name <- "1700-2000_100y"

# set path to data and scripts
path_data <- "~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/raw_data/"
path_out_temp <- "~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/"
path_loc <- "~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/species_locations/"

# create directory to store the data for this time set
# first check if directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(path_out_temp, time_name))), dir.create(file.path(paste0(path_out_temp, time_name))), FALSE)
# set directory for results to be sent to
path_out <- paste0(path_out_temp, time_name, "/")


########################################################################################
################################ Extinction data prep ##################################


# read in vertebrate extinction data
vertex_full <- read.csv(paste0(path_data, "vert_extinctions.csv"))

# combine species and genus column to create column to use with location data
vertex_full <- vertex_full %>% unite('Taxa', Genus:Species, sep = " ", remove = F)

# remove uncertain extinction time intervals
vertex <- vertex_full[which(!(vertex_full$EX.Last.seen. == "1700-1750")),]
vertex <- vertex[which(!(vertex$EX.Last.seen == "1500-1600")),]
vertex <- vertex[!grepl("s", vertex$EX.Last.seen),]
vertex <- vertex[!grepl("century", vertex$EX.Last.seen),]
vertex <- vertex[!grepl("before", vertex$EX.Last.seen),]
vertex <- vertex[!grepl("approx", vertex$EX.Last.seen),]

# set as date as numeric
vertex$EX.Last.seen. <- as.numeric(as.character(vertex$EX.Last.seen.))

# species locations
species_loc <- read_sf(paste0(path_loc, "database_SPECIES.shp"))

# combine species exintction and location data
vertex <- merge(vertex, species_loc, by.x='Taxa', by.y='specsNm')

# group by specified time periods
vertex$Year_Block_Var <- NA
for (t in 2:length(time_periods)-1) {
  for (spec in 1:nrow(vertex)) {
    if (time_periods[t] < vertex$EX.Last.seen.[spec] && vertex$EX.Last.seen.[spec] <= time_periods[t+1]) {
      vertex$Year_Block_Var[spec] <- time_periods[t+1]  
    }
  }
}

# remove NAs, which corresponds to any:
# species extinct after 2000
# species without exact years of exinction
# species that went extinct before first time in time periods
vertex <- vertex[which(!is.na(vertex$Year_Block_Var)),]

# save for mapping
saveRDS(vertex, paste0(path_out, time_name, "_vertex.rds"))

# compute different summaries
# sum over required columns
vertex_block <- vertex %>% group_by(Class.x, Year_Block_Var, Class_Num) %>% summarize(No_Ex_Spec = n())
# include island/continent in summary
vertex_block_contis <- vertex %>% group_by(Class.x, Cont.2..Island.1., Year_Block_Var, Class_Num) %>% summarize(No_Ex_Spec = n())
# total extinctions per year
vertex_tot <- vertex_block %>% group_by(Year_Block_Var) %>% summarize(No_Ex_Spec_tot = sum(No_Ex_Spec))

# remove amphibians
vertex_noamph <- vertex_block[which(!(vertex_block$Class == "Amphibia")),]
# total extinctions without amphibians
vertex_noamph_tot <- vertex_noamph %>% group_by(Year_Block_Var) %>% summarize(No_Ex_Spec_tot = sum(No_Ex_Spec))

# save results to import into later scripts
saveRDS(vertex_block, paste0(path_out, time_name, "_vertex_class.rds"))
saveRDS(vertex_block_contis, paste0(path_out, time_name, "_vertex_contis.rds"))
saveRDS(vertex_tot, paste0(path_out, time_name, "_vertex_tot.rds"))
saveRDS(vertex_noamph, paste0(path_out, time_name, "_noamph_vertex_class.rds"))
saveRDS(vertex_noamph_tot, paste0(path_out, time_name, "_noamph_vertex_tot.rds"))


#########################################################################################
############################# Population density data prep ##############################


# read in population density data
popden_full <- read.csv(paste0(path_data, "human_pop_den.csv"), stringsAsFactors = F)

# rename population density column
names(popden_full)[4] <- "PopDen"

# filter by time periods
popden <- popden_full[which(popden_full$Year %in% time_periods),]
popden <- popden[which(!(popden$PopDen == 0)),]


####### Summarise over mean population density across countries per time period ########


# load world map data
world <- ne_countries(scale = 'large', returnclass = 'sf')

# find area of countries
world$area_km2 <- st_area(world$geometry)/1e6
world_sub <- as.data.frame(cbind(world$name, world$area_km2))
names(world_sub) <- c("name", "area_km2")
world_sub$area_km2 <- as.numeric(as.character(world_sub$area_km2))

# fix major naming discrepancies
popden$Entity[which(popden$Entity == "United States")] <- "United States of America"
popden$Entity[which(popden$Entity == "Cote d'Ivoire")] <- "Côte d'Ivoire"
popden$Entity[which(popden$Entity == "Democratic Republic of Congo")] <- "Dem. Rep. Congo"
popden$Entity[which(popden$Entity == "Equatorial Guinea")] <- "Eq. Guinea"
popden$Entity[which(popden$Entity == "Central African Republic")] <- "Central African Rep."

# merge area data and pop den data
popden_world_area <- merge(popden, world_sub, by.x="Entity", by.y="name")

# calculate mean population density across countries per time period,
# accounting for area weighting
popden_tot <- popden_world_area %>% group_by(Year) %>% summarize(Popden_tot = sum(PopDen * area_km2)/sum(area_km2))

# save result for use in later scripts
saveRDS(popden_tot, paste0(path_out, time_name, "_popden_tot.rds"))


######### Summarise over mean change in population density across time periods ##########


# calculate change in population density between time periods
popden$Prop_change <- NA
popden_final <- data.frame()
for (area in unique(popden$Entity)) {
  popden_sub <- popden[which(popden$Entity == area),]
  if (nrow(popden_sub) > 1) {
    for (time in (2:nrow(popden_sub))) {
      popden_sub$Prop_change[time] <- (popden_sub$PopDen[time] - popden_sub$PopDen[time-1])/popden_sub$PopDen[time-1]
    }
  } else {popden_sub$Prop_change[1] <- NA}
  popden_final <- rbind(popden_final, popden_sub)
}

# remove any NA proportional increases
popden_final <- popden_final[which(!(is.na(popden_final$Prop_change))),]

# save this result for mapping
saveRDS(popden_final, paste0(path_out, time_name, "_popden_by_country.rds"))

# summarise by the mean proportional change across countries per year
# note this accounts for population density weighting
prop_change_tot <- popden_final %>% group_by(Year) %>% summarize(Prop_change_tot = sum(Prop_change * PopDen)/sum(PopDen))

# save result for use in later scripts
saveRDS(prop_change_tot, paste0(path_out, time_name, "_popden_change_tot.rds"))

## manually check summarise results for mean proportional change across countries per year

# # create temporary pop den dataframe
# popden_final_temp <- popden_final
# # set temporary column for weighted numerator in mean calculation
# popden_final_temp$Prop_change_temp <- NA
# # set column to sum pop densities over each time period
# popden_final_temp$time_sum <- NA
# # calculate proportional change multiplied by pop den for use as weighted numerator
# for (row in 1:nrow(popden_final_temp)) {
#   popden_final_temp$Prop_change_temp[row] <- popden_final_temp$Prop_change[row] * popden_final_temp$PopDen[row]
# }
# # sum population densities over each time period but repeatedly store per relevant row for later calculation  
# for (time in unique(popden_final_temp$Year)) {
#   popden_final_temp$time_sum[which(popden_final_temp$Year == time)] <- sum(popden_final_temp$PopDen[which(popden_final_temp$Year == time)])
# }
# # sum across weighted numerator
# prop_temp_change_global <- popden_final_temp %>% group_by(Year, time_sum) %>% summarize(Prop_change_tot = sum(Prop_change_temp))
# # set final proprtional change column
# prop_temp_change_global$prop_final <- NA
# # divide weighted numerator by sum of pop densities per year
# for (row in (1:nrow(prop_temp_change_global))) {
#   prop_temp_change_global$prop_final[row] <- prop_temp_change_global$Prop_change_tot[row]/prop_temp_change_global$time_sum[row]
# }
# # view results to compare to prop_change_tot
# prop_temp_change_global
# prop_change_tot


#########################################################################################
###################### Combine vertex and population datasets ###########################


# combine total extinctions and total pop den changes into one dataframe
expop <- merge(vertex_tot, prop_change_tot, by.x = "Year_Block_Var", by.y = "Year")
names(expop) <- c("Year", "NoExSpec", "PopDenChange")

# combine total extinctions w/o amphibians and total pop den changes into one dataframe
expop_noamph <- merge(vertex_noamph_tot, prop_change_tot, by.x = "Year_Block_Var", by.y = "Year")
names(expop_noamph) <- c("Year", "NoExSpec", "PopDenChange")

# combine total extinctions and total pop densities into one dataframe
expopden <- merge(vertex_tot, popden_tot, by.x = "Year_Block_Var", by.y = "Year")
names(expopden) <- c("Year", "NoExSpec", "PopDen")

# combine total extinctions w/o amphibians and total pop den changes into one dataframe
expopden_noamph <- merge(vertex_noamph_tot, popden_tot, by.x = "Year_Block_Var", by.y = "Year")
names(expopden_noamph) <- c("Year", "NoExSpec", "PopDen")

# save results to import into later scripts
saveRDS(expop, paste0(path_out, time_name, "_expop.rds"))
saveRDS(expop_noamph, paste0(path_out, time_name, "_expop_noamph.rds"))
saveRDS(expopden, paste0(path_out, time_name, "_expopden.rds"))
saveRDS(expopden_noamph, paste0(path_out, time_name, "_expopden_noamph.rds"))


## end of script
