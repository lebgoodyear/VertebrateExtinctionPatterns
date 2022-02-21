###########################################################################################
################################## Explore Data  ##########################################
###########################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Created: Feb 2022
# Last updated: Feb 2022

# clear workspace
rm(list=ls())


########################################################################################
################################### Set up #############################################


# load required packages
library("dplyr")
library("ggplot2")
library("gridExtra")
#library("sf")
#library("rnaturalearth")
#install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
#library("rnaturalearthhires")

# select time name to analyse data from set time periods
time_name <- "1460-2000_20y"

# set path to data and scripts
path_data <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/", time_name, "/")
path_out <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/outputs/", time_name, "/")

# load data
vertex <- readRDS(paste0(path_data, time_name, "_vertex.rds"))
vertex_class <- readRDS(paste0(path_data, time_name, "_vertex_class.rds"))
vertex_tot <- readRDS(paste0(path_data, time_name, "_vertex_tot.rds"))
vertex_tot_noamph <- readRDS(paste0(path_data, time_name, "_noamph_vertex_tot.rds"))

popdenchange <- readRDS(paste0(path_data, time_name, "_popden_change_tot.rds"))

expopden <- readRDS(paste0(path_data, time_name, "_expopden.rds"))
expop <- readRDS(paste0(path_data, time_name, "_expop.rds"))


##################################################################################
######################### Extinction data only ###################################


hist(vertex_tot$No_Ex_Spec_tot)
hist(vertex_tot_noamph$No_Ex_Spec_tot)

mean(vertex_tot$No_Ex_Spec_tot)
var(vertex_tot$No_Ex_Spec_tot)

mean(vertex_tot_noamph$No_Ex_Spec_tot)
var(vertex_tot_noamph$No_Ex_Spec_tot)

# note that mean is not equal to variance so we have overdispersion
# and therefore we need to use quasipoisson or negative binomial model


############################# Extinctions over time ###############################


# histogram of extinctions by time
hist(vertex$EX.Last.seen.)
vertex_noamph <- vertex[which(!(vertex$Class.x == "Amphibia")),]
hist(vertex_noamph$EX.Last.seen.)

# plot extinctions over time intervals as bar and scatter plots
barplot(vertex_tot$No_Ex_Spec_tot, vertex_tot$Year_Block_Var)
plot(vertex_tot$Year_Block_Var, vertex_tot$No_Ex_Spec_tot)

# plot extinctions with amphibians
barplot(vertex_tot_noamph$No_Ex_Spec_tot, vertex_tot_noamph$Year_Block_Var)
plot(vertex_tot_noamph$Year_Block_Var, vertex_tot_noamph$No_Ex_Spec_tot)


#########################################################################################
#################################### Initial plots ######################################


######################## Separate by year block #########################################


# create plotting object for pop density data
plot1 <- ggplot() +
  geom_point(data = popdenchange, aes(x=Year, y=Prop_change_tot)) +
  theme_bw()
# create plotting object for total extinctions, split by vertebrate class
plot2 <- ggplot() +
  geom_point(data = vertex_class, aes(x=Year_Block_Var, y=No_Ex_Spec, col=Class.x)) +
  theme_bw()
# create plotting object for total extinctions across all classes
plot3 <- ggplot() +
  geom_point(data = vertex_tot, aes(x=Year_Block_Var, y=No_Ex_Spec_tot)) +
  theme_bw()
# create plotting object for total extinctions without amphibians
plot4 <- ggplot() +
  geom_point(data = vertex_tot_noamph, aes(x=Year_Block_Var, y=No_Ex_Spec_tot)) +
  theme_bw()

# plot above plot on a single grid
grid.arrange(plot1, plot2, plot3, plot4, ncol=2) 


ggplot(data = expopden, aes(x=Year, y=NoExSpec)) +
  geom_bar(stat='identity') +
  theme_bw()

ggplot(data = expop, aes(x=Year, y=PopDenChange)) +
  geom_bar(stat='identity') +
  theme_bw()

ggplot() +
  geom_point(data = expopden, aes(x=Year, y=PopDen)) +
  theme_bw()


#################################### Combined #########################################


ggplot() +
  geom_point(data = expopden, aes(x=PopDen, y=NoExSpec)) +
  theme_bw()

ggplot() +
  geom_point(data = expop, aes(x=PopDenChange, y=NoExSpec)) +
  theme_bw()


################################## Correlations ########################################


# calculate correlation between number of extinctions and change in population density
cor(expop$NoExSpec, expop$PopDenChange)
# plot correlations
pairs(expop)

# calculate correlation between number of extinctions and population density
cor(expopden$NoExSpec, expopden$PopDen)
# plot correlations
pairs(expopden)


## end of script
