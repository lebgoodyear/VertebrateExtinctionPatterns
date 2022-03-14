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
time_name <- "1805-2000_5y"

# set path to data and scripts
path_data <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/", time_name, "/")
path_out <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/outputs/", time_name, "/")

# load data
# extinction data
vertex <- readRDS(paste0(path_data, time_name, "_vertex.rds"))
vertex_class <- readRDS(paste0(path_data, time_name, "_vertex_class.rds"))
vertex_tot <- readRDS(paste0(path_data, time_name, "_vertex_tot.rds"))
vertex_tot_noamph <- readRDS(paste0(path_data, time_name, "_noamph_vertex_tot.rds"))

# population data
poptot <- readRDS(paste0(path_data, time_name, "_pop_tot.rds"))

# extinction/population data
expop <- readRDS(paste0(path_data, time_name, "_expoptot.rds"))


##################################################################################
######################### Extinction data only ###################################


# view distribution of number of extinctions
hist(vertex_tot$No_Ex_Spec_tot)
hist(vertex_tot_noamph$No_Ex_Spec_tot)

# we can see a rough Poisson distribution pattern here

# view mean and variance
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

# plot extinctions without amphibians
barplot(vertex_tot_noamph$No_Ex_Spec_tot, vertex_tot_noamph$Year_Block_Var)
plot(vertex_tot_noamph$Year_Block_Var, vertex_tot_noamph$No_Ex_Spec_tot)


#########################################################################################
#################################### Initial plots ######################################


######################## Separate by year block #########################################


# create plotting object for pop data
plot1 <- ggplot() +
  geom_point(data = poptot, aes(x=Year, y=Population)) +
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

# plot as bar chart
ggplot(data = expop, aes(x=Year, y=NoExSpec)) +
  geom_bar(stat='identity') +
  theme_bw()

# plot population
ggplot() +
  geom_point(data = expop, aes(x=Year, y=Pop)) +
  theme_bw()


#################################### Combined #########################################


# plot population vs number of extinctions
ggplot() +
  geom_point(data = expop, aes(x=Pop, y=NoExSpec)) +
  theme_bw()


################################## Correlations ########################################


# calculate correlation between number of extinctions and change in population density
corr <- cor.test(expop$NoExSpec, expop$Pop)
cat(capture.output(corr), file=paste0(path_out, "correlations.csv"))
# plot correlations
pairs(expop)


## end of script
