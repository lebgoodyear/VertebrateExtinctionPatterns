###########################################################################################
############################ Statistical analysis for #####################################
######################### Change in human population density ##############################
###################################### and ################################################
######################## Number of vertebrate extinctions #################################
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
library("MASS") # negative binomial GLM

# select time name to analyse data from set time periods
time_name <- "1460-2000_20y"

# set path to data and scripts
path_data <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/", time_name, "/")
path_out <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/outputs/", time_name, "/")

# load data
vertex <- readRDS(paste0(path_data, time_name, "_vertex.rds"))
vertex_tot <- readRDS(paste0(path_data, time_name, "_vertex_tot.rds"))
vertex_tot_noamph <- readRDS(paste0(path_data, time_name, "_noamph_vertex_tot.rds"))

expopden <- readRDS(paste0(path_data, time_name, "_expopden.rds"))
expop <- readRDS(paste0(path_data, time_name, "_expop.rds"))


#####################################################################################
#################################### Models #########################################


################################ Population density #################################


# we know mean is not equal to variance so quasipoisson or
# negative binomial distributions must be used

m1 <- glm(NoExSpec ~ PopDen, family=quasipoisson(link="log"), data=expopden)
summary(m1)

m2 <- glm.nb(NoExSpec ~ PopDen, data=expopden)
summary(m2)

par(mfrow=c(1,1))
plot(expopden$PopDen,expopden$NoExSpec,pch=16)
lines(expopden$PopDen,predict(m1, type="response"))
lines(expopden$PopDen,predict(m2, type = "response"),lty=1,col="red")

#plot(m1)


############### overdispersion notes


# model using poisson
#m0 <- glm(NoExSpec ~ PopDen, family=poisson(link="log"), data=expopden)
#summary(m0)

# plot mean vs variance to view overdispersion
#plot(log(fitted(m1)),
#     log((expopden$NoExSpec-fitted(m1))^2),
#     xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),
#     pch=20,col="blue")

#abline(0,1) ## 'variance = mean' line

# calculate dispersion parameter
#dp = sum(residuals(m1,type ="pearson")^2)/m1$df.residual
#dp

# how much are coefficient estimates affected by the overdispersion?
# this is the same as using quasipoisson
#summary(m1,dispersion = dp)


######################### Change in Population density ##################################


m3 <- glm(NoExSpec ~ PopDenChange, family=quasipoisson(link='log'), data=expop)
summary(m3)

m4 <- glm.nb(NoExSpec ~ PopDenChange, data=expop)
summary(m4)

par(mfrow=c(1,1))
plot(expop$PopDenChange,expop$NoExSpec,pch=16)
lines(expop$PopDenChange,predict(m3, type="response"))
lines(expop$PopDenChange,predict(m4, type = "response"),lty=1,col="red")

#plot(m3)


########################################################################################
############################ Plots with regression lines ###############################


# plot human population density change (per time interval as specified above)
# vs number of extinctions, including regression line calculated by Poisson family GLM
# using log link function
ggplot(data = expop, aes(x = PopDenChange, y = NoExSpec)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson")) +
  theme_bw()

# plot with predictions
ggplot(data = expop, aes(x = PopDenChange, y = NoExSpec)) +
  geom_point() +
  xlim(0, 2.0) +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), fullrange=TRUE) +
  theme_bw()

# plot human population density (per time interval as specified above)
# vs number of extinctions, including regression line calculated by Poisson family GLM
# using log link function
ggplot(data = expopden, aes(x = PopDen, y = NoExSpec)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson")) +
  theme_bw()

# plot with predictions
ggplot(data = expopden, aes(x = PopDen, y = NoExSpec)) +
  geom_point() +
  xlim(0, 70) +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), fullrange=TRUE) +
  theme_bw()

# plot predictions for different scenarios
# using data from UN
# create separate plots for no extinctions and another for pop den data
# plot on same grid to directly compare


## end of script
