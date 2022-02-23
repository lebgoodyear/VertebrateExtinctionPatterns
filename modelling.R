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
time_name <- "1800-2000_5y"

# set path to data and scripts
path_data <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/", time_name, "/")
path_rawdata <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/raw_data/")
path_out <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/outputs/", time_name, "/")

# load data
expoptot <- readRDS(paste0(path_data, time_name, "_expoptot.rds"))


#####################################################################################
#################################### Models #########################################


################################ Population density #################################


# we know mean is not equal to variance so quasipoisson or
# negative binomial distributions must be used

m1 <- glm(NoExSpec ~ Pop, family=quasipoisson(link="log"), data=expoptot)
summary(m1)

m2 <- glm.nb(NoExSpec ~ Pop, data=expoptot)
summary(m2)

par(mfrow=c(1,1))
plot(expoptot$Pop,expoptot$NoExSpec,pch=16)
lines(expoptot$Pop,predict(m1, type="response"))
lines(expoptot$Pop,predict(m2, type = "response"),lty=1,col="red")

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


#m3 <- glm(NoExSpec ~ PopDenChange, family=quasipoisson(link='log'), data=expop)
#summary(m3)

#m4 <- glm.nb(NoExSpec ~ PopDenChange, data=expop)
#summary(m4)

#par(mfrow=c(1,1))
#plot(expop$PopDenChange,expop$NoExSpec,pch=16)
#lines(expop$PopDenChange,predict(m3, type="response"))
#lines(expop$PopDenChange,predict(m4, type = "response"),lty=1,col="red")

#plot(m3)


########################################################################################
############################ Plots with regression lines ###############################


# plot human population density change (per time interval as specified above)
# vs number of extinctions, including regression line calculated by Poisson family GLM
# using log link function
#ggplot(data = expop, aes(x = PopDenChange, y = NoExSpec)) +
#  geom_point() +
#  geom_smooth(method = "glm", method.args = list(family = "quasipoisson")) +
#  theme_bw()

# plot with predictions
#ggplot(data = expop, aes(x = PopDenChange, y = NoExSpec)) +
#  geom_point() +
#  xlim(0, 2.0) +
#  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), fullrange=TRUE) +
#  theme_bw()

# plot human population density (per time interval as specified above)
# vs number of extinctions, including regression line calculated by Poisson family GLM
# using log link function
ggplot(data = expoptot, aes(x = Pop, y = NoExSpec)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson")) +
  theme_bw()

# plot with predictions
ggplot(data = expoptot, aes(x = Pop, y = NoExSpec)) +
  geom_point() +
  xlim(0, 13e9) +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), fullrange=TRUE) +
  theme_bw()

# plot with predictions
#ggplot(data = expoptot, aes(x = Pop, y = NoExSpec)) +
#  geom_point() +
#  xlim(0, 13e9) +
#  geom_smooth(method = "glm.nb", fullrange=TRUE) +
#  theme_bw()

# plot predictions for different scenarios
# using data from UN
# create separate plots for no extinctions and another for pop den data
# plot on same grid to directly compare


###########################################################################


un_pred <- read.csv(paste0(path_rawdata, "un_pop_pred_world.csv"))

un_pred <- as.data.frame(un_pred[,c("X2050", "X2100")])
un_2050 <- as.data.frame(un_pred$X2050*1000)
names(un_2050) <- "Pop"
un_2100 <- as.data.frame(un_pred$X2100*1000)
names(un_2100) <- "Pop"

pred_2050 <- predict(m1, un_2050, type = "response")
pred_2100 <- predict(m1, un_2100, type = "response")

# finding the confidence intervals
ci_2050 <- predict(m1, un_2050, se.fit=TRUE, type = "link")
ci_2100 <- predict(m1, un_2100, se.fit=TRUE, type = "link")

exp(ci_2050$fit - 2*ci_2050$se.fit)
exp(ci_2050$fit + 2*ci_2050$se.fit)

exp(ci_2100$fit - 2*ci_2100$se.fit)
exp(ci_2100$fit + 2*ci_2100$se.fit)

######## use negative binomial and quasipoisson since mean is not 
# equal to variance (overdispersion)
# this is possibly due to clustering in space

# also note that rate at which events occur is not constant (increases
# over time) - note sure how to account for this

# make estimate of lambda (mean) for poisson/quasi-poiisson/nb
lambda <- sum(vertex_tot$No_Ex_Spec_tot/length(time_periods))

# run a model for years since 1700 or something (e.g. 1700->0, 1705->1, 1710->2 etc.)

# years are unlikely to be correlated with each other 
# so no need for time series analysis? Number of extinctions does not depend
# on previous number of extinctions

# could look at rate of extinctions per year by dividing time periods accordingly...


## end of script
