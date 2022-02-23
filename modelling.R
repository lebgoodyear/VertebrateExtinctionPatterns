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
time_name <- "1700-2000_10y"

# set path to data and scripts
path_data <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/", time_name, "/")
path_rawdata <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/raw_data/")
path_out <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/outputs/", time_name, "/")

# load data
expoptot <- readRDS(paste0(path_data, time_name, "_expoptot.rds"))


#####################################################################################
#################################### Models #########################################


#################################### Population ######################################


# we know mean is not equal to variance so quasipoisson or
# negative binomial distributions must be used

# quasipoisson
m1 <- glm(NoExSpec ~ Pop, family=quasipoisson(link="log"), data=expoptot)
summary(m1)

# negative binomial
m2 <- glm.nb(NoExSpec ~ Pop, data=expoptot)
summary(m2)

# plot two models over data
par(mfrow=c(1,1))
plot(expoptot$Pop,expoptot$NoExSpec,pch=16)
lines(expoptot$Pop,predict(m1, type="response"))
lines(expoptot$Pop,predict(m2, type = "response"),lty=1,col="red")

# negative binomial model will produce much larger predictions
# than quasipoisson since curve has a larger growth rate

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

# plot human population (per time interval as specified above)
# vs number of extinctions, including regression line calculated by quasioisson family GLM
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

# plot with predictions for negative binomial
#ggplot(data = expoptot, aes(x = Pop, y = NoExSpec)) +
#  geom_point() +
#  xlim(0, 13e9) +
#  geom_smooth(method = "glm.nb", fullrange=TRUE) +
#  theme_bw()
# as suspected, these estimates are much larger

# plot predictions for different scenarios
# using data from UN
# create separate plots for no extinctions and another for pop den data
# plot on same grid to directly compare


########################################################################################
################################## Predictions #########################################

pop_full <- read.csv(paste0(path_rawdata, "population_past_future.csv"), stringsAsFactors = F)
names(pop_full)[4] <- "Pop"
#mil2 <- pop_full[which(pop_full$Year %in% c(2005, 2010, 2015)),]
mil2 <- pop_full[which(pop_full$Year == 2010),]
mil2 <- as.data.frame(mil2[which(mil2$Entity=="World"),])
mil2 <- subset(mil2, select=c("Year", "Pop"))

pred_mil2 <- predict(m1, mil2, type = "response")
pred_mil2_df <- as.data.frame(cbind(mil2$Year, pred_mil2, mil2$Pop))
names(pred_mil2_df) <- c("Year", "NoExSpec", "Pop")
 

 
un_pred_full <- read.csv(paste0(path_rawdata, "un_pop_pred_world.csv"))
un_pred_full[,2:ncol(un_pred_full)] <- un_pred_full[,2:ncol(un_pred_full)]*1000

#un_pred <- as.data.frame(un_pred_full[,c("X2050", "X2100")])
#names(un_2050) <- "Pop"
#names(un_2100) <- "Pop"

#pred_2050 <- predict(m1, un_2050, type = "response")
#pred_2100 <- predict(m1, un_2100, type = "response")

# finding the confidence intervals
#ci_2050 <- predict(m1, un_2050, se.fit=TRUE, type = "link")
#ci_2100 <- predict(m1, un_2100, se.fit=TRUE, type = "link")

#exp(ci_2050$fit - 2*ci_2050$se.fit)
#exp(ci_2050$fit + 2*ci_2050$se.fit)

#exp(ci_2100$fit - 2*ci_2100$se.fit)
#exp(ci_2100$fit + 2*ci_2100$se.fit)


######## use negative binomial and quasipoisson since mean is not 
# equal to variance (overdispersion)
# this is possibly due to clustering in space

# also note that rate at which events occur is not constant (increases
# over time) - note sure how to account for this

# make estimate of lambda (mean) for poisson/quasi-poisson/nb
#lambda <- sum(expoptot$NoExSpec/nrow(expoptot))

# run a model for years since 1700 or something (e.g. 1700->0, 1705->1, 1710->2 etc.)

# years are unlikely to be correlated with each other 
# so no need for time series analysis? Number of extinctions does not depend
# on previous number of extinctions

# could look at rate of extinctions per year by dividing time periods accordingly...


#########################################################################################
########################### Plots for different scenarios ###############################


years <- colnames(un_pred_full[2:ncol(un_pred_full)])
years <- gsub("X", "", years)
years <- as.numeric(years)

years <- seq(2020, 2100, 10)
un10 <- un_pred_full[,which(names(un_pred_full) %in% c("X2020", "X2030", "X2040", "X2050", "X2060", "X2070", "X2080", "X2090", "X2100"))]

create_pred_df <- function(i) {
  sub <- as.data.frame(t(un10[i,1:ncol(un10)]))
  names(sub) <- "Pop"
  pred_sub <- predict(m1, sub, type = "response")
  predt_sub <- cbind(years, pred_sub, sub)
  names(predt_sub) <- c("Year", "NoExSpec", "Pop")
  
  return(predt_sub)
}

l95 <- create_pred_df(1)
med <- create_pred_df(3)
u95 <- create_pred_df(5)

expop_l95 <- rbind(expoptot, pred_mil2_df, l95)
expop_med <- rbind(expoptot, pred_mil2_df, med)
expop_u95 <- rbind(expoptot, pred_mil2_df, u95)

# plot with predictions
ggplot() +
  geom_point(data = expop_l95, aes(x = Year, y = NoExSpec), col="blue") +
  geom_point(data = expop_med, aes(x = Year, y = NoExSpec), col="red") +
  geom_point(data = expop_u95, aes(x = Year, y = NoExSpec), col="black") +
  geom_smooth(mapping = aes(x = Year, y = NoExSpec), data = expop_l95, 
              method = "gam", method.args = list(family = "quasipoisson"), fullrange=TRUE, col="blue") +
  geom_smooth(mapping = aes(x = Year, y = NoExSpec), data = expop_med, 
              method = "gam", method.args = list(family = "quasipoisson"), fullrange=TRUE, col="red") +
  geom_smooth(mapping = aes(x = Year, y = NoExSpec), data = expop_u95, 
              method = "gam", method.args = list(family = "quasipoisson"), fullrange=TRUE, col="black") +
  theme_bw()


## end of script
