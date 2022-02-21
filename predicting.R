###########################################################################################
######################### Predictions based on UN projections  ############################
###########################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Created: Feb 2022
# Last updated: Feb 2022

# clear workspace
rm(list=ls())


########################################################################################
################################### Set up #############################################


# read in population density data
pop_full <- read.csv(paste0(path, "population.csv"), stringsAsFactors = F)

# rename population density column
names(pop_full)[4] <- "Pop"

# filter by time periods
pop <- pop_full[which(pop_full$Year %in% time_periods),]
pop <- pop[which(!(pop$Pop == 0)),]

# remove continent data
pop <- pop[which(!(pop$Code == "")),]
pop <- pop[which(!(pop$Entity == "World")),]


####### Summarise over mean population density across countries per time period ########


# calculate mean population density across countries per time period,
pop_tot <- pop %>% group_by(Year) %>% summarize(Pop_tot = sum(Pop))


######### Summarise over mean change in population density across time periods ##########


# calculate change in population density between time periods
pop$Prop_change <- NA
pop_final <- data.frame()
for (area in unique(pop$Entity)) {
  pop_sub <- pop[which(pop$Entity == area),]
  for (time in (2:nrow(pop_sub))) {
    pop_sub$Prop_change[time] <- (pop_sub$Pop[time] - pop_sub$Pop[time-1])/pop_sub$Pop[time-1]
  }
  pop_final <- rbind(pop_final, pop_sub)
}

# remove any NA proportional increases
pop_final <- pop_final[which(!(is.na(pop_final$Prop_change))),]

# summarise by the mean proportional change across countries per year
# note this accounts for population density weighting
pop_change_tot <- pop_final %>% group_by(Year) %>% summarize(Prop_change_tot = sum(Prop_change * Pop)/sum(Pop))


# combine total extinctions and total pop den changes into one dataframe
expop <- merge(vertex_tot, pop_change_tot, by.x = "Year_Block_Var", by.y = "Year")
names(expop) <- c("Year", "NoExSpec", "PopDenChange")

# combine total extinctions w/o amphibians and total pop den changes into one dataframe
expop_noamph <- merge(vertex_noamph_tot, pop_change_tot, by.x = "Year_Block_Var", by.y = "Year")
names(expop_noamph) <- c("Year", "NoExSpec", "PopDenChange")

# combine total extinctions and total pop densities into one dataframe
expopden <- merge(vertex_tot, pop_tot, by.x = "Year_Block_Var", by.y = "Year")
names(expopden) <- c("Year", "NoExSpec", "PopDen")

# combine total extinctions w/o amphibians and total pop den changes into one dataframe
expopden_noamph <- merge(vertex_noamph_tot, pop_tot, by.x = "Year_Block_Var", by.y = "Year")
names(expopden_noamph) <- c("Year", "NoExSpec", "PopDen")



################################################################################


# plot human population density (per time interval as specified above)
# vs number of extinctions, including regression line calculated by Poisson family GLM
# using log link function
ggplot(data = expopden, aes(x = PopDen, y = NoExSpec)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) +
  theme_bw()

# plot with predictions
ggplot(data = expopden, aes(x = PopDen, y = NoExSpec)) +
  geom_point() +
  xlim(0, 13e9) +
  geom_smooth(method = "glm", method.args = list(family = quasipoisson(link="log")), fullrange=TRUE) +
  geom_smooth(method = "glm", method.args = list(family = poisson(link="log")), fullrange=TRUE) +
  geom_smooth(method = "glm", method.args = list(family = quasipoisson(link="identity")), fullrange=TRUE) +
  theme_bw()

# plot with predictions
ggplot(data = expopden, aes(x = PopDen, y = NoExSpec)) +
  geom_point() +
  xlim(0, 13e9) +
  geom_smooth(method = "glm", method.args = list(family = quasipoisson(link="log")), fullrange=TRUE) +
  theme_bw()


###########################################################################


un_pred <- read.csv(paste0(path, "un_pop_pred_world.csv"))

un_pred <- as.data.frame(un_pred[,c("X2050", "X2100")])
un_2050 <- as.data.frame(un_pred$X2050*1000)
names(un_2050) <- "PopDen"
un_2100 <- as.data.frame(un_pred$X2100*1000)
names(un_2100) <- "PopDen"

pred_2050 <- predict(m1, un_2050, interval = "prediction")
pred_2100 <- predict(m1, un_2100, interval = "prediction")

pred_2100 <- predict(m4, un_2100, interval = "prediction")

# plot with predictions
ggplot(data = expopden, aes(x = PopDen, y = NoExSpec)) +
  geom_point() +
  xlim(0, 13e9) +
  geom_smooth(method = "glm.nb", fullrange=TRUE) +
  theme_bw()

# plot with predictions
ggplot(data = expop, aes(x = PopDenChange, y = NoExSpec)) +
  geom_point() +
  xlim(0, 3) +
  geom_smooth(method = "glm.nb", fullrange=TRUE) +
  theme_bw()



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




