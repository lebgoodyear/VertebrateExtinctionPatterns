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
library("car") # Anova() for a model
library("mgcv") # for gam
library("caret") # for cross validation

# select time name to analyse data from set time periods
time_name <- "1805-2000_5y"
time_periods <- seq(1805, 2000, 5)
year_split <- 5

# set path to data and scripts
path_data <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/", time_name, "/")
path_rawdata <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/raw_data/")
path_out_temp <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/outputs/")

# create directory to store the data for this time set
# first check if directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(path_out_temp, time_name))), dir.create(file.path(paste0(path_out_temp, time_name))), FALSE)
# set directory for results to be sent to
path_out <- paste0(path_out_temp, time_name, "/")

# load data
# extinction/population data
expoptot <- readRDS(paste0(path_data, time_name, "_expoptot.rds"))
# un population prediction data
un_pred_full <- read.csv(paste0(path_rawdata, "un_pop_pred_world.csv"))
# un data is in 1000s so need to multiply by 1000
un_pred_full[,2:ncol(un_pred_full)] <- un_pred_full[,2:ncol(un_pred_full)]*1000
# load population raw data to fill in missing years
pop_full <- read.csv(paste0(path_rawdata, "population_past_future.csv"), stringsAsFactors = F)
names(pop_full)[4] <- "Pop"

# define 3-fold cross validation with 3 repeats
train_control <- trainControl(method="repeatedcv", number=5, repeats=5)


#####################################################################################
#################################### Models #########################################


# from data exploration, we know data has rough Poisson distribution and that
# mean is not equal to variance so quasipoisson or negative binomial distributions must be used

# Models used:
# BASELINE and with CROSS-VALIDATION
# 1) m1 is a GLM using quasipoisson family with log link function
# 2) m2 is a GLM using quasipoisson family with log link function, with outliers removed
# 3) m3 is a negative binomial GLM
# 4) m4 is a negative binomial GLM, with outliers removed
# 5) m5 is a GAM using quasipoisson family with log link function
# 6) m6 is a GAM using quasipoisson family with log link function, with outliers removed
# 7) m7 is a negative binomial GAM
# 8) m8 is a negative binomial GAM, with outliers removed

# all these models are combined into one function with outputs stored for comparison


###########################################################################################
################################### Generate model function ###############################


# function for running and saving model results
run_models <- function(dataset, x, y, outliers=NA, output_name) {
  
  # INPUTS
  # dataset must be dataframe
  # x must be name of predictor column as string
  # y must be name of response column as string
  # outliers must be vector containing outlying years or NA
  # output_name is string to save outputs under
  
  # create dataset without outliers
  dataset_woo <- dataset[which(!(dataset$Year %in% outliers)),]
  
  # set formula for model creation
  form <- as.formula(paste(y, '~', x))
  
  # set max iterations for negative binomial models
  maxset = 1000
  
  
  # MODEL 1: quasipoisson GLM
  
  
  ## base
  m1 <- glm(form, family=quasipoisson(link="log"), data=dataset)
  r2_m1 <- with(summary(m1), 1 - deviance/null.deviance) # calculate r^2 value
  # save outputs
  cat("\n\n Quasipoisson GLM\n", 
      capture.output(summary(m1)), 
      file=paste0(path_out, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nR^2", capture.output(r2_m1), 
      file=paste0(path_out, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\n\nAnova \n", capture.output(Anova(m1)), 
      file=paste0(path_out, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  ## cross-validation
  m1a <- train(form,
              method = "glm", 
              family = "quasipoisson",
              data=dataset,
              trControl=train_control)
  r2_m1a <- with(summary(m1a), 1 - deviance/null.deviance) # calculate r^2 value
  # save outputs
  cat("\n\n Quasipoisson GLM with cross-validation\n", 
      capture.output(summary(m1a)), 
      file=paste0(path_out, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nR^2", capture.output(r2_m1a), 
      file=paste0(path_out, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nCV summary\n", capture.output(print(m1a)), 
      file=paste0(path_out, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  
  
  # MODEL 2: quasipoisson GLM without outliers
  
  
  if (!is.na(outliers)) {
    ## base
    m2 <- glm(form, family=quasipoisson(link="log"), data=dataset_woo)
    r2_m2 <- with(summary(m2), 1 - deviance/null.deviance) # calculate r^2 value
    # save outputs
    cat("\n\n Quasipoisson GLM without outliers\n", 
        capture.output(summary(m2)), 
        file=paste0(path_out, output_name, "_GLM.txt"), 
        sep = "\n", append=TRUE)
    cat("\n\nR^2", capture.output(r2_m2), 
        file=paste0(path_out, output_name, "_GLM.txt"), 
        sep = "\n", append=TRUE)
    cat("\n\n\nAnova \n", capture.output(Anova(m2)), 
        file=paste0(path_out, output_name, "_GLM.txt"), 
        sep = "\n", append=TRUE)
    ## cross-validation
    m2a <- train(form, 
                 method = "glm", 
                 family = "quasipoisson",
                 data=dataset_woo,
                 trControl=train_control)
    r2_m2a <- with(summary(m2a), 1 - deviance/null.deviance) # calculate r^2 value
    # save outputs
    cat("\n\n Quasipoisson GLM with cross-validation without outliers\n", 
        capture.output(summary(m2a)), 
        file=paste0(path_out, output_name, "_GLM.txt"), 
        sep = "\n", append=TRUE)
    cat("\n\nR^2", capture.output(r2_m2a), 
        file=paste0(path_out, output_name, "_GLM.txt"), 
        sep = "\n", append=TRUE)
  }
    
  # MODEL 3: negative binomial GLM
  
  
  ## base
  m3 <- glm.nb(form, data=dataset, maxit=maxset)
  r2_m3 <- with(summary(m3), 1 - deviance/null.deviance) # calculate r^2 value
  # save outputs
  cat("\n\n Negative Binomial GLM\n", 
      capture.output(summary(m3)), 
      file=paste0(path_out, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nR^2", capture.output(r2_m3), 
      file=paste0(path_out, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\n\nAnova \n", capture.output(Anova(m3)), 
      file=paste0(path_out, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  ## cross-validation
  m3a <- train(form, 
               method = "glm.nb", 
               maxit = maxset,
               data=dataset,
               trControl=train_control)
  r2_m3a <- with(summary(m1a), 1 - deviance/null.deviance) # calculate r^2 value
  # save outputs
  cat("\n\n Negative Binomial GLM with cross-validation\n", 
      capture.output(summary(m3a)), 
      file=paste0(path_out, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nR^2", capture.output(r2_m3a), 
      file=paste0(path_out, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  
  
  # MODEL 4: negative binomial GLM without outliers
  
  
  if (!is.na(outliers)) {
    ## base
    m4 <- glm.nb(form, data=dataset_woo, maxit=maxset)
    r2_m4 <- with(summary(m4), 1 - deviance/null.deviance) # calculate r^2 value
    # save outputs
    cat("\n\n Negative binomial GLM without outliers\n", 
        capture.output(summary(m4)), 
        file=paste0(path_out, output_name, "_GLMnb.txt"), 
        sep = "\n", append=TRUE)
    cat("\n\nR^2", capture.output(r2_m4), 
        file=paste0(path_out, output_name, "_GLMnb.txt"), 
        sep = "\n", append=TRUE)
    cat("\n\n\nAnova \n", capture.output(Anova(m4)), 
        file=paste0(path_out, output_name, "_GLMnb.txt"), 
        sep = "\n", append=TRUE)
    ## cross-validation
    m4a <- train(form, 
                 method = "glm.nb",
                 maxit=maxset,
                 data=dataset_woo,
                 trControl=train_control)
    r2_m4a <- with(summary(m1a), 1 - deviance/null.deviance) # calculate r^2 value
    # save outputs
    cat("\n\n Negative Binomial GLM with cross-validation without outliers\n", 
        capture.output(summary(m4a)), 
        file=paste0(path_out, output_name, "_GLMnb.txt"), 
        sep = "\n", append=TRUE)
    cat("\n\nR^2", capture.output(r2_m4a), 
        file=paste0(path_out, output_name, "_GLMnb.txt"), 
        sep = "\n", append=TRUE)
  }
  
  
  # MODEL 5: quasipoisson GAM
  
  
  ## base
  m5 <- gam(form, family=quasipoisson(link="log"), data=dataset)
  # save outputs
  cat("\n\n Quasipoisson GAM\n", 
      capture.output(summary(m5)), 
      file=paste0(path_out, output_name, "_GAM.txt"), 
      sep = "\n", append=TRUE)
  ## cross-validation
  m5a <- train(form, 
               method = "gam", 
               family = "quasipoisson",
               data=dataset,
               trControl=train_control)
  # save outputs
  cat("\n\n Quasipoisson GAM with cross-validation\n", 
      capture.output(summary(m5a)), 
      file=paste0(path_out, output_name, "_GAM.txt"), 
      sep = "\n", append=TRUE)

  
  # MODEL 6: quasipoisson GAM without outliers
  
  
  if (!is.na(outliers)) {
    ## base
    m6 <- gam(form, family=quasipoisson(link="log"), data=dataset_woo)
    # save outputs
    cat("\n\n Quasipoisson GAM without outliers\n", 
        capture.output(summary(m6)), 
        file=paste0(path_out, output_name, "_GAM.txt"), 
        sep = "\n", append=TRUE)
    ## cross-validation
    m6a <- train(form, 
                 method = "gam", 
                 family = "quasipoisson",
                 data=dataset_woo,
                 trControl=train_control)
    # save outputs
    cat("\n\n Quasipoisson GAM with cross-validation without outliers\n", 
        capture.output(summary(m6a)), 
        file=paste0(path_out, output_name, "_GAM.txt"), 
        sep = "\n", append=TRUE)
  }
    
    
  # MODEL 7: negative binomial GAM
  
  
  ## base
  m7 <- gam(form, family=nb(), data=dataset)
  # save outputs
  cat("\n\n Negative Binomial GAM\n", 
      capture.output(summary(m7)), 
      file=paste0(path_out, output_name, "_GAMnb.txt"), 
      sep = "\n", append=TRUE)
  ## cross-validation
  m7a <- train(form, 
               method = "gam", 
               family = nb(),
               data=dataset,
               trControl=train_control)
  # save outputs
  cat("\n\n Negative Binomial GAM with cross-validation\n", 
      capture.output(summary(m7a)), 
      file=paste0(path_out, output_name, "_GAMnb.txt"), 
      sep = "\n", append=TRUE)
  
  
  # MODEL 8: negative binomial GAM without outliers
  
  
  if (!is.na(outliers)) {
    ## base
    m8 <- gam(form, family=nb(), data=dataset_woo)
    # save outputs
    cat("\n\n Negative Binomial GAM without outliers\n", 
        capture.output(summary(m8)), 
        file=paste0(path_out, output_name, "_GAMnb.txt"), 
        sep = "\n", append=TRUE)
    ## cross-validation
    m8a <- train(form, 
                 method = "gam", 
                 family = nb(),
                 data=dataset_woo,
                 trControl=train_control)
    # save outputs
    cat("\n\n Negative Binomial GAM with cross-validation without outliers\n", 
        capture.output(summary(m8a)), 
        file=paste0(path_out, output_name, "_GAMnb.txt"), 
        sep = "\n", append=TRUE)
  }
  
  # save models as list
  if (!is.na(outliers)) {
    models <- list(m1,m1a,m2,m2a,m3,m3a,m4,m4a,m5,m5a,m6,m6a,m7,m7a,m8,m8a)
  } else {models <- list(m1,m1a,m3,m3a,m5,m5a,m7,m7a)}
  
  #return models
  return(models)
  
}


############### overdispersion notes


# model using poisson
#m0 <- glm(NoExSpec ~ Pop, family=poisson(link="log"), data=expoptot)
#summary(m0)

# plot mean vs variance to view overdispersion
#plot(log(fitted(m1)),
#     log((expoptot$NoExSpec-fitted(m1))^2),
#     xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),
#     pch=20,col="blue")

#abline(0,1) ## 'variance = mean' line

# calculate dispersion parameter
#dp = sum(residuals(m1,type ="pearson")^2)/m1$df.residual
#dp

# how much are coefficient estimates affected by the overdispersion?
# this is the same as using quasipoisson
#summary(m1,dispersion = dp)


#######################################################################################
################################### Run models ########################################


# function will be run over:
# a) TOTAL EXTINCTIONS, b) EACH CLASS SEPARATELY, AND c) ISLAND/CONTINENT
# 1) 1805-2000, split by 5 years
# 2) 1805-2000, split by 5 years, with 1900 and 1955 outliers removed
# 3) 1805-2000, split by 5 years, with 1900, 1955 outliers removed and 1995, 2000 removed to remove downward trend in GAMs
# 4) 1710-2000, split by 10 years
# 5) 1710-2000, split by 10 years, with 1900 outlier removed
# 6) 1710-2000, split by 10 yearsm with 1900 outlier removed and 2000 removed to remove downward trend in GAMs

# note outliers were chosen from the plots in explore_data.R
# note correlation for 1805-2000 (5 years) was 0.76 so quite high (explore_data.R)
# and correlation for 1710-2000 (10 years) was0.89 so very high (explore_data.R)

# set outliers depending on time period
if (year_split == 10) {
  outliersp <- c(1900)
  outlierspf <- c(1900, 2000)
} 
if (year_split == 5) {
  outliersp <- c(1900, 1955)
  outlierspf <- c(1900, 1955, 1995, 2000)
}


# a) TOTAL EXTINCTIONS

extot <- run_models(expoptot, "Pop", "NoExSpec", NA, "extot_na")
extot_outp <- run_models(expoptot, "Pop", "NoExSpec", outliersp, "extot_outp")
extot_outpf <- run_models(expoptot, "Pop", "NoExSpec", outlierspf, "extot_outpf")


# b) BY CLASS




########################################################################################
############################ Plots with regression lines ###############################


# plot human population (per time interval as specified above)
# vs number of extinctions, including regression line calculated by quasioisson family GLM
# using log link function

#ggplot(data = expoptot, aes(x = Pop, y = NoExSpec)) +
#  geom_point() +
#  geom_smooth(method = "glm", method.args = list(family = "quasipoisson")) +
#  theme_bw()

# plot with predictions
#ggplot(data = expoptot, aes(x = Pop, y = NoExSpec)) +
#  geom_point() +
#  xlim(0, 13e9) +
#  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), fullrange=TRUE) +
#  theme_bw()

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


##################################### Functions ########################################


# function to process data in order to make predictions
process_data <- function(dataset, year_split, scenario) {
  
  # INPUTS
  
  # dataset is dataframe of year ("Year") and predictor ("Pop")
  # model is the model you want to use to predict response
  # year_split is either 5 or 10, depending on time period you want to predict over
  # scenario is the prediction curve you want: "Upper 95",  "Lower 95" or "Median"
  
  # OUTPUTS
  
  # a dataframe containing population data for all past and future time points,
  # accroding to chosen UN scenario
  
  # only keep world pop data
  dataworld <- as.data.frame(dataset[which(dataset$Entity=="World"),])
  # subset by year and population
  dataworld <- subset(dataworld, select=c("Year", "Pop"))
  # subset by years in time_periods
  datasub <- as.data.frame(dataworld[which(dataworld$Year %in% time_periods),])
  
  # fill year gaps
  if (year_split == 10) {
    gap <- dataworld[which(dataworld$Year == 2010),]
  } 
  if (year_split == 5) {
    gap <- dataworld[which(dataworld$Year %in% c(2005, 2010, 2015)),]
  }
  
  # combine with dataset
  pops_past <- as.data.frame(rbind(datasub, gap))
  
  # store scenario bound names
  bounds_names <- un_pred_full$PredictionBound
  
  # split by required years
  years <- seq(2020, 2100, year_split)
  years_char <- as.character(years)
  years_char <- paste0("X", years_char) # set to match colnames
  
  # subset by years
  unpred <- un_pred_full[,which(names(un_pred_full) %in% years_char)]
  unpred <- as.data.frame(t(unpred)) # transpose
  names(unpred) <- bounds_names # set column names to scenario boundaries
  
  # add years to pops
  pops_future <- as.data.frame(cbind(years, unpred[[scenario]]))
  names(pops_future) <- c("Year", "Pop")
  
  # conbine past and future pops into one dataframe
  allpops <- rbind(pops_past, pops_future)
  
  return(allpops)
  
}
  
# function to make predictions on a per model basis
make_predictions <- function(dataset, model) {
  
  # extract populations
  pops <- as.data.frame(dataset$Pop)
  names(pops) <- "Pop"
  
  # make precitions
  if (class(model)[1] == "train") {
    preds <- predict(model, pops, type = "raw")
  } else {
    preds <- predict(model, pops, type = "response")
  }
  
  # combine pop dataset with predictions
  preddf <- cbind(dataset, preds)
  names(preddf) <- c("Year", "Pop", "NoExSpec")
  
  # finding the confidence intervals
  if (class(model)[1] == "train") {
    ci <- predict(model$finalModel, preddf, se.fit=TRUE, type = "link")
  } else {
    ci <- predict(model, preddf, se.fit=TRUE, type = "link")
  }
  
  # lower limit
  cil <- (ci$fit - 2*ci$se.fit)^2
  # upper limit
  ciu <- (ci$fit + 2*ci$se.fit)^2
  
  # combine all data and confidence intervals for plotting
  fulldf <- cbind(dataset, preds, cil, ciu)
  names(fulldf) <- c("Year", "Pop", "NoExSpec", "UpperCI", "LowerCI")
  
  return(fulldf)
  
}


########################################################################################
######################### Generate prediction datasets #################################


# create datasets for predictions
popsu95 <- process_data(pop_full, year_split, "Upper 95")
popsmed <- process_data(pop_full, year_split, "Median")
popsl95 <- process_data(pop_full, year_split, "Lower 95")

# create predictions for each model
dfs <- list(popsu95, popsmed, popsl95)

# list of models to predict
model_ls <- extot

# make predictions for each model for each dataset, save predictions as csv
# and plot each model on one plot for all three datasets
for (mod in 1:length(model_ls)) {
  pred <- vector(mode="list", length=length(dfs))
  for (df in 1:length(dfs)) {
    full_preds <- make_predictions(dfs[[df]], m3a)#model_ls[[mod]])
    write.csv(full_preds, paste0(path_out, "extot_m3a_", df, "_preds.csv"))
    pred[[df]] <- full_preds
  }
  pdf(file=paste0(path_out, "m3a_extot_preds.pdf"))
  print(ggplot() +
          geom_point(data = expoptot, aes(x = Year, y = NoExSpec), col="black") +
          geom_line(data = pred[[1]], aes(x = Year, y = NoExSpec), col="red") +
          geom_ribbon(data = pred[[1]], aes(x = Year, ymin=LowerCI, ymax=UpperCI), alpha=0.2, fill="red") +
          geom_line(data = pred[[2]], aes(x = Year, y = NoExSpec), col="blue") +
          geom_ribbon(data = pred[[2]], aes(x = Year, ymin=LowerCI, ymax=UpperCI), alpha=0.2, fill="blue") +
          geom_line(data = pred[[3]], aes(x = Year, y = NoExSpec), col="green")) +
          geom_ribbon(data = pred[[3]], aes(x = Year, ymin=LowerCI, ymax=UpperCI), alpha=0.2, fill="green") +
          theme_bw()
  dev.off()
}


















# set best model to use for predictions
model <- m3a#extot[[4]]


pred_2050 <- predict(model, popsl95, type = "raw")

# finding the confidence intervals

ci_2100 <- as.data.frame(predict(model, popsl95, se.fit=TRUE, type = "raw"))

ci_2050 <- as.data.frame(predict(model$finalModel, popsl95, se.fit=TRUE, type = "link"))

(ci_2050$fit - 2*ci_2050$se.fit)^2
(ci_2050$fit + 2*ci_2050$se.fit)^2


exp(ci_2100$fit - 2*ci_2100$se.fit)
exp(ci_2100$fit + 2*ci_2100$se.fit)



# extract populations
pops <- as.data.frame(dataset$Pop)
names(pops) <- "Pop"

# make precitions
if (class(model)[1] == "train") {
  preds <- predict(model, pops, type = "raw")
} else {
  preds <- predict(model, pops, type = "response")
}

# combine pop dataset with predictions
preddf <- cbind(dataset, preds)
names(preddf) <- c("Year", "Pop", "NoExSpec")

# finding the confidence intervals
if (class(model)[1] == "train") {
  ci <- predict(model$finalModel, preddf, se.fit=TRUE, type = "link")
} else {
  ci <- predict(model, preddf, se.fit=TRUE, type = "link")
}

# lower limit
cil <- (ci$fit - 2*ci$se.fit)^2
# upper limit
ciu <- (ci$fit + 2*ci$se.fit)^2

# combine all data and confidence intervals for plotting
fulldf <- cbind(dataset, preds, cil, ciu)
names(fulldf) <- c("Year", "Pop", "NoExSpec", "UpperCI", "LowerCI")

return(fulldf)

pred <- vector(mode="list", length=length(dfs))
for (df in 1:length(dfs)) {
  full_preds <- make_predictions(dfs[[df]], model_ls[[mod]])
  write.csv(full_preds, paste0(path_out, "extot_", mod, "_", df, "_preds.csv"))
  pred[[df]] <- full_preds
}
pdf(file=paste0(path_out, mod, "_extot_preds.pdf"))
print(ggplot() +
        geom_point(data = expoptot, aes(x = Year, y = NoExSpec), col="black") +
        geom_line(data = pred[[1]], aes(x = Year, y = NoExSpec), col="red") +
        geom_line(data = pred[[2]], aes(x = Year, y = NoExSpec), col="blue") +
        geom_line(data = pred[[3]], aes(x = Year, y = NoExSpec), col="green"))
dev.off()
######## notes

# use negative binomial and quasipoisson since mean is not 
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


## end of script

