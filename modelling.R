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
theme_set(theme_bw())
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             panel.border=element_blank(),
             axis.line = element_line(colour = "black"))
library("gridExtra")
library("MASS") # negative binomial GLM
library("car") # Anova() for a model
library("mgcv") # for gam
library("caret") # for cross validation

# set seed for reproducibility
set.seed(26)

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
# by class
expopam <- readRDS(paste0(path_data, time_name, "_expopam.rds"))
expopav <- readRDS(paste0(path_data, time_name, "_expopav.rds"))
expopma <- readRDS(paste0(path_data, time_name, "_expopma.rds"))
expopre <- readRDS(paste0(path_data, time_name, "_expopre.rds"))
# un population prediction data
un_pred_full <- read.csv(paste0(path_rawdata, "un_pop_pred_world.csv"))
# un data is in 1000s so need to multiply by 1000
un_pred_full[,2:ncol(un_pred_full)] <- un_pred_full[,2:ncol(un_pred_full)]*1000
# load population raw data to fill in missing years
pop_full <- read.csv(paste0(path_rawdata, "population_past_future.csv"), stringsAsFactors = F)
names(pop_full)[4] <- "Pop"

# define 3-fold cross validation with 3 repeats
train_control <- trainControl(method="repeatedcv", number=5, repeats=5, seeds=as.list(54:80))
seq <- seq(54:(26*3**3))
train_control_nb <- trainControl(method="repeatedcv", number=5, repeats=5, seeds=split(seq, ceiling(seq_along(seq)/3)))


#####################################################################################
#################################### Models #########################################


# from data exploration, we know data has rough Poisson distribution and that
# mean is not equal to variance so quasipoisson or negative binomial distributions must be used

# Models used:
# BASELINE and with CROSS-VALIDATION fro GLMs
# 1) mglm and mglmcv are GLMs using quasipoisson family with log link function
# 2) mglmnb and mglmnbcv are negative binomial GLMs
# BASELINE only for GAMs (as difficult to use smoothers)
# 3) mgam is a GAM using quasipoisson family with log link function
# 4) mgamnb is a negative binomial GAM

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
  
  # create directory to store the data for this time set
  # first check if directory exists and if not, create it
  ifelse(!dir.exists(file.path(paste0(path_out, output_name))), dir.create(file.path(paste0(path_out, output_name))), FALSE)
  # set directory for results to be sent to
  path_out_mod <- paste0(path_out, output_name, "/")
  
  # create dataset without outliers
  if (is.na(outliers)) {
    df <- dataset
  } else {df <- dataset[which(!(dataset$Year %in% outliers)),]}
  
  # set formula for model creation
  form <- as.formula(paste(y, '~', x))
  forms <- as.formula(paste0(y, ' ~ s(' , x, ')'))
  
  # set max iterations for negative binomial models
  maxset = 1000
  
  
  # MODEL 1: Quasipoisson GLM
  
  
  ## cross-validation for identity link
  mglmcvid <- train(form,
                  method = "glm",
                  family = quasipoisson(link = "identity"),
                  data=df,
                  maxit = maxset,
                  trControl=train_control)
  # save outputs
  cat("\n\n Quasipoisson GLM (identity link) with cross-validation for ", outliers, " outliers\n", 
      capture.output(summary(mglmcvid)), 
      file=paste0(path_out_mod, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nCV summary\n", capture.output(print(mglmcvid)), 
      file=paste0(path_out_mod, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  ## cross-validation for log link
  mglmcvlog <- train(form,
              method = "glm",
              family = quasipoisson(link = "log"),
              data=df,
              maxit = maxset,
              trControl=train_control)
  # save outputs
  cat("\n\n Quasipoisson GLM (log link) with cross-validation for ", outliers, " outliers\n", 
      capture.output(summary(mglmcvlog)), 
      file=paste0(path_out_mod, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nCV summary\n", capture.output(print(mglmcvlog)), 
      file=paste0(path_out_mod, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  
  
  # MODEL 2: Negative Binomial GLM
  
  
  ## base
  mglmnb <- glm.nb(form, data=df, link="identity", maxit=maxset)
  r2_mglmnb <- with(summary(mglmnb), 1 - deviance/null.deviance) # calculate r^2 value
  # save outputs
  cat("\n\n Negative Binomial GLM for ", outliers, " outliers\n", 
      capture.output(summary(mglmnb)), 
      file=paste0(path_out_mod, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nR^2", capture.output(r2_mglmnb), 
      file=paste0(path_out_mod, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  ## cross-validation
  mglmnbcv <- train(form, 
               method = "glm.nb", 
               data=df,
               maxit = maxset,
               trControl=train_control_nb)
  r2_mglmnbcv <- with(summary(mglmnbcv), 1 - deviance/null.deviance) # calculate r^2 value
  # save outputs
  cat("\n\n Negative Binomial GLM with cross-validation for ", outliers, " outliers\n", 
      capture.output(summary(mglmnbcv)), 
      file=paste0(path_out_mod, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nR^2", capture.output(r2_mglmnbcv), 
      file=paste0(path_out_mod, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
  cat("\n\nCV summary\n", capture.output(print(mglmnbcv)), 
      file=paste0(path_out_mod, output_name, "_GLMnb.txt"), 
      sep = "\n", append=TRUE)
        
  
  # MODEL 3: Quasipoisson GAM
  
  
  # ## base
  # mgam <- gam(forms, family=quasipoisson(link="log"), data=df)
  # # save outputs
  # cat("\n\n Quasipoisson GAM for ", outliers, " outliers\n", 
  #     capture.output(summary(mgam)), 
  #     file=paste0(path_out_mod, output_name, "_GAM.txt"), 
  #     sep = "\n", append=TRUE)
    
    
  # MODEL 4: Negative Binomial GAM
  
  
  # ## base
  # mgamnb <- gam(forms, family=nb(), data=df)
  # save outputs
  # cat("\n\n Negative Binomial GAM for ", outliers, " outliers\n", 
  #    capture.output(summary(mgamnb)), 
  #    file=paste0(path_out_mod, output_name, "_GAMnb.txt"), 
  #    sep = "\n", append=TRUE)
  
  # save models as list
  models <- list("mglmcvid" = mglmcvid,
                 "mglmcvlog" = mglmcvlog,
                 #mglmnb = mglmnb,
                 "mglmnbcv" = mglmnbcv)
                 #mgam = mgam,
                 #mgamnb = mgamnb)
  
  # generate csv of RMSE and R^2 values
  metrics <- data.frame(matrix(ncol=3,nrow=length(models)))
  names(metrics) <- c("Model", "RMSE", "R^2")
  for (model in 1:length(models)) {
    metrics$Model[model] <- names(models)[[model]]
    if (!is.null(models[[model]]$results$parameter)) {
      metrics$RMSE[model] <- models[[model]]$results[[2]]
      metrics$`R^2`[[model]] <- models[[model]]$results[[3]]
    } else {
      metrics$RMSE[model] <- models[[model]]$results[[2]][[which(models[[model]]$results$link==models[[model]]$finalModel$family$link)]]
      metrics$`R^2`[[model]] <- models[[model]]$results[[3]][[which(models[[model]]$results$link==models[[model]]$finalModel$family$link)]]
    }
  }
  
  # save as csv
  write.csv(metrics, paste0(path_out_mod, "model_comparison.csv"))
  
  # add extra nb model to view graphically
  models[["mglmnb"]] = mglmnb
  
  # return models
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
# a) TOTAL EXTINCTIONS AND b) EACH CLASS SEPARATELY
# 1) 1805-2000, split by 5 years
# 2) 1805-2000, split by 5 years, with 1900 and 1955 outliers removed
# 3) 1805-2000, split by 5 years, with 1900, 1955 outliers removed and 1995, 2000 removed to remove downward trend in GAMs
# 4) 1710-2000, split by 10 years
# 5) 1710-2000, split by 10 years, with 1900 outlier removed
# 6) 1710-2000, split by 10 yearsm with 1900 outlier removed and 2000 removed to remove downward trend in GAMs

# note outliers were chosen from the plots in explore_data.R
# note correlation for 1805-2000 (5 years) was 0.76 so quite high (explore_data.R)
# and correlation for 1710-2000 (10 years) was 0.89 so very high (explore_data.R)

# set outliers depending on time period
if (year_split == 10) {
  outliersp <- c(1900)
  outliersr <- c(1900, 2000)
} 
if (year_split == 5) {
  outliersp <- c(1900, 1955)
  outliersr <- c(1900, 1955, 1995, 2000)
}


# a) TOTAL EXTINCTIONS

extot_outna <- run_models(expoptot, "Pop", "NoExSpec", NA, "extot_outna")
extot_outp <- run_models(expoptot, "Pop", "NoExSpec", outliersp, "extot_outp")
#extot_outr <- run_models(expoptot, "Pop", "NoExSpec", outliersr, "extot_outr")


# b) BY CLASS

## Amphibians
expopam_outna <- run_models(expopam, "Pop", "NoExSpec", NA, "expopam_outna")
expopam_outp <- run_models(expopam, "Pop", "NoExSpec", outliersp, "expopam_outp")

## Aves
expopav_outna <- run_models(expopav, "Pop", "NoExSpec", NA, "expopav_outna")
expopav_outp <- run_models(expopav, "Pop", "NoExSpec", outliersp, "expopav_outp")

## Mammals
expopma_outna <- run_models(expopma, "Pop", "NoExSpec", NA, "expopma_outna")
expopma_outp <- run_models(expopma, "Pop", "NoExSpec", outliersp, "expopma_outp")

## Reptiles
expopre_outna <- run_models(expopre, "Pop", "NoExSpec", NA, "expopre_outna")
expopre_outp <- run_models(expopre, "Pop", "NoExSpec", outliersp, "expopre_outp")



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
  cil <- ci$fit - 2*ci$se.fit
  # upper limit
  ciu <- ci$fit + 2*ci$se.fit
  
  if (class(model)[1] == "train") {
    cilf <- model$finalModel$family$linkinv(cil)
    ciuf <- model$finalModel$family$linkinv(ciu)
  } else {
    cilf <- model$family$linkinv(cil)
    ciuf <- model$family$linkinv(ciu)
  }
  
  # calculate confidence intervals based on link function
  #if (!is.null(model$bestTune[[1]])) {
  # if (model$bestTune[[1]] == "sqrt") {
  #   # lower limit
  #   cil <- (ci$fit - 2*ci$se.fit)^2
  #   # upper limit
  #   ciu <- (ci$fit + 2*ci$se.fit)^2
  # }
  # if (model$bestTune[[1]] == "identity") {
  #   # lower limit
  #   cil <- ci$fit - 2*ci$se.fit
  #   # upper limit
  #   ciu <- ci$fit + 2*ci$se.fit
  # }
  #if ((model$bestTune[[1]] == "log") | (model$bestTune[[1]] == "none")) {
  #   # lower limit
  #   cil <- exp(ci$fit - 2*ci$se.fit)
  #   # upper limit
  #   ciu <- exp(ci$fit + 2*ci$se.fit)
  # }
  #}
  #if (is.null(model$bestTune[[1]])) { # all non-cv models use log link
  #  # lower limit
  #  cil <- exp(ci$fit - 2*ci$se.fit)
  #  # upper limit
  #  ciu <- exp(ci$fit + 2*ci$se.fit)
  #}
  
  # combine all data and confidence intervals for plotting
  fulldf <- cbind(dataset, preds, cilf, ciuf)
  names(fulldf) <- c("Year", "Pop", "NoExSpec", "LowerCI", "UpperCI")
  
  return(fulldf)
  
}


########################################################################################
######################### Generate prediction datasets #################################


# create datasets for predictions
popsu95 <- process_data(pop_full, year_split, "Upper 95")
popsmed <- process_data(pop_full, year_split, "Median")
popsl95 <- process_data(pop_full, year_split, "Lower 95")

# create list of datasets
dfs <- list(u95=popsu95, med=popsmed, l95=popsl95)

# list of models to predict
model_ls <- list("expopav_outna" = expopav_outna, "expopav_outp" = expopav_outp,
                 "expopre_outna" = expopre_outna, "expopre_outp" = expopre_outp)

plot_data <- list(expopav, expopav, expopre, expopre)

# make predictions for each model for each dataset, save predictions as csv
# and plot each model on one plot for all three datasets
for (ds in (1:length(model_ls))) {
  # set directory for results to be sent to
  path_out_mod <- paste0(path_out, names(model_ls)[[ds]], "/")
  # set up correct data to plot
  if (grepl("na", names(model_ls)[[ds]], fixed=TRUE) == TRUE) {
    dataset <- plot_data[[ds]]
  }
  if (grepl("outp", names(model_ls)[[ds]], fixed = TRUE) == TRUE) {
    dataset <- plot_data[[ds]][which(!(plot_data[[ds]]$Year %in% outliersp)),]
  }
  if (grepl("outr", names(model_ls)[[ds]], fixed = TRUE) == TRUE) {
    dataset <- plot_data[[ds]][which(!((plot_data[[ds]]$Year %in% outliersr))),]
  }
  # run over all models in dataset
  for (mod in 1:length(model_ls[[ds]])) {
    pred <- vector(mode="list", length=length(dfs))
    for (df in 1:length(dfs)) {
      full_preds <- make_predictions(dfs[[df]], model_ls[[ds]][[mod]])
      write.csv(full_preds, paste0(path_out_mod, names(model_ls)[[ds]], "_", names(model_ls[[ds]])[mod], "_", names(dfs)[[df]], "_preds.csv"))
      pred[[df]] <- full_preds
    }
    pdf(file=paste0(path_out_mod, names(model_ls[[ds]])[mod], "_", names(model_ls)[[ds]], "_preds.pdf"), width=10, height=7)
    print(ggplot() +
            geom_point(data = dataset, aes(x = Year, y = NoExSpec), colour="black") +
            geom_line(data = pred[[1]], aes(x = Year, y = NoExSpec, colour="Upper 95% Population Prediction Interval")) +
            geom_ribbon(data = pred[[1]], aes(x = Year, ymin=LowerCI, ymax=UpperCI), alpha=0.15) +
            geom_line(data = pred[[2]], aes(x = Year, y = NoExSpec, colour="Median Population Prediction")) +
            geom_ribbon(data = pred[[2]], aes(x = Year, ymin=LowerCI, ymax=UpperCI), alpha=0.2) +
            geom_line(data = pred[[3]], aes(x = Year, y = NoExSpec, colour="Lower 95% Population Prediction Interval")) +
            geom_ribbon(data = pred[[3]], aes(x = Year, ymin=LowerCI, ymax=UpperCI), alpha=0.15) +
            labs(y=paste0("Number of extinct species per ", year_split, " 5 year period")) +
            scale_color_manual(name="Extinctions based on UN Human Population Predictions",
                               breaks=c("Lower 95% Population Prediction Interval", 
                                        "Median Population Prediction", 
                                        "Upper 95% Population Prediction Interval"),
                               values=c("Lower 95% Population Prediction Interval"="darkgreen", 
                                        "Median Population Prediction"="blue", 
                                        "Upper 95% Population Prediction Interval"="red")) +
            scale_x_continuous(limits=c(time_periods[1]-year_split,2100), expand=c(0,0)) +
            theme(legend.position = c(0.3, 0.8)) +
            theme(plot.margin = margin(0.5,1,0.5,0.5, "cm")))
    dev.off()
  }
}


#### confidence intervals work for GAMs??!?!
#### think of way to extract model comparison variables in csv for easy comparison

#### once code is completed, subset data for class and cont/island and run for each different thing
#### also run over the two different time period options


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

