###########################################################################################
############################ Statistical analysis for #####################################
####################### Number of vertebrate extinctions #################################
###########################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Created: Feb 2022
# Last updated: Mar 2022
# based on script by Jack V Johnson

# clear workspace
rm(list=ls())


########################################################################################
################################### Set up #############################################


# load required packages
library("dplyr")
library("tidyr")
library("MASS") # for negative binomial
library("caret") # for cross-validation
library("pscl") # for zero-inflated models
library("Kendall") # for Mann-Kendall test

# set seed for reproducibility
set.seed(26)

# select time name to analyse data from set time periods
time_name <- "1460-2000_20y"
binned_years <- seq(1460, 2000, 20)
#year_split <- 5

# set path to data and scripts
path_data <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/r_data_objects/")
path_rawdata <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/data/raw_data/")
path_out_temp <- paste0("~/Dropbox/luke/documents/academia/phd/papers/2022_global_extinctions/outputs/")

# create directory to store the data for this time set
# first check if directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(path_out_temp, time_name))), dir.create(file.path(paste0(path_out_temp, time_name))), FALSE)
# set directory for results to be sent to
path_out_time <- paste0(path_out_temp, time_name, "/")

# create directory to store the data for these models
# first check if directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(path_out_time, "ex_time_glms/"))), dir.create(file.path(paste0(path_out_time, "ex_time_glms/"))), FALSE)
# set directory for results to be sent to
path_out <- paste0(path_out_time, "ex_time_glms/")

# define 3-fold cross validation with 3 repeats
train_control <- trainControl(method="repeatedcv", number=5, repeats=5, seeds=as.list(54:80))
seq <- seq(54:(26*3**3))
train_control_nb <- trainControl(method="repeatedcv", number=5, repeats=5, seeds=split(seq, ceiling(seq_along(seq)/3)))

  
###################################################################################
############################## Data wrangling ######################################


# load data
# extinction/population data
vertex_full <- read.csv(paste0(path_rawdata, "vert_extinctions.csv"))

# remove uncertain extinction time intervals
vertex <- vertex_full[which(!(vertex_full$Century == 2000)),]
vertex <- vertex[which(!(vertex$EX.Last.seen. == "1700-1750")),]
vertex <- vertex[which(!(vertex$EX.Last.seen == "1500-1600")),]
vertex <- vertex[!grepl("s", vertex$EX.Last.seen),]
vertex <- vertex[!grepl("century", vertex$EX.Last.seen),]
vertex <- vertex[!grepl("before", vertex$EX.Last.seen),]
vertex <- vertex[!grepl("approx", vertex$EX.Last.seen),]

# convert extinctions to numeric
vertex$EX.Last.seen.<- as.numeric(vertex$EX.Last.seen.)

# group by specified time periods
vertex$Year_Block_Var <- NA
for (t in 2:length(binned_years)-1) {
  for (spec in 1:nrow(vertex)) {
    if (binned_years[t] < vertex$EX.Last.seen.[spec] && vertex$EX.Last.seen.[spec] <= binned_years[t+1]) {
      vertex$Year_Block_Var[spec] <- binned_years[t+1]  
    }
  }
}

# remove species that went extinct outside time period bounds
vertex <- vertex[which(!is.na(vertex$Year_Block_Var)),]

# set to match Jack's code
vert_extinct <- vertex

# replace continent/island with continent since island effects won't be the same
table(vert_extinct$Cont.2..Island.1.)
vert_extinct[] <- lapply(vert_extinct, gsub, pattern='Continent/Island', replacement='Continent')

# bin over year breaks
#binned_years <- .bincode(vert_extinct$EX.Last.seen., binned_years, right = T, include.lowest = T)
#binned_years <- as.vector(binned_years) # set as vector
#vert_extinct <- cbind(vert_extinct, binned_years) # combine vertex data with binned years
#summary(vert_extinct$binned_years)


###################################################################################
############################### Functions #########################################


# function to add missing time periods with 0 extinctions
add_zeroes <- function(data) {
  # initiliase data frame with binned years
  data0 <- data.frame(binned_years)
  # set count to zero for all years
  data0$n <- 0
  # add extinctions to dataframe
  for (i in (1:nrow(data0))) {
    if (length(data$n[which(data$Year_Block_Var == data0[i,1])]) > 0) {
      data0$n[i] <- data$n[which(data$Year_Block_Var == data0[i,1])]
    }
  }
  # rescale year bins
  years <- seq(1, nrow(data0), 1)
  data0$years <- years
  return(data0)
}

# function for running and saving model results
run_models <- function(dataset, x, y, outliers=NA, output_name, region, taxa) {
  
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
  
  # set max iterations for negative binomial models
  maxset = 1000
  
  
  # MODEL 1: Quasipoisson GLM
  
  
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
  
  
  # save models as list
  models <- list("mglmcvlog" = mglmcvlog,
                 "mglmnbcv" = mglmnbcv)
  
  # generate csv of RMSE and R^2 values
  metrics <- data.frame(matrix(ncol=10,nrow=length(models)))
  names(metrics) <- c("Region", "Taxa", "Model", "Link", "Estimate", "Std. Error", "z/t value", "Pr(>|t|)", "RMSE", "R^2")
  
  for (model in 1:length(models)) {
    
    metrics$Region <- region
    metrics$Taxa <- taxa
    
    if (names(models)[[model]] == "mglmcvlog") {
      metrics$Model[model] <- "Quasipoisson GLM"
      metrics$Link[model] <- "log"
    }
    if (names(models)[[model]] == "mglmnbcv") {
      metrics$Model[model] <- "Negative Binomial GLM"
      metrics$Link[model] <- models[[model]]$results$link
    }
  
    metrics$Estimate[[model]] <- coef(summary(models[[model]]))[2,1]
    metrics$`Std. Error`[[model]] <- coef(summary(models[[model]]))[2,2]
    metrics$`z/t value`[[model]] <- coef(summary(models[[model]]))[2,3]
    metrics$`Pr(>|t|)`[[model]] <- coef(summary(models[[model]]))[2,4]
    if (!is.null(models[[model]]$results$parameter)) {
      metrics$RMSE[model] <- models[[model]]$results[[2]]
      metrics$`R^2`[[model]] <- models[[model]]$results[[3]]
      metrics$Estimate[[model]] <- models[[model]]$finalModel$coefficients[["years"]]
    } else {
      metrics$RMSE[model] <- models[[model]]$results[[2]][[which(models[[model]]$results$link==models[[model]]$finalModel$family$link)]]
      metrics$`R^2`[[model]] <- models[[model]]$results[[3]][[which(models[[model]]$results$link==models[[model]]$finalModel$family$link)]]
    }
  }
  
  # save as csv
  write.csv(metrics, paste0(path_out_mod, output_name, "_model_comparison.csv"))
  
  # return models
  return(metrics)
  
}

# function to run zero-inflated models and save results
run_zeroinf <- function(dataset, x, y, outliers=NA, output_name, region, taxa) {
  
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
  
  # run zero-inflated model
  z <- zeroinfl(form, data=dataset, dist="negbin")
  
  # save outputs
  cat("\n\n Zero-inflated negative binomial for ", outliers, " outliers\n", 
      capture.output(summary(z)), 
      file=paste0(path_out_mod, output_name, "_GLM.txt"), 
      sep = "\n", append=TRUE)
  
  # generate csv of RMSE and R^2 values
  metrics <- data.frame(matrix(ncol=11,nrow=1))
  names(metrics) <- c("Region", "Taxa", "Model", "Estimate", "Std. Error", "z/t value", "Pr(>|t|)", "Estimate0", "Std. Error0", "z/t value0", "Pr(>|t|)0")
  
  metrics$Region <- region
  metrics$Taxa <- taxa
  metrics$Model <- "Zero-inflated Negative Binomial"
    
  metrics$Estimate <- coef(summary(z))[[1]][2,1]
  metrics$`Std. Error` <- coef(summary(z))[[1]][2,2]
  metrics$`z/t value` <- coef(summary(z))[[1]][2,3]
  metrics$`Pr(>|t|)` <- coef(summary(z))[[1]][2,4]
  
  metrics$Estimate0 <- coef(summary(z))[[2]][2,1]
  metrics$`Std. Error0` <- coef(summary(z))[[2]][2,2]
  metrics$`z/t value0` <- coef(summary(z))[[2]][2,3]
  metrics$`Pr(>|t|)0` <- coef(summary(z))[[2]][2,4]
  
  # return models
  return(metrics)
  
}

# function to run Mann-Kendall test
run_mk <- function(dataset, y, region, taxa) {
  
  # INPUTS
  # dataset must be dataframe
  # y must be name of response column as string
  
  str(dataset)
  # convert to time series
  data_time <- ts(dataset[[y]], start=1)
  str(data_time)
  
  # check for autocorrelation and partial autocorrelation
  # (spikes outside blue dotted band are significant autocorrelations and partial autocorrelations)
  par(mfrow=c(2,1))
  acf(data_time)
  pacf(data_time)
  # no significant autocorrelation or partial autocorrelation
  
  # run Mann-Kendall test
  mk <- MannKendall(data_time)
  summary(mk)
  
  # store results
  mkres <- c(region, taxa, mk[["tau"]][[1]], mk[["sl"]][[1]])
  
  return(mkres)
  
}


################################################################################
########################## Total extinctions ###################################


#  group by year
vertex_tot <- vert_extinct %>% group_by(Year_Block_Var) %>% summarise(n = n())
# add missing zeroes
vertex0 <- add_zeroes(vertex_tot)

# view data
#hist(vertex0$n)
#table(vertex0$n)
#plot(vertex0$binned_years, vertex0$n)

# run models
#totalex <- run_models(vertex0, "years", "n", outliers=NA, "totalex", "Global", "All Taxa")
#totalex0 <- run_zeroinf(vertex0, "years", "n", outliers=NA, "totalex0", "Global", "All Taxa")
totalex_mk0 <- run_mk(vertex0, "n", "Global", "All Taxa")


################################################################################
######################### Island vs continental ################################


# summarise over island/continental
vertex_tot_contis <- vert_extinct %>% group_by(Year_Block_Var, Cont.2..Island.1.) %>% summarise(n = n())

# remove any NA records
vertex_tot_contis <- vertex_tot_contis[!is.na(vertex_tot_contis$Cont.2..Island.1.),]


################################# Island ###########################################


# subset by island
df_island <- subset(vertex_tot_contis, Cont.2..Island.1.== 'Island', select =c(Year_Block_Var, n))
# add mising zeroes
island0 <- add_zeroes(df_island)

# view data
#hist(island0$n)
#table(island0$n)
#plot(island0$binned_years, island0$n)

# run models
#totalex_is <- run_models(island0, "years", "n", outliers=NA, "totalex_is", "Island", "All Taxa")
#totalex_is0 <- run_zeroinf(island0, "years", "n", outliers=NA, "totalex_is0", "Island", "All Taxa")
totalex_is_mk0 <- run_mk(island0, "n", "Island", "All Taxa")


############################## Continent ##########################################


# subset by continent
df_cont <- subset(vertex_tot_contis, Cont.2..Island.1.== 'Continent', select =c(Year_Block_Var, n))                    
# add mising zeroes
cont0 <- add_zeroes(df_cont)

# view data
#hist(cont0$n)
#table(cont0$n)
#plot(cont0$binned_years, cont0$n)

# run models
#totalex_cont <- run_models(cont0, "years", "n", outliers=NA, "totalex_cont", "Continent", "All Taxa")
#totalex_cont0 <- run_zeroinf(cont0, "years", "n", outliers=NA, "totalex_cont0", "Continent", "All Taxa")
totalex_cont_mk0 <- run_mk(cont0, "n", "Continent", "All Taxa")


####################################################################################
################################# Class ############################################


# group by class
vertex_class <- vert_extinct %>% group_by(Year_Block_Var, Class) %>% summarise(n = n())

# group by class and island/continent
vertex_class_contis <- vert_extinct %>% group_by(Year_Block_Var, Cont.2..Island.1., Class) %>% summarise(n = n())


############################## Amphibians ##########################################


# subset by amphibians
df_amphib <- subset(vertex_class, Class == 'Amphibia', select =c(Year_Block_Var, n))
# add mising zeroes
amph0 <- add_zeroes(df_amphib)

# view data
#hist(amph0$n)
#table(amph0$n)
#plot(amph0$binned_years, amph0$n)

# subset amphibians by island
df_amphib_island <- subset(vertex_class_contis, Cont.2..Island.1. == 'Island'& Class =='Amphibia', select =c(Year_Block_Var, n))
# add mising zeroes
amphis0 <- add_zeroes(df_amphib_island)
# view data
#hist(amphis0$n)
#table(amphis0$n)
#plot(amphis0$binned_years, amphis0$n)

# subset amphibians by continent
df_amphib_cont <- subset(vertex_class_contis, Cont.2..Island.1.== 'Continent' & Class =='Amphibia', select =c(Year_Block_Var, n))
# add mising zeroes
amphcont0 <- add_zeroes(df_amphib_cont)
# view data
#hist(amphcont0$n)
#table(amphcont0$n)
#plot(amphcont0$binned_years, amphcont0$n)

# run models
# all amphibians
#amph <- run_models(amph0, "years", "n", outliers=NA, "amph", "Global", "Amphibia")
#amphmod0 <- run_zeroinf(amph0, "years", "n", outliers=NA, "amph0", "Global", "Amphibia")
amphmod_mk0 <- run_mk(amph0, "n", "Global", "Amphibia")
# island only
#amphis <- run_models(amphis0, "years", "n", outliers=NA, "amphis", "Island", "Amphibia")
#amphismod0 <- run_zeroinf(amphis0, "years", "n", outliers=NA, "amphis0", "Island", "Amphibia")
amphismod_mk0 <- run_mk(amphis0, "n", "Island", "Amphibia")
# continent only
#amphcont <- run_models(amphcont0, "years", "n", outliers=NA, "amphcont", "Continent", "Amphibia")
#amphcontmod0 <- run_zeroinf(amphcont0, "years", "n", outliers=NA, "amphcont0", "Continent", "Amphibia")
amphcontmod_mk0 <- run_mk(amphcont0, "n", "Continent", "Amphibia")


################################## Aves ############################################


# subset by aves
df_aves <- subset(vertex_class, Class == 'AVES', select =c(Year_Block_Var, n))
# add mising zeroes
aves0 <- add_zeroes(df_aves)
# view data
#hist(aves0$n)
#table(aves0$n)
#plot(aves0$binned_years, aves0$n)

# subset aves by island
df_aves_island <- subset(vertex_class_contis, Cont.2..Island.1.== 'Island' & Class =='AVES', select =c(Year_Block_Var, n))
# add missing zeroes
avesis0 <- add_zeroes(df_aves_island)
# view data
#hist(avesis0$n)
#table(avesis0$n)
#plot(avesis0$binned_years, avesis0$n)

# subset aves by continent
df_aves_cont <- subset(vertex_class_contis, Cont.2..Island.1.== 'Continent' & Class == 'AVES', select =c(Year_Block_Var, n))
# add missing zeroes
avescont0 <- add_zeroes(df_aves_cont)
# view data
#hist(avescont0$n)
#table(avescont0$n)
#plot(avescont0$binned_years, avescont0$n)

# run models
#aves <- run_models(aves0, "years", "n", outliers=NA, "aves", "Global", "Aves")
#avesmod0 <- run_zeroinf(aves0, "years", "n", outliers=NA, "aves0", "Global", "Aves")
avesmod_mk0 <- run_mk(aves0, "n", "Global", "Aves")
# island only
#avesis <- run_models(avesis0, "years", "n", outliers=NA, "avesis", "Island", "Aves")
#avesismod0 <- run_zeroinf(avesis0, "years", "n", outliers=NA, "avesis0", "Island", "Aves")
avesismod_mk0 <- run_mk(avesis0, "n", "Island", "Aves")
# continent only
#avescont <- run_models(avescont0, "years", "n", outliers=NA, "avescont", "Continent", "Aves")
#avescontmod0 <- run_zeroinf(avescont0, "years", "n", outliers=NA, "avescont0", "Continent", "Aves")
avescontmod_mk0 <- run_mk(avescont0, "n", "Continent", "Aves")


############################### Mammalia ###########################################


# subset by mammals
df_mammal <- subset(vertex_class, Class == 'MAMMALIA', select =c(Year_Block_Var, n))
# add missing zeroes
mamm0 <- add_zeroes(df_mammal)

# view data
#hist(mamm0$n)
#table(mamm0$n)
#plot(mamm0$binned_years, mamm0$n)

# subset mammals by island
df_mammal_island <- subset(vertex_class_contis, Cont.2..Island.1.== 'Island' & Class == 'MAMMALIA', select =c(Year_Block_Var, n))
# add missing zeroes
mammis0 <- add_zeroes(df_mammal_island)
# view data
#hist(mammis0$n)
#table(mammis0$n)
#plot(mammis0$binned_years, mammis0$n)

# subset mammals by continent
df_mammal_cont <- subset(vertex_class_contis, Cont.2..Island.1.== 'Continent' & Class =='MAMMALIA', select =c(Year_Block_Var, n))
# add missing zeroes
mammcont0 <- add_zeroes(df_mammal_cont)
# view data
#hist(mammcont0$n)
#table(mammcont0$n)
#plot(mammcont0$binned_years, mammcont0$n)

# run models
#mamm <- run_models(mamm0, "years", "n", outliers=NA, "mamm", "Global", "Mammalia")
#mammmod0 <- run_zeroinf(mamm0, "years", "n", outliers=NA, "mamm0", "Global", "Mammalia")
mammmod_mk0 <- run_mk(mamm0, "n", "Global", "Mammalia")
# island only
#mammis <- run_models(mammis0, "years", "n", outliers=NA, "mammis", "Island", "Mammalia")
#mammismod0 <- run_zeroinf(mammis0, "years", "n", outliers=NA, "mammis0", "Island", "Mammalia")
mammismod_mk0 <- run_mk(mammis0, "n", "Island", "Mammalia")
# continent only
#mammcont <- run_models(mammcont0, "years", "n", outliers=NA, "mammcont", "Continent", "Mammalia")
#mammcontmod0 <- run_zeroinf(mammcont0, "years", "n", outliers=NA, "mammcont0", "Continent", "Mammalia")
mammcontmod_mk0 <- run_mk(mammcont0, "n", "Continent", "Mammalia")


############################## Reptilia ############################################


# subset by reptiles
df_reptile <- subset(vertex_class, Class == 'Reptilia', select =c(Year_Block_Var, n))
# add missing zeroes
rept0 <- add_zeroes(df_reptile)

# view data
#hist(rept0$n)
#table(rept0$n)
#plot(rept0$binned_years, rept0$n)

# subset reptiles by island
df_reptile_island <- subset(vertex_class_contis, Cont.2..Island.1.== 'Island' & Class == 'Reptilia', select =c(Year_Block_Var, n))          
# add missing zeroes
reptis0 <- add_zeroes(df_reptile_island)
# view data
#hist(reptis0$n)
#table(reptis0$n)
#plot(reptis0$binned_years, reptis0$n)

# subset reptiles by continent
df_reptile_cont  <- subset(vertex_class_contis, Cont.2..Island.1.== 'Continent' & Class == 'Reptilia', select =c(Year_Block_Var, n))
# add missing zeroes
reptcont0 <- add_zeroes(df_reptile_cont)
# view data
#hist(reptcont0$n)
#table(reptcont0$n)
#plot(reptcont0$binned_years, reptcont0$n)

# run models
#rept <- run_models(rept0, "years", "n", outliers=NA, "rept", "Global", "Reptilia")
#reptmod0 <- run_zeroinf(rept0, "years", "n", outliers=NA, "rept0", "Global", "Reptilia")
reptmod_mk0 <- run_mk(rept0, "n", "Global", "Reptilia")
# island only
#reptis <- run_models(reptis0, "years", "n", outliers=NA, "reptis", "Island", "Reptilia")
#reptismod0 <- run_zeroinf(reptis0, "years", "n", outliers=NA, "reptis0", "Island", "Reptilia")
reptismod_mk0 <- run_mk(reptis0, "n", "Island", "Reptilia")
# continent only
#reptcont <- run_models(reptcont0, "years", "n", outliers=NA, "reptcont", "Continent", "Reptilia")
#reptcontmod0 <- run_zeroinf(reptcont0, "years", "n", outliers=NA, "reptcont0", "Continent", "Reptilia")
reptcontmod_mk0 <- run_mk(reptcont0, "n", "Continent", "Reptilia")


#####################################################################################
############################## Combine results ######################################


#final_results <- rbind(totalex, amph, aves, mamm, rept,
#                       totalex_is, amphis, avesis, mammis, reptis,
#                       totalex_cont, amphcont, avescont, mammcont, reptcont)
#final_results <- rbind(totalex, amph, aves, mamm, rept,
#                       totalex_is, amphis, avesis, avesis, reptis,
#                       totalex_cont, amphcont)
#write.csv(final_results, paste0(path_out, "final_GLM_results.csv"))


#final_results0 <- rbind(totalex0,amphmod0, avesmod0, mammmod0, reptmod0,
#                       totalex_is0, amphismod0, avesismod0, mammismod0, reptismod0,
#                       totalex_cont0, amphcontmod0, avescontmod0, mammcontmod0, reptcontmod0)
#write.csv(final_results0, paste0(path_out, "final_zeroinf_results.csv"))

final_results_mk <- as.data.frame(rbind(totalex_mk0,amphmod_mk0, avesmod_mk0, mammmod_mk0, reptmod_mk0,
                        totalex_is_mk0, amphismod_mk0, avesismod_mk0, mammismod_mk0, reptismod_mk0,
                        totalex_cont_mk0, amphcontmod_mk0, avescontmod_mk0, mammcontmod_mk0, reptcontmod_mk0))
# create dataframe to save Mann-Kendall results
#mk_final <- data.frame(matrix(ncol=4))
names(final_results_mk) <- c("Region", "Taxa", "Tau", "2-sided P-value")
write.csv(final_results_mk, paste0(path_out, "final_mannkendall_results.csv"))

# 1460-1760 since many zeroes (0 amph ex, 0 rept ex, 0 cont bird ex)
#final_results_mk_pre1760 <- as.data.frame(rbind(totalex_mk0, avesmod_mk0, mammmod_mk0,
#                                        totalex_is_mk0, avesismod_mk0, mammismod_mk0,
#                                        totalex_cont_mk0, mammcontmod_mk0))
#names(final_results_mk_pre1760) <- c("Region", "Taxa", "Tau", "2-sided P-value")
#write.csv(final_results_mk_pre1760, paste0(path_out, "final_mannkendall_results.csv"))


#####################################################################################
######################## Pre- and post-industrial age onset #########################


############################## Total extinctions ####################################


# find scaled year for 1760
indrev <- vertex0$years[which(vertex0$binned_years == 1760)]

# subset by pre- and po-st industrial age onset
df_glob_pre_1760 <- subset(vertex0, years <= indrev, select =c(years, n))
df_glob_post_1760 <- subset(vertex0, years > indrev, select =c(years, n))

# set up table for chi-squared test
ind_table <- data.frame(matrix(ncol=2,nrow=2))
# if two groups are the same then the expected means are the means of the means
names(ind_table) <- c("Observed", "Expected")
ind_table[1,1] <- mean(df_glob_pre_1760$n)
ind_table[2,1] <- mean(df_glob_post_1760$n)
ind_table[1,2] <- (mean(df_glob_pre_1760$n) + mean(df_glob_post_1760$n))/2
ind_table[2,2] <- (mean(df_glob_pre_1760$n) + mean(df_glob_post_1760$n))/2

# save chi squared table as csv
write.csv(ind_table, paste0(path_out, "chi_sq_tab.csv"))

# run chisquare test
chi_tot <- chisq.test(t(ind_table), correct=FALSE)

# save output
cat("Chi-squared for global: ", capture.output(chi_tot), 
    file=paste0(path_out, "Totalex_chisq_indrev.txt"), 
    sep = "\n", append=TRUE)


################################ Island ###########################################


#names(cont0)[2] <- "Continent"
#names(island0)[2] <- "Island"
#df_contis <- merge(cont0, island0, by.x = "binned_years", by.y="binned_years")

# subset by pre- and po-st industrial age onset
df_is_pre_1760 <- subset(island0, years <= indrev, select =c(years, n))
df_is_post_1760 <- subset(island0, years > indrev, select =c(years, n))

# set up table for chi-squared test
ind_table_is <- data.frame(matrix(ncol=2,nrow=2))
# if two groups are the same then the expected means are the means of the means
names(ind_table_is) <- c("Observed", "Expected")
ind_table_is[1,1] <- mean(df_is_pre_1760$n)
ind_table_is[2,1] <- mean(df_is_post_1760$n)
ind_table_is[1,2] <- (mean(df_is_pre_1760$n) + mean(df_is_post_1760$n))/2
ind_table_is[2,2] <- (mean(df_is_pre_1760$n) + mean(df_is_post_1760$n))/2

# save chi squared table as csv
write.csv(ind_table_is, paste0(path_out, "chi_sq_is_tab.csv"))

# run chisquare test
chi_is <- chisq.test(t(ind_table_is), correct=FALSE)

# save output
cat("Chi-squared for islands: ", capture.output(chi_is), 
    file=paste0(path_out, "Totalex_chisq_indrev.txt"), 
    sep = "\n", append=TRUE)


############################## Continent #########################################


# subset by pre- and po-st industrial age onset
df_cont_pre_1760 <- subset(cont0, years <= indrev, select =c(years, n))
df_cont_post_1760 <- subset(cont0, years > indrev, select =c(years, n))

# set up table for chi-squared test
ind_table_cont <- data.frame(matrix(ncol=2,nrow=2))
# if two groups are the same then the expected means are the means of the means
names(ind_table_cont) <- c("Observed", "Expected")
ind_table_cont[1,1] <- mean(df_cont_pre_1760$n)
ind_table_cont[2,1] <- mean(df_cont_post_1760$n)
ind_table_cont[1,2] <- (mean(df_cont_pre_1760$n) + mean(df_cont_post_1760$n))/2
ind_table_cont[2,2] <- (mean(df_cont_pre_1760$n) + mean(df_cont_post_1760$n))/2

# save chi squared table as csv
write.csv(ind_table_cont, paste0(path_out, "chi_sq_cont_tab.csv"))

# run chisquare test
chi_cont <- chisq.test(ind_table_cont, correct=FALSE)

# save output
cat("Chi-squared for continents: ", capture.output(chi_cont), 
    file=paste0(path_out, "Totalex_chisq_indrev.txt"), 
    sep = "\n", append=TRUE)


##################################################################################
####################### Recent amphibian extinctions #############################


# find scaled year for 1950 and 1975
pre <- amph0$years[which(amph0$binned_years == 1950)]
post <- amph0$years[which(amph0$binned_years == 1975)]

# subset by 1950-1975 and 1975-2000 where 1975 is start of amphibian spike
df_1860 <- subset(amph0, ((pre < years) & (years <= post)), select =c(years, n))
df_1975 <- subset(amph0, ((post < years) & (years < length(amph0$years))), select =c(years, n))

# set up table for chi-squared test
ind_table_amph <- data.frame(matrix(ncol=2,nrow=2))
# if two groups are the same then the expected means are the means of the means
names(ind_table_amph) <- c("Observed", "Expected")
ind_table_amph[1,1] <- mean(df_1860$n)
ind_table_amph[2,1] <- mean(df_1975$n)
ind_table_amph[1,2] <- (mean(df_1860$n) + mean(df_1975$n))/2
ind_table_amph[2,2] <- (mean(df_1860$n) + mean(df_1975$n))/2

# save chi squared table as csv
write.csv(ind_table_amph, paste0(path_out, "chi_sq_amph_tab.csv"))

# run chisquare test
chi_amph <- chisq.test(ind_table_amph, correct=FALSE)

# save output
cat("Chi-squared for Amphibians: ", capture.output(chi_amph), 
    file=paste0(path_out, "Amph_chisq.txt"), 
    sep = "\n", append=TRUE)


## end of script
