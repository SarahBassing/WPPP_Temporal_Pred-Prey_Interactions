  #'  ------------------------------------------
  #'  Ungulate effects on latency of site use
  #'  WA Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Nov. 2022
  #'  ------------------------------------------
  #'  Read in time-between-detection data (times between detection of an ungulate
  #'  followed by a different ungulate species) and build generalized linear model 
  #'  w\ exponential distribution in JAGS to estimate latency of camera site use 
  #'  following detection of different ungulate. To be compared to mean latency of
  #'  prey detections following predators. 
  #'  
  #'  Time-between-detections data generated in TimeBetweenDetections.R script
  #'  -------------------------------------------
  
  #'  Clear workspace & load libraries
  rm(list = ls())
  
  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  
  #'  Read in data
  load("./Outputs/tbd_ungulate_2022-11-18.RData") #tbd_conspif_5min_2022-11-18 5-min interval
  tbd_ungulate <- tbd_ungulate %>%
    #'  Change units of time
    mutate(tbd_min_round = round(TimeSinceLastDet, 0),
           tbd_min = TimeSinceLastDet,
           tbd_hour = tbd_min/60,
           tbd_day = tbd_hour/24) %>%
    dplyr::select(-TimeSinceLastDet)
  
  #'  Filter images by ungulate species
  tbd_md_ung <- filter(tbd_ungulate, Species == "Mule Deer")
  tbd_elk_ung <- filter(tbd_ungulate, Species == "Elk")
  tbd_moose_ung <- filter(tbd_ungulate, Species == "Moose")
  tbd_wtd_ung <- filter(tbd_ungulate, Species == "White-tailed Deer")
  
  #'  Function to identify time cutoff between detections
  tbd_summary <- function(tbd, spp, quant) {
    #'  Plot frequency of time-between-detections (should look exponential)
    hist(tbd$tbd_day, breaks = 50, main = paste("Number of days between detections for", spp))
    
    #'  Review range of TBD values
    print("Quantiles of days between detections of ungualted followed by conspecific")
    print(quantile(tbd$tbd_day))
    #'  Review 90 - 99th quantiles- where are there gaps and outliers in the distribution?
    print("90th, 95th, and 99th quantiles of days between detections of ungulate followed by conspecific") 
    print(quantile(tbd$tbd_day, c(0.9, 0.95, 0.97, 0.99)))
    
    #'  Re-plot frequency of time-btwn-detections after removing extreme values
    short_tbd <- filter(tbd, tbd_day <= quantile(tbd$tbd_day, c(quant)))
    hist(short_tbd$tbd_day, breaks = 25, main = paste("Number of days between detections for", spp, "\nup to quantile =", quant))
    
    #'  Mean tbd
    print("Mean TBD for summer, fall, winter, spring")
    print(mean(short_tbd$tbd_day[short_tbd$Season == "Summer"]))
    print(mean(short_tbd$tbd_day[short_tbd$Season == "Fall"]))
    print(mean(short_tbd$tbd_day[short_tbd$Season == "Winter"]))
    print(mean(short_tbd$tbd_day[short_tbd$Season == "Spring"]))
    #'  Summary of observations in each season
    print("Total TBDs for each season")
    print(table(short_tbd$Season))
    
    #'  Return dataset after removing extreme values
    return(short_tbd)
  }
  tbd_md_ung_short <- tbd_summary(tbd_md_ung, spp = "mule deer", quant = 0.97)
  tbd_elk_ung_short <- tbd_summary(tbd_elk_ung, spp = "elk", quant = 0.97)
  tbd_moose_ung_short <- tbd_summary(tbd_moose_ung, spp = "moose", quant = 0.95)
  tbd_wtd_ung_short <- tbd_summary(tbd_wtd_ung, spp = "white-tailed deer", quant = 0.97)
  
  tbd_all_ung_short <- rbind(tbd_md_ung_short, tbd_elk_ung_short, tbd_moose_ung_short, tbd_wtd_ung_short)
  # write.csv(tbd_all_ung_short, "./Outputs/tbd_ungulate_NoOutliers.csv")
  
  #'  Quick summary stats
  #'  Range of values
  round(range(tbd_all_ung_short$tbd_min), 2) # minutes
  round(range(tbd_all_ung_short$tbd_min)/60, 2) # hours
  round(range(tbd_all_ung_short$tbd_min)/1440, 2) # days
  #'  Mean
  mean(tbd_all_ung_short$tbd_min) # minutes
  mean(tbd_all_ung_short$tbd_min)/60 # hours
  #'  Median
  median(tbd_all_ung_short$tbd_min) # minutes
  median(tbd_all_ung_short$tbd_min)/60 # hours
  #'  Standard error
  sd(tbd_all_ung_short$tbd_min)/(sqrt(nrow(tbd_all_ung_short))) # minutes
  sd(tbd_all_ung_short$tbd_min)/(sqrt(nrow(tbd_all_ung_short)))/60 # hours
  
  #'  Are there study area differences for the ungulate mean tbd?
  (mulie_meanOK <- mean(tbd_md_ung_short$tbd_min[tbd_md_ung_short$Study_Area == "OK"]))
  (mulie_meanNE <- mean(tbd_md_ung_short$tbd_min[tbd_md_ung_short$Study_Area == "NE"]))
  (wtd_meanOK <- mean(tbd_wtd_ung_short$tbd_min[tbd_wtd_ung_short$Study_Area == "OK"]))
  (wtd_meanNE <- mean(tbd_wtd_ung_short$tbd_min[tbd_wtd_ung_short$Study_Area == "NE"]))
  (moose_meanOK <- mean(tbd_moose_ung_short$tbd_min[tbd_moose_ung_short$Study_Area == "OK"]))
  (moose_meanNE <- mean(tbd_moose_ung_short$tbd_min[tbd_moose_ung_short$Study_Area == "NE"]))
  (elk_meanOK <- mean(tbd_elk_ung_short$tbd_min[tbd_elk_ung_short$Study_Area == "OK"]))
  (elk_meanNE <- mean(tbd_elk_ung_short$tbd_min[tbd_elk_ung_short$Study_Area == "NE"]))
  
  (mulie_meanOK <- mean(tbd_md_ung_short$tbd_min[tbd_md_ung_short$Monitoring == "Dirt road"]))
  (mulie_meanNE <- mean(tbd_md_ung_short$tbd_min[tbd_md_ung_short$Monitoring == "Trail"]))
  (wtd_meanOK <- mean(tbd_wtd_ung_short$tbd_min[tbd_wtd_ung_short$Monitoring == "Dirt road"]))
  (wtd_meanNE <- mean(tbd_wtd_ung_short$tbd_min[tbd_wtd_ung_short$Monitoring == "Trail"]))
  (moose_meanOK <- mean(tbd_moose_ung_short$tbd_min[tbd_moose_ung_short$Monitoring == "Dirt road"]))
  (moose_meanNE <- mean(tbd_moose_ung_short$tbd_min[tbd_moose_ung_short$Monitoring == "Trail"]))
  (elk_meanOK <- mean(tbd_elk_ung_short$tbd_min[tbd_elk_ung_short$Monitoring == "Dirt road"]))
  (elk_meanNE <- mean(tbd_elk_ung_short$tbd_min[tbd_elk_ung_short$Monitoring == "Trail"]))
  
  ####  Setup data & MCMC specifications for JAGS  ####
  #'  ----------------------------------------------
  #'  MCMC settings
  # nc <- 3; ni <- 100000; nb <- 75000; nt <- 1; na <- 20000
  nc <- 3; ni <- 75000; nb <- 25000; nt <- 1; na <- 5000
  # nc <- 3; ni <- 7500; nb <- 2000; nt <- 10; na <- 1000
  
  #'  Function to define and bundle data
  bundle_dat <- function(dat, ungulateID) {
    #'  Number of time-btwn-detection observations
    ntbd <- nrow(dat)
    #'  Number of unique camera locations
    ncams <- length(unique(dat$CameraLocation))
    #'  Format covariate data
    tbd_dat <- dplyr::select(dat, c(tbd_min, tbd_hour, tbd_day, CameraLocation, 
                                    Season, Complexity_index1, TRI, TRI_250m,
                                    PercForest, Species, PreviousDet, Study_Area, Monitoring)) %>%
      mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
             Season = as.numeric(factor(Season, levels = c("Summer", "Fall", "Winter", "Spring"))), # levels must be 1-4 (not 0-3) for nested indexing
             PreviousDet = as.numeric(factor(PreviousDet, levels = ungulateID)), # actual UngulateID per index varies based on which ungulate is the focal species
             Complexity_index1 = scale(Complexity_index1),
             TRI = scale(TRI),
             PercForest = scale(PercForest),
             Study_Area = ifelse(Study_Area == "NE", 0, 1),
             Monitoring = ifelse(Monitoring == "Dirt road", 0, 1))
    
    print(summary(tbd_dat))
    print(head(tbd_dat))
    
    #'  Covariate matrix for JAGS
    covs <- matrix(NA, ncol = 7, nrow = ntbd)
    covs[,1] <- tbd_dat$Season
    covs[,2] <- tbd_dat$PreviousDet
    covs[,3] <- tbd_dat$TRI
    covs[,4] <- tbd_dat$PercForest
    covs[,5] <- tbd_dat$Complexity_index1
    covs[,6] <- tbd_dat$Study_Area
    covs[,7] <- tbd_dat$Monitoring
    print(head(covs))
    
    #'  Generate range of continuous covariate values to predict across
    print(minmax_tri <- range(covs[,2])); print(minmax_for <- range(covs[,3]))
    newTRI <- seq(from = minmax_tri[1], to = minmax_tri[2], length.out = 100)
    newFor <- seq(from = minmax_for[1], to = minmax_for[2], length.out = 100)
    newcovs <- as.matrix(cbind(newTRI, newFor))
    
    #'  Number of covariates
    ncovs <- ncol(covs)
    
    #'  Time-between-detections
    tbd <- dat$tbd_min
    print(summary(tbd))
    hist(tbd)
    
    bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                    site = tbd_dat$cams, newcovs = newcovs)
    return(bundled)
  }
  #'  Provide specific order for ungulateID levels - will differ for each species
  #'  Order generally goes mule deer, white-tailed deer, elk, moose
  md_ung_bundled <- bundle_dat(tbd_md_ung_short, ungulateID = c("White-tailed Deer", "Elk", "Moose"))
  tbd_elk_ung_short <- filter(tbd_elk_ung_short, Study_Area == "NE")
  elk_ung_bundled <- bundle_dat(tbd_elk_ung_short, ungulateID = c("Mule Deer", "White-tailed Deer", "Moose"))
  moose_ung_bundled <- bundle_dat(tbd_moose_ung_short, ungulateID = c("Mule Deer", "White-tailed Deer", "Elk"))
  wtd_ung_bundled <- bundle_dat(tbd_wtd_ung_short, ungulateID = c("Mule Deer", "Elk", "Moose"))
  
  #'  Save for making figures
  ung_md_bundled <- md_ung_bundled; save(ung_md_bundled, file = "./Data/ung_md_bundled.RData")
  ung_elk_bundled <- elk_ung_bundled; save(ung_elk_bundled, file = "./Data/ung_elk_bundled.RData")
  ung_moose_bundled <- moose_ung_bundled; save(ung_moose_bundled, file = "./Data/ung_moose_bundled.RData")
  ung_wtd_bundled <- wtd_ung_bundled; save(ung_wtd_bundled, file = "./Data/ung_wtd_bundled.RData")
  
  
  #'  ------------------------
  #####  MULE DEER Analysis  ####
  #'  ------------------------
  #'  Source JAGS model
  #'  Make sure inits and parameters being monitored match up with sourced model
  #'  Make sure model parameterization matches order of covariates in bundled data
  source("./Scripts/JAGS_models/JAGS_tbdungulates_season_ungID_X_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(md_ung_bundled$y, list(md_ung_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "beta3", "sigma", "season.tbd", 
             "ung.tbd", "mu.tbd")   #"sa.tbd", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.ung.md <- jags(md_ung_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_X_habitat.txt',
                     inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                     n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.ung.md)
  mcmcplot(tbd.ung.md$samples)
  save(tbd.ung.md, file = "./Outputs/TimeBtwnDetections/tbd.ung.md-season_ungID_X_habitat.RData")
  #'  Keep in mind UngulateID levels are white-tailed deer [1], elk [2], moose [3]
  
  
  #'  ------------------
  #####  ELK Analysis  ####
  #'  ------------------
  #'  Source JAGS model
  #'  NOTE: dropping the interaction term b/c model failed to converge with it
  #'  Still isn't amazing - lots of autocorrelation in some variables but probably OK for now
  source("./Scripts/JAGS_models/JAGS_tbdungulates_season_ungID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(elk_ung_bundled$y, list(elk_ung_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", 
              "ung.tbd", "mu.tbd")   #"sa.tbd", "beta3", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.ung.elk <- jags(elk_ung_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat.txt', #tbd_season_ungID_X_habitat
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.ung.elk)
  mcmcplot(tbd.ung.elk$samples)
  save(tbd.ung.elk, file = "./Outputs/TimeBtwnDetections/tbd.ung.elk-season_ungID_habitat.RData") #tbd.ung.elk-season_ungID_X_habitat
  #'  Keep in mind UngulateID levels are mule deer [1], white-tailed deer [2], moose [3]
  
  
  #'  --------------------
  #####  MOOSE Analysis  ####
  #'  --------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdungulates_season_ungID_X_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(moose_ung_bundled$y, list(moose_ung_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "beta3", "sigma", "season.tbd", 
              "ung.tbd", "mu.tbd")   #"sa.tbd", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.ung.moose <- jags(moose_ung_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_X_habitat.txt',
                    inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                    n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.ung.moose)
  mcmcplot(tbd.ung.moose$samples)
  save(tbd.ung.moose, file = "./Outputs/TimeBtwnDetections/tbd.ung.moose-season_ungID_X_habitat.RData")
  #'  Keep in mind UngulateID levels are mule deer [1], white-tailed deer [2], elk [3]
  
  
  #'  --------------------------------
  #####  WHITE-TAILED DEER Analysis  ####
  #'  --------------------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdungulates_season_ungID_X_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(wtd_ung_bundled$y, list(wtd_ung_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "beta3", "sigma", "season.tbd", 
              "ung.tbd", "mu.tbd")   #"sa.tbd",   
  
  #'  Run model
  start.time <- Sys.time()
  tbd.ung.wtd <- jags(wtd_ung_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_X_habitat.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.ung.wtd)
  mcmcplot(tbd.ung.wtd$samples)
  save(tbd.ung.wtd, file = "./Outputs/TimeBtwnDetections/tbd.ung.wtd-season_ungID_X_habitat.RData")
  #'  Keep in mind UngulateID levels are mule deer [1], elk [2], moose [3]
  
  
  
  

  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT  ####
  
  