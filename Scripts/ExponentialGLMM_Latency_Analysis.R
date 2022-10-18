  #'  ---------------------------------------------
  #'  Predation risk effects on latency of site use
  #'  WA Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Sept. 2022
  #'  ---------------------------------------------
  #'  Read in time-between-detection data (times between detection of a predator
  #'  followed by a prey species) and run generalized linear model with an 
  #'  exponential distribution in JAGS to estimate effect of predation risk on 
  #'  latency of prey use of camera sites following detection of a predator. 
  #'  Predation risk defined as "risky habitat" (e.g., habitat complexity) 
  #'  associated with different predator hunting modes at a camera site.
  #'  
  #'  Handy website for Poisson and exponential lambda interpretation:
  #'  https://towardsdatascience.com/what-is-exponential-distribution-7bdd08590e2a#:~:text=The%20definition%20of%20exponential%20distribution,Poisson%20(X%3D0)
  #'  The definition of exponential distribution is the probability distribution 
  #'  of the time *between* the events in a Poisson process.
  #'  Poisson distribution assumes that events occur independent of one another.
  #'  y ~ Exp(lambda) --> lambda is the reciprocal (1/λ) of the rate (λ) in a Poisson distribution.
  #'  lambda is often expressed in the amount of time before an event occurs.
  #'  
  #'  Time-between-detections data generated in TimeBetweenDetections.R script
  #'  Source individual JAGS models depending on ungulate species.
  #'  ---------------------------------------------
  
  #'  Clear workspace & load libraries
  rm(list = ls())
  
  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  
  #'  Read in data
  load("./Outputs/tbd_pred.prey_2022-10-05.RData") #2022-09-23
  #'  Remove observations involving lynx due to too few lynx detections
  tbd_pred.prey <- tbd_pred.prey %>%
    filter(PredatorID != "Lynx") %>%
    #'  Change units of time
    mutate(tbd_min_round = round(TimeSinceLastDet, 0),
           tbd_min = TimeSinceLastDet,
           tbd_hour = tbd_min/60,
           tbd_day = tbd_hour/24) %>%
    dplyr::select(-TimeSinceLastDet)
  
  #'  Filter images by ungulate species
  tbd_md <- filter(tbd_pred.prey, Species == "Mule Deer")
  tbd_elk <- filter(tbd_pred.prey, Species == "Elk")
  tbd_moose <- filter(tbd_pred.prey, Species == "Moose")
  tbd_wtd <- filter(tbd_pred.prey, Species == "White-tailed Deer")
  
  #'  Function to identify time cutoff after which tbd observations are outliers
  tbd_summary <- function(tbd, spp, quant) {
    #'  Plot frequency of time-between-detections (should look exponential)
    hist(tbd$tbd_day, breaks = 50, main = paste("Number of days between detections for\n", spp))
    boxplot(tbd$tbd_day, ylab = "Days", main = paste("Number of days between detections for\n", spp))
    
    #'  Review range of TBD values
    print("Quantiles of days between detections of predator followed by prey")
    print(quantile(tbd$tbd_day))
    #'  Review 90 - 99th quartiles- where are there gaps and outliers in the distribution?
    print("90th, 95th, and 99th quantiles of days between detections of predator followed by prey") 
    print(quantile(tbd$tbd_day, c(0.9, 0.95, 0.97, 0.99)))
    
    #'  Re-plot frequency of time-btwn-detections after removing extreme values
    short_tbd <- filter(tbd, tbd_day <= quantile(tbd$tbd_day, c(quant)))
    hist(short_tbd$tbd_day, breaks = 25, main = paste("Number of days between detections for\n", spp, "up to quantile =", quant))
    boxplot(short_tbd$tbd_day, ylab = "Days", main = paste("Number of days between detections for\n", spp, "up to quantile =", quant))
    
    #'  Summary of observations with each predator species
    print("Total TBDs with each predator species")
    print(table(short_tbd$PredatorID))
    #'  Summary of observations in each season
    print("Total TBDs for each season")
    print(table(short_tbd$Season))
    
    #'  Return dataset after removing extreme values
    return(short_tbd)
  }
  tbd_md_short <- tbd_summary(tbd_md, spp = "mule deer", quant = 0.97)
  tbd_elk_short <- tbd_summary(tbd_elk, spp = "elk", quant = 1.0)
  tbd_moose_short <- tbd_summary(tbd_moose, spp = "moose", quant = 0.97)
  tbd_wtd_short <- tbd_summary(tbd_wtd, spp = "white-tailed deer", quant = 0.97)
  
  tbd_all_short <- rbind(tbd_md_short, tbd_elk_short, tbd_moose_short, tbd_wtd_short)
  # write.csv(tbd_all_short, "./Outputs/tbd_pred.prey_NoOutliers.csv")
  
  #'  Are there study area differences for the ungulate mean tbd?
  (mulie_meanOK <- mean(tbd_md_short$tbd_min[tbd_md_short$Study_Area == "OK"]))
  (mulie_meanNE <- mean(tbd_md_short$tbd_min[tbd_md_short$Study_Area == "NE"]))
  (wtd_meanOK <- mean(tbd_wtd_short$tbd_min[tbd_wtd_short$Study_Area == "OK"]))
  (wtd_meanNE <- mean(tbd_wtd_short$tbd_min[tbd_wtd_short$Study_Area == "NE"]))
  (moose_meanOK <- mean(tbd_moose_short$tbd_min[tbd_moose_short$Study_Area == "OK"]))
  (moose_meanNE <- mean(tbd_moose_short$tbd_min[tbd_moose_short$Study_Area == "NE"]))
  (elk_meanOK <- mean(tbd_elk_short$tbd_min[tbd_elk_short$Study_Area == "OK"]))
  (elk_meanNE <- mean(tbd_elk_short$tbd_min[tbd_elk_short$Study_Area == "NE"]))
  
  (mulie_meanOK <- mean(tbd_md_short$tbd_min[tbd_md_short$Monitoring == "Dirt road"]))
  (mulie_meanNE <- mean(tbd_md_short$tbd_min[tbd_md_short$Monitoring == "Trail"]))
  (wtd_meanOK <- mean(tbd_wtd_short$tbd_min[tbd_wtd_short$Monitoring == "Dirt road"]))
  (wtd_meanNE <- mean(tbd_wtd_short$tbd_min[tbd_wtd_short$Monitoring == "Trail"]))
  (moose_meanOK <- mean(tbd_moose_short$tbd_min[tbd_moose_short$Monitoring == "Dirt road"]))
  (moose_meanNE <- mean(tbd_moose_short$tbd_min[tbd_moose_short$Monitoring == "Trail"]))
  (elk_meanOK <- mean(tbd_elk_short$tbd_min[tbd_elk_short$Monitoring == "Dirt road"]))
  (elk_meanNE <- mean(tbd_elk_short$tbd_min[tbd_elk_short$Monitoring == "Trail"]))
  
  ####  Setup data & MCMC specifications for JAGS  ####
  #'  ----------------------------------------------
  #'  MCMC settings
  # nc <- 3; ni <- 100000; nb <- 75000; nt <- 10; na <- 20000
  nc <- 3; ni <- 75000; nb <- 25000; nt <- 1; na <- 5000
  # nc <- 3; ni <- 7500; nb <- 2000; nt <- 10; na <- 1000

  #'  Function to define and bundle data
  bundle_dat <- function(dat) { 
    #'  Number of time-btwn-detection observations
    ntbd <- nrow(dat)
    #'  Number of unique camera locations
    ncams <- length(unique(dat$CameraLocation))
    #'  Format covariate data
    tbd_dat <- dplyr::select(dat, c(tbd_min, tbd_hour, tbd_day, CameraLocation, 
                                    Season, PredatorID, HuntingMode, TrophicLevel,
                                    Complexity_index1, backgroundRisk_HCI,
                                    backgroundRisk_TRI, backgroundRisk_For, 
                                    TRI, TRI_250m, PercForest, Monitoring, 
                                    Species, spp_pair, Study_Area)) %>%
      mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
             Season = as.numeric(factor(Season, levels = c("Summer", "Fall", "Winter", "Spring"))), # levels must be 1-4 (not 0-3) for nested indexing
             PredatorID = as.numeric(factor(PredatorID, levels = c("Bobcat", "Coyote", "Black Bear", "Cougar", "Wolf"))), # levels must be 1-5 for nested indexing
             Complexity_index1 = scale(Complexity_index1),
             TRI = scale(TRI),
             PercForest = scale(PercForest),
             Monitoring = ifelse(Monitoring == "Dirt road", 0, 1),
             TRI_level = ifelse(backgroundRisk_TRI == "Low", 0, 1),
             PercFor_level = ifelse(backgroundRisk_For == "Low", 0, 1),
             Study_Area = ifelse(Study_Area == "NE", 0, 1))
    print(summary(tbd_dat))
    print(head(tbd_dat))
    
    #'  Covariate matrix for JAGS
    covs <- matrix(NA, ncol = 9, nrow = ntbd)
    covs[,1] <- tbd_dat$Season
    covs[,2] <- tbd_dat$PredatorID
    covs[,3] <- tbd_dat$Complexity_index1
    covs[,4] <- tbd_dat$TRI
    covs[,5] <- tbd_dat$PercForest
    covs[,6] <- tbd_dat$Monitoring
    covs[,7] <- tbd_dat$TRI_level
    covs[,8] <- tbd_dat$PercFor_level
    covs[,9] <- tbd_dat$Study_Area
    head(covs)
    
    #'  Generate range of continuous covariate values to predict across
    print(minmax_tri <- range(covs[,4])); print(minmax_for <- range(covs[,5]))
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
  md_bundled <- bundle_dat(tbd_md_short) 
  wtd_bundled <- bundle_dat(tbd_wtd_short)
  #'  Remove single wolf-elk observation
  tbd_elk_shorter <- tbd_elk_short %>% filter(PredatorID != "Wolf")
  elk_bundled <- bundle_dat(tbd_elk_shorter)
  moose_bundled <- bundle_dat(tbd_moose_short)
  all_bundled <- bundle_dat(tbd_all_short)
  
  #'  Save for making figures
  pred_md_bundled <- md_bundled; save(pred_md_bundled, file = "./Data/pred_md_bundled.RData")
  pred_elk_bundled <- md_bundled; save(pred_elk_bundled, file = "./Data/pred_elk_bundled.RData")
  pred_moose_bundled <- md_bundled; save(pred_moose_bundled, file = "./Data/pred_moose_bundled.RData")
  pred_wtd_bundled <- md_bundled; save(pred_wtd_bundled, file = "./Data/pred_wtd_bundled.RData")
  
  #'  ==========================
  ####  No Interaction Models  ####
  #'  ==========================
  
  #'  -----------------------------------
  #####  Predator - MULE DEER Analysis  ####
  #'  -----------------------------------
  #'  Source JAGS model
  #'  Make sure inits and parameters being monitored match up with sourced model
  #'  Make sure model parameterization matches order of covariates in bundled data
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(md_bundled$y, list(md_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", 
              "mu.tbd", "pred.tbd.tri", "pred.tbd.for") 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.md <- jags(md_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.md)
  mcmcplot(tbd.pred.md$samples[,1:34])
  save(tbd.pred.md, file = "./Outputs/TimeBtwnDetections/tbd.pred.md-season_predID_habitat.RData")
  
  
  #'  -----------------------------
  #####  Predator - ELK Analysis  ####
  #'  -----------------------------
  #'  Source JAGS model
  #'  NOTE: NO wolf observations in this model so slightly different version of model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat_nowolf.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(elk_bundled$y, list(elk_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", 
              "mu.tbd")   #"sa.tbd", , "pred.tbd.tri", "pred.tbd.for"
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.elk <- jags(elk_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat_nowolf.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.elk)
  mcmcplot(tbd.pred.elk$samples)  
  save(tbd.pred.elk, file = "./Outputs/TimeBtwnDetections/tbd.pred.elk-season_predID_habitat.RData")
  
  
  #'  -------------------------------
  #####  Predator - MOOSE Analysis  ####
  #'  -------------------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(moose_bundled$y, list(moose_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", 
              "mu.tbd", "pred.tbd.tri", "pred.tbd.for")   
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.moose <- jags(moose_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.moose)
  mcmcplot(tbd.pred.moose$samples[,1:34])  
  save(tbd.pred.moose, file = "./Outputs/TimeBtwnDetections/tbd.pred.moose-season_predID_habitat.RData")
  
  
  #'  -------------------------------------------
  #####  Predator - WHITE-TAILED DEER Analysis  ####
  #'  -------------------------------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(wtd_bundled$y, list(wtd_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", 
              "mu.tbd", "pred.tbd.tri", "pred.tbd.for")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.wtd <- jags(wtd_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.wtd)
  mcmcplot(tbd.pred.wtd$samples[,1:34])
  save(tbd.pred.wtd, file = "./Outputs/TimeBtwnDetections/tbd.pred.wtd-season_predID_habitat.RData")
  
  
  #'  -----------------------------------
  #####  Predator - ALL PREY Analysis  ####
  #'  -----------------------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(all_bundled$y, list(all_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", 
              "mu.tbd", "pred.tbd.tri", "pred.tbd.for") 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.all <- jags(all_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat.txt',
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.all)
  mcmcplot(tbd.pred.all$samples[,1:34])
  save(tbd.pred.all, file = "./Outputs/TimeBtwnDetections/tbd.pred.all-season_predID_habitat.RData")
  
  
  
  #'  ========================================
  ####  Predator-Habitat Interaction Models  ####
  #'  ========================================
  
  #'  -----------------------------------
  #####  Predator - MULE DEER Analysis  ####
  #'  -----------------------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_X_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(md_bundled$y, list(md_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "beta3", "sigma", 
              "season.tbd", "pred.tbd", "mu.tbd") #"beta4", "sa.tbd",  , "pred.tbd.tri", "pred.tbd.for"
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.md <- jags(md_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_X_habitat.txt',
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.md)
  mcmcplot(tbd.pred.md$samples)
  save(tbd.pred.md, file = "./Outputs/TimeBtwnDetections/tbd.pred.md-season_predID_X_habitat.RData")
  
  
  #'  -----------------------------
  #####  Predator - ELK Analysis  ####
  #'  -----------------------------
  #'  Source JAGS model
  #'  NOTE: NO wolf observations in this model so slightly different version of model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_X_habitat_nowolf.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(elk_bundled$y, list(elk_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "beta3", "sigma", 
              "season.tbd", "pred.tbd", "mu.tbd")  #"beta4", "sa.tbd", , "pred.tbd.tri", "pred.tbd.for"
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.elk <- jags(elk_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_X_habitat_nowolf.txt',
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.elk)
  mcmcplot(tbd.pred.elk$samples)   # Rhat & traceplots indicate model isn't converging well even with 150,000 iterations - revert to model w/o interactions
  save(tbd.pred.elk, file = "./Outputs/TimeBtwnDetections/tbd.pred.elk-season_predID_X_habitat.RData")
  
  
  #'  -------------------------------
  #####  Predator - MOOSE Analysis  ####
  #'  -------------------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_X_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(moose_bundled$y, list(moose_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "beta3", "sigma", 
              "season.tbd", "pred.tbd", "mu.tbd")   #"beta4", "sa.tbd", , "pred.tbd.tri", "pred.tbd.for"
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.moose <- jags(moose_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_X_habitat.txt',
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                         n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.moose)
  mcmcplot(tbd.pred.moose$samples)  #'  Rhat & traceplots look good compared to model w/o interactions - use this one
  save(tbd.pred.moose, file = "./Outputs/TimeBtwnDetections/tbd.pred.moose-season_predID_X_habitat.RData")
  
  
  #'  -------------------------------------------
  #####  Predator - WHITE-TAILED DEER Analysis  ####
  #'  -------------------------------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_X_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(wtd_bundled$y, list(wtd_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "beta3", "sigma", 
              "season.tbd", "pred.tbd", "mu.tbd") #"beta4", "sa.tbd", , "pred.tbd.tri", "pred.tbd.for"
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.wtd <- jags(wtd_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_X_habitat.txt',
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.wtd)
  mcmcplot(tbd.pred.wtd$samples)
  save(tbd.pred.wtd, file = "./Outputs/TimeBtwnDetections/tbd.pred.wtd-season_predID_X_habitat.RData")
  
  
  
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT - X^2 test  ####
  
  
  