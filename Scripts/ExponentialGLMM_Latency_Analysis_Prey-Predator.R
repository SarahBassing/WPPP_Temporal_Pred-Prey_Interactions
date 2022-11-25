  #'  -------------------------------------------------
  #'  Prey effects on latency of site use for predators
  #'  WA Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Nov. 2022
  #'  -------------------------------------------------
  #'  Read in time-between-detection data (times between detection of a prey species
  #'  followed by a predator) and run generalized linear model with an 
  #'  exponential distribution in JAGS to estimate effect of recent prey presence 
  #'  on timing of predator use of camera sites following detection of an ungulate. 
  #'  
  #'  Time-between-detections data generated in TimeBetweenDetections.R script
  #'  Source individual JAGS models depending on predator species.
  #'  ---------------------------------------------
  
  #'  Clear workspace & load libraries
  rm(list = ls())
  
  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  
  #'  Read in data
  load("./Outputs/tbd_prey.pred_2022-11-23.RData") 
  #'  Remove observations involving lynx due to too few lynx detections
  tbd_prey.pred <- tbd_prey.pred %>%
    #'  Change units of time
    mutate(tbd_min_round = round(TimeSinceLastDet, 0),
           tbd_min = TimeSinceLastDet,
           tbd_hour = tbd_min/60,
           tbd_day = tbd_hour/24) %>%
    dplyr::select(-c(TimeSinceLastDet, HuntingMode, TrophicLevel)) %>%
    rename(UngulateID = PredatorID)
  
  #'  Filter images by ungulate species
  tbd_bear <- filter(tbd_prey.pred, Species == "Black Bear")
  tbd_bob <- filter(tbd_prey.pred, Species == "Bobcat")
  tbd_coug <- filter(tbd_prey.pred, Species == "Cougar")
  tbd_coy <- filter(tbd_prey.pred, Species == "Coyote")
  tbd_wolf <- filter(tbd_prey.pred, Species == "Wolf")
  
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
    
    #'  Summary of observations with each ungulate species
    print("Total TBDs with each ungulate species")
    print(table(short_tbd$UngulateID))
    #'  Summary of observations in each season
    print("Total TBDs for each season")
    print(table(short_tbd$Season))
    
    #'  Return dataset after removing extreme values
    return(short_tbd)
  }
  tbd_bear_short <- tbd_summary(tbd_bear, spp = "black bear", quant = 0.97)
  tbd_bob_short <- tbd_summary(tbd_bob, spp = "bobcat", quant = 0.97)
  tbd_coug_short <- tbd_summary(tbd_coug, spp = "cougar", quant = 0.97)
  tbd_coy_short <- tbd_summary(tbd_coy, spp = "coyote", quant = 0.97)
  tbd_wolf_short <- tbd_summary(tbd_wolf, spp = "wolf", quant = 0.95)
  
  tbd_all_short <- rbind(tbd_bear_short, tbd_bob_short, tbd_coug_short, tbd_coy_short, tbd_wolf_short)
  # write.csv(tbd_all_short, "./Outputs/tbd_prey.pred_NoOutliers.csv")
  
  #'  Are there study area differences for the predator mean tbd?
  (bear_meanOK <- mean(tbd_bear_short$tbd_min[tbd_bear_short$Study_Area == "OK"]))
  (bear_meanNE <- mean(tbd_bear_short$tbd_min[tbd_bear_short$Study_Area == "NE"]))
  (bob_meanOK <- mean(tbd_bob_short$tbd_min[tbd_bob_short$Study_Area == "OK"]))
  (bob_meanNE <- mean(tbd_bob_short$tbd_min[tbd_bob_short$Study_Area == "NE"]))
  (coug_meanOK <- mean(tbd_coug_short$tbd_min[tbd_coug_short$Study_Area == "OK"]))
  (coug_meanNE <- mean(tbd_coug_short$tbd_min[tbd_coug_short$Study_Area == "NE"]))
  (coy_meanOK <- mean(tbd_coy_short$tbd_min[tbd_coy_short$Study_Area == "OK"]))
  (coy_meanNE <- mean(tbd_coy_short$tbd_min[tbd_coy_short$Study_Area == "NE"]))
  (wolf_meanOK <- mean(tbd_wolf_short$tbd_min[tbd_wolf_short$Study_Area == "OK"]))
  (wolf_meanNE <- mean(tbd_wolf_short$tbd_min[tbd_wolf_short$Study_Area == "NE"]))
  #'  Are there differences based on trail type?
  (bear_meanRd <- mean(tbd_bear_short$tbd_min[tbd_bear_short$Monitoring == "Dirt road"]))
  (bear_meanTrl <- mean(tbd_bear_short$tbd_min[tbd_bear_short$Monitoring == "Trail"]))
  (bob_meanRd <- mean(tbd_bob_short$tbd_min[tbd_bob_short$Monitoring == "Dirt road"]))
  (bob_meanTrl <- mean(tbd_bob_short$tbd_min[tbd_bob_short$Monitoring == "Trail"]))
  (coug_meanRd <- mean(tbd_coug_short$tbd_min[tbd_coug_short$Monitoring == "Dirt road"]))
  (coug_meanTrl <- mean(tbd_coug_short$tbd_min[tbd_coug_short$Monitoring == "Trail"]))
  (coy_meanRd <- mean(tbd_coy_short$tbd_min[tbd_coy_short$Monitoring == "Dirt road"]))
  (coy_meanTrl <- mean(tbd_coy_short$tbd_min[tbd_coy_short$Monitoring == "Trail"]))
  (wolf_meanRd <- mean(tbd_wolf_short$tbd_min[tbd_wolf_short$Monitoring == "Dirt road"]))
  (wolf_meanTrl <- mean(tbd_wolf_short$tbd_min[tbd_wolf_short$Monitoring == "Trail"]))
  
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
                                    Season, UngulateID, Complexity_index1, 
                                    backgroundRisk_HCI, backgroundRisk_TRI, 
                                    backgroundRisk_For, TRI, TRI_250m, PercForest, 
                                    Monitoring, Species, spp_pair, Study_Area)) %>%
      mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
             Season = as.numeric(factor(Season, levels = c("Summer", "Fall", "Winter", "Spring"))), # levels must be 1-4 (not 0-3) for nested indexing
             UngulateID = as.numeric(factor(UngulateID, levels = c("Elk", "Moose", "Mule Deer", "White-tailed Deer"))), # levels must be 1-4 for nested indexing
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
    covs[,2] <- tbd_dat$UngulateID
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
  #'  Remove the two winter bear observation
  tbd_bear_shorter <- tbd_bear_short %>% filter(Season != "Winter") 
  bear_bundled <- bundle_dat(tbd_bear_shorter)
  #'  Convert seasonal categories in covs matrix from 1-4 to 1-3 since the few 
  #'  winter [3] observations were excluded
  bear_bundled[[2]][,1][bear_bundled[[2]][,1] == 4] <- 3 #' Spring now season 3
  bob_bundled <- bundle_dat(tbd_bob_short) 
  coug_bundled <- bundle_dat(tbd_coug_short)
  coy_bundled <- bundle_dat(tbd_coy_short) 
  #'  Remove the few spring, elk, and moose observations from the wolf data
  tbd_wolf_shorter <- tbd_wolf_short %>% filter(Season != "Spring") %>%
    filter(UngulateID == "Mule Deer" | UngulateID == "White-tailed Deer")
  wolf_bundled <- bundle_dat(tbd_wolf_shorter)
  #'  Convert ungulateID categories in covs matrix from 1-4 to 1-2 since elk [1] 
  #'  and moose [2] observations were excluded (don't need to worry about spring [4] 
  #'  b/c already the last (4th) category so just drops off)
  wolf_bundled[[2]][,2][wolf_bundled[[2]][,2] == 3] <- 1 #' Mule deer now ung 1
  wolf_bundled[[2]][,2][wolf_bundled[[2]][,2] == 4] <- 2 #' White-tailed deer now ung 2
  
  #'  Save for making figures
  prey_bear_bundled <- bear_bundled; save(prey_bear_bundled, file = "./Data/prey_bear_bundled.RData")
  prey_bob_bundled <- bob_bundled; save(prey_bob_bundled, file = "./Data/prey_bob_bundled.RData")
  prey_coug_bundled <- coug_bundled; save(prey_coug_bundled, file = "./Data/prey_coug_bundled.RData")
  prey_coy_bundled <- coy_bundled; save(prey_coy_bundled, file = "./Data/prey_coy_bundled.RData")
  prey_wolf_bundled <- wolf_bundled; save(prey_wolf_bundled, file = "./Data/prey_wolf_bundled.RData")
  
  #'  ==========================
  ####  No Interaction Models  ####
  #'  ==========================
  #'  Including interaction terms resulted in over-parameterization of models and
  #'  poor convergence for most of the predator species. Chose to exclude the 
  #'  interaction across all models.
  
  #'  -----------------------------------
  #####  Prey - BLACK BEAR Analysis  ####
  #'  -----------------------------------
  #'  Source JAGS model
  #'  Note: winter observations excluded so model is slightly different
  source("./Scripts/JAGS_models/JAGS_tbdpreypred_season_ungID_habitat_noWtr.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(bear_bundled$y, list(bear_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "prey.tbd", "mu.tbd") 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.prey.bear <- jags(bear_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat_noWtr.txt',
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.prey.bear)
  mcmcplot(tbd.prey.bear$samples)
  save(tbd.prey.bear, file = "./Outputs/TimeBtwnDetections/tbd.prey.bear-season_ungID_habitat.RData")
  
  
  #'  -----------------------------
  #####  Prey - BOBCAT Analysis  ####
  #'  -----------------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdpreypred_season_ungID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(bob_bundled$y, list(bob_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "prey.tbd", "mu.tbd") 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.prey.bob <- jags(bob_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat.txt',
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.prey.bob)
  mcmcplot(tbd.prey.bob$samples)   
  save(tbd.prey.bob, file = "./Outputs/TimeBtwnDetections/tbd.prey.bob-season_ungID_habitat.RData")
  
  
  #'  -------------------------------
  #####  Prey - COUGAR Analysis  ####
  #'  -------------------------------
  #'  Source JAGS model
  #'  Note: Dropped the interaction - overparameterized & poor convergence with interaction
  source("./Scripts/JAGS_models/JAGS_tbdpreypred_season_ungID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(coug_bundled$y, list(coug_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "prey.tbd", "mu.tbd")   
  
  #'  Run model
  start.time <- Sys.time()
  tbd.prey.coug <- jags(coug_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat.txt',
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                         n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.prey.coug)
  mcmcplot(tbd.prey.coug$samples)  
  save(tbd.prey.coug, file = "./Outputs/TimeBtwnDetections/tbd.prey.coug-season_ungID_habitat.RData")
  
  
  #'  -------------------------------------------
  #####  Prey - COYOTE Analysis  ####
  #'  -------------------------------------------
  #'  Source JAGS model
  #'  Note: excluding the interaction term - converges better without it
  source("./Scripts/JAGS_models/JAGS_tbdpreypred_season_ungID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(coy_bundled$y, list(coy_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "prey.tbd", "mu.tbd") 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.prey.coy <- jags(coy_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat.txt',
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.prey.coy)
  mcmcplot(tbd.prey.coy$samples)
  save(tbd.prey.coy, file = "./Outputs/TimeBtwnDetections/tbd.prey.coy-season_ungID_habitat.RData")
  
  
  
  #'  -------------------------------------------
  #####  Prey - WOLF Analysis  ####
  #'  -------------------------------------------
  #'  Source JAGS model
  #'  Note: removes spring, elk, and moose observations so model is slightly different
  source("./Scripts/JAGS_models/JAGS_tbdpreypred_season_ungID_habitat_noSprgElkMoose.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(wolf_bundled$y, list(wolf_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(1,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "prey.tbd", "mu.tbd") 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.prey.wolf <- jags(wolf_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat_noSprgElkMoose.txt',
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.prey.wolf)
  mcmcplot(tbd.prey.wolf$samples)
  save(tbd.prey.wolf, file = "./Outputs/TimeBtwnDetections/tbd.prey.wolf-season_ungID_habitat.RData")
  
  
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT - X^2 test  ####
  
  
