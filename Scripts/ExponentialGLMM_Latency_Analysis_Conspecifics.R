  #'  ------------------------------------------
  #'  Conspecific effects on latency of site use
  #'  WA Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Sept. 2022
  #'  ------------------------------------------
  #'  Read in time-between-detection data (times between detection of an ungulate
  #'  followed by a conspecific) and build generalized linear model with an 
  #'  exponential distribution in JAGS to estimate latency of camera site use 
  #'  following detection of the same species. To be compared to mean latency of
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
  load("./Outputs/tbd_conspif_2022-09-20.RData") 
  #'  Remove observations involving lynx due to too few lynx detections
  tbd_conspif <- tbd_conspif %>%
    #'  Change units of time
    mutate(tbd_min_round = round(TimeSinceLastDet, 0),
           tbd_min = TimeSinceLastDet,
           tbd_hour = tbd_min/60,
           tbd_day = tbd_hour/24) %>%
    dplyr::select(-TimeSinceLastDet)
  
  #'  Filter images by ungulate species
  tbd_md_con <- filter(tbd_conspif, Species == "Mule Deer")
  tbd_elk_con <- filter(tbd_conspif, Species == "Elk")
  tbd_moose_con <- filter(tbd_conspif, Species == "Moose")
  tbd_wtd_con <- filter(tbd_conspif, Species == "White-tailed Deer")
  
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
  tbd_md_con_short <- tbd_summary(tbd_md_con, spp = "mule deer", quant = 0.97)
  tbd_elk_con_short <- tbd_summary(tbd_elk_con, spp = "elk", quant = 0.97)
  tbd_moose_con_short <- tbd_summary(tbd_moose_con, spp = "moose", quant = 0.95)
  tbd_wtd_con_short <- tbd_summary(tbd_wtd_con, spp = "white-tailed deer", quant = 0.97)
  
  tbd_all_con_short <- rbind(tbd_md_con_short, tbd_elk_con_short, tbd_moose_con_short, tbd_wtd_con_short)
  
  ####  Setup data & MCMC specifications for JAGS  ####
  #'  ----------------------------------------------
  #'  MCMC settings
  nc <- 3; ni <- 50000; nb <- 10000; nt <- 5; na <- 2000
  
  #'  Function to define and bundle data
  bundle_dat <- function(dat) {
    #'  Number of time-btwn-detection observations
    ntbd <- nrow(dat)
    #'  Number of unique camera locations
    ncams <- length(unique(dat$CameraLocation))
    #'  Format covariate data
    tbd_dat <- dplyr::select(dat, c(tbd_min, tbd_hour, tbd_day, CameraLocation, 
                                    Season, Complexity_index1, TRI, TRI_250m,
                                    PercForest, Species)) %>%
      mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
             Season = as.numeric(factor(Season, levels = c("Summer", "Fall", "Winter", "Spring"))), # levels must be 1-4 (not 0-3) for nested indexing
             Complexity_index1 = scale(Complexity_index1),
             TRI = scale(TRI),
             PercForest = scale(PercForest))
    
    print(summary(tbd_dat))
    print(head(tbd_dat))
    
    #'  Covariate matrix for JAGS
    covs <- matrix(NA, ncol = 4, nrow = ntbd)
    covs[,1] <- tbd_dat$Season
    covs[,2] <- tbd_dat$TRI
    covs[,3] <- tbd_dat$PercForest
    covs[,4] <- tbd_dat$Complexity_index1
    print(head(covs))
    
    #'  Number of covariates
    ncovs <- ncol(covs)
    
    #'  Time-between-detections
    tbd <- dat$tbd_min
    print(summary(tbd))
    hist(tbd)
    
    bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                    site = tbd_dat$cams)
    return(bundled)
  }
  md_con_bundled <- bundle_dat(tbd_md_con_short)
  elk_con_bundled <- bundle_dat(tbd_elk_con_short)
  moose_con_bundled <- bundle_dat(tbd_moose_con_short)
  wtd_con_bundled <- bundle_dat(tbd_wtd_con_short)
  
  
  #'  ------------------------
  #####  MULE DEER Analysis  ####
  #'  ------------------------
  #'  Source JAGS model
  #'  Make sure inits and parameters being monitored match up with sourced model
  #'  Make sure model parameterization matches order of covariates in bundled data
  source("./Scripts/JAGS_models/JAGS_tbdconspecific_season_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(md_con_bundled$y, list(md_con_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "sigma", "mu.mu") #"beta3", "beta4", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.md <- jags(md_con_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_habitat.txt',
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.md)
  mcmcplot(tbd.md$samples)
  save(tbd.md, file = "./Outputs/TimeBtwnDetections/tbd.md-season_habitat.RData")
  
  
  #'  ------------------
  #####  ELK Analysis  ####
  #'  ------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdconspecific_season_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(elk_con_bundled$y, list(elk_con_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "sigma", "mu.mu") #"beta3", "beta4", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.elk <- jags(elk_con_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_habitat.txt',
                 inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                 n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.elk)
  mcmcplot(tbd.elk$samples)
  save(tbd.elk, file = "./Outputs/TimeBtwnDetections/tbd.elk-season_habitat.RData")
  
  
  #'  --------------------
  #####  MOOSE Analysis  ####
  #'  --------------------
  #'  Source JAGS model
  source("./Scripts/JAGS_models/JAGS_tbdconspecific_season_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(moose_con_bundled$y, list(moose_con_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "sigma", "mu.mu") #"beta3", "beta4", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.moose <- jags(moose_con_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_habitat.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.moose)
  mcmcplot(tbd.moose$samples)
  save(tbd.moose, file = "./Outputs/TimeBtwnDetections/tbd.moose-season_habitat.RData")
  
  
  #'  --------------------------------
  #####  WHITE-TAILED DEER Analysis  ####
  #'  --------------------------------
  #'  Source JAGS model
  #'  Make sure inits and parameters being monitored match up with sourced model
  #'  Make sure model parameterization matches order of covariates in bundled data
  source("./Scripts/JAGS_models/JAGS_tbdconspecific_season_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(wtd_con_bundled$y, list(wtd_con_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "sigma", "mu.mu") #"beta3", "beta4", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wtd <- jags(wtd_con_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_habitat.txt',
                 inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                 n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wtd)
  mcmcplot(tbd.wtd$samples)
  save(tbd.wtd, file = "./Outputs/TimeBtwnDetections/tbd.wtd-season_habitat.RData")
  
  
  
 
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT  ####
  
