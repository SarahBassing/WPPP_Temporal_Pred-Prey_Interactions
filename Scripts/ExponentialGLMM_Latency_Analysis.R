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
  load("./Outputs/tbd_pred.prey_2022-09-23.RData") #2022-09-20
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
  
  ####  Setup data & MCMC specifications for JAGS  ####
  #'  ----------------------------------------------
  #'  MCMC settings
  nc <- 3; ni <- 5000; nb <- 1000; nt <- 5; na <- 200
  
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
                                    Species, spp_pair)) %>%
      mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
             Season = as.numeric(factor(Season, levels = c("Summer", "Fall", "Winter", "Spring"))), # levels must be 1-4 (not 0-3) for nested indexing
             PredatorID = as.numeric(factor(PredatorID, levels = c("Bobcat", "Coyote", "Black Bear", "Cougar", "Wolf"))), # levels must be 1-5 for nested indexing
             Complexity_index1 = scale(Complexity_index1),
             TRI = scale(TRI),
             PercForest = scale(PercForest),
             Monitoring = ifelse(Monitoring == "Trail", 1, 2),
             TRI_level = ifelse(backgroundRisk_TRI == "Low", 0, 1),
             PercFor_level = ifelse(backgroundRisk_For == "Low", 0, 1))
    print(summary(tbd_dat))
    print(head(tbd_dat))
    
    #'  Covariate matrix for JAGS
    covs <- matrix(NA, ncol = 8, nrow = ntbd)
    covs[,1] <- tbd_dat$Season
    covs[,2] <- tbd_dat$PredatorID
    covs[,3] <- tbd_dat$Complexity_index1
    covs[,4] <- tbd_dat$TRI
    covs[,5] <- tbd_dat$PercForest
    covs[,6] <- tbd_dat$Monitoring
    covs[,7] <- tbd_dat$TRI_level
    covs[,8] <- tbd_dat$PercFor_level
    head(covs)
    
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
  md_bundled <- bundle_dat(tbd_md_short)
  wtd_bundled <- bundle_dat(tbd_wtd_short)
  #'  Remove single wolf-elk observation
  tbd_elk_shorter <- tbd_elk_short %>% filter(PredatorID != "Wolf")
  elk_bundled <- bundle_dat(tbd_elk_shorter)
  moose_bundled <- bundle_dat(tbd_moose_short)
  #' #'  Remove single wolf-elk observation and bobcat observations since not a main 
  #' #'  predator of elk (keeping coyotes b/c prey on neonates)
  #' tbd_elk_shorter <- tbd_elk_short %>% filter(PredatorID != "Wolf") %>%
  #'   filter(PredatorID != "Bobcat")
  #' elk_bundled <- bundle_dat(tbd_elk_shorter)
  #' #'  Remove bobcat and coyote observations - not main predators of moose
  #' tbd_moose_shorter <- tbd_moose_short %>% filter(PredatorID != "Bobcat") %>%
  #'   filter(PredatorID != "Coyote")
  #' moose_bundled <- bundle_dat(tbd_moose_shorter)
  all_bundled <- bundle_dat(tbd_all_short)
  
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
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", "mu.tbd", "mu.mu") #"beta3", "beta4", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.md <- jags(md_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.md)
  mcmcplot(tbd.pred.md$samples)
  save(tbd.pred.md, file = "./Outputs/TimeBtwnDetections/tbd.pred.md-season_predID_habitat.RData")
  
  
  #'  -----------------------------
  #####  Predator - ELK Analysis  ####
  #'  -----------------------------
  #'  Source JAGS model
  #'  NOTE: NO bobcat or wolf observations in this model (consider dropping coyote too)
  #'  PredID 1 = black bear, PredID2 = cougar, PredID3 = coyote
  #'  Dropping interactions owing to relatively small sample size
  #'  Changes parameterization of JAGS model - source slightly different model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat_nowolf.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(elk_bundled$y, list(elk_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", "mu.tbd", "mu.mu")  #"beta3", "beta4", 
  
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
  #'  NOTE: NO bobcat or coyote observations in this model
  #'  PredID 1 = black bear, PredID2 = cougar, PredID3 = wolf
  #'  Dropping interactions owing to relatively small sample size
  #'  Changes parameterization of JAGS model - source slightly different model
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(moose_bundled$y, list(moose_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", "mu.tbd", "mu.mu")  #"beta3", "beta4", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.moose <- jags(moose_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.moose)
  mcmcplot(tbd.pred.moose$samples)  # sigma looks a little unhappy
  save(tbd.pred.moose, file = "./Outputs/TimeBtwnDetections/tbd.pred.moose-season_predID_habitat.RData")
  
  
  #'  -------------------------------------------
  #####  Predator - WHITE-TAILED DEER Analysis  ####
  #'  -------------------------------------------
  #'  Source JAGS model
  #'  Make sure inits and parameters being monitored match up with sourced model
  #'  Make sure model parameterization matches order of covariates in bundled data
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(wtd_bundled$y, list(wtd_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", "mu.tbd", "mu.mu") #"beta3", "beta4", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.wtd <- jags(wtd_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.wtd)
  mcmcplot(tbd.pred.wtd$samples)
  save(tbd.pred.wtd, file = "./Outputs/TimeBtwnDetections/tbd.pred.wtd-season_predID_habitat.RData")
  
  
  #'  -----------------------------------
  #####  Predator - ALL PREY Analysis  ####
  #'  -----------------------------------
  #'  Source JAGS model
  #'  Make sure inits and parameters being monitored match up with sourced model
  #'  Make sure model parameterization matches order of covariates in bundled data
  source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_habitat.R")
  
  #'  Set up initial values
  alpha.init <- log(aggregate(all_bundled$y, list(all_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))} 
  
  #'  Parameters to be monitored
  params <- c("alpha0", "beta", "beta1", "beta2", "sigma", "season.tbd", "pred.tbd", "mu.tbd", "mu.mu") #"beta3", "beta4", 
  
  #'  Run model
  start.time <- Sys.time()
  tbd.pred.all <- jags(all_bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_habitat.txt',
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.pred.all)
  mcmcplot(tbd.pred.all$samples)
  save(tbd.pred.all, file = "./Outputs/TimeBtwnDetections/tbd.pred.all-season_predID_habitat.RData")
  
  
  
  
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT - X^2 test  ####
  
  
  
  
 
  #' #####  Set up model in BUGS language  ####
  #' #'  -----------------------------------
  #' cat(file = './Outputs/TimeBtwnDetections/tbd_season_predID_HCI.txt', "
  #'     model{
  #'     
  #'     #'  Define priors
  #'     #'  -------------
  #'     #'  Prior for intercept
  #'     alpha0 ~ dnorm(0, 0.001)
  #'     
  #'     #' #'  Prior for beta coefficient (habitat complexity index)
  #'     #' beta ~ dnorm(0, 0.01)  # dunif(-10, 10)
  #'     
  #'     #'  Priors for categorical beta coefficients
  #'     #'  Season
  #'     beta1[1] <- 0
  #'     for(hh in 2:4){
  #'       beta1[hh] ~ dnorm(0, 0.01)  # dunif(-10,10)
  #'     }
  #'     #'  Predator species ID
  #'     beta2[1] <- 0
  #'     for(jj in 2:5){
  #'       beta2[jj] ~ dnorm(0, 0.01)  # dunif(-10,10)
  #'     }
  #'     
  #'     #'  Interaction for TRI * predator species ID
  #'     beta3[1] <- 0
  #'     for(jj in 2:5){
  #'       beta3[jj] ~ dnorm(0, 0.01)  # dunif(-10,10)
  #'     }
  #'     
  #'     #'  Interaction for PercForest * predator species ID
  #'     beta4[1] <- 0
  #'     for(jj in 2:5){
  #'       beta4[jj] ~ dnorm(0, 0.01)  # dunif(-10,10)
  #'     }
  #'     
  #'     #'  Prior for random effect for each camera location
  #'     for(j in 1:ncams){
  #'       alpha[j] ~ dnorm(0, tau.alpha) # mu.alpha for mean if no intercept (alpha0)
  #'     } 
  #'     
  #'     #'  Priors for TRI and PercForest
  #'     for(k in 1:2){  #ncovs
  #'       beta[k] ~ dnorm(0, 0.0001)
  #'       # beta[k] ~ dunif(-10, 10)
  #'     }
  #'     
  #'     #'  Hyperpriors for random effect
  #'     #'  -----------------------------
  #'     #'  Precision = 1/variance
  #'     # mu.alpha ~ dnorm(0, 0.001)
  #'     sigma ~ dunif(0, 10)   # sigma ~ dgamma(0.001, 0.001)  # gamma(shape, rate)
  #'     tau.alpha <- pow(sigma, -2) 
  #'     
  #'     #'  Define likelihood
  #'     #'  -----------------
  #'     for(i in 1:ntbd){
  #'       y[i] ~ dexp(lambda[i])   
  #'       #y[i] ~ dgamma(1, lambda[i]) # different parameterization, same result
  #' 
  #'       lambda[i] <- 1/mu[i]
  #'     
  #'       log(mu[i]) <- alpha0 + beta[1]*covs[i, 5] + beta[2]*covs[i, 6] + beta1[covs[i,1]] + beta2[covs[i,2]] + beta3[covs[i,2]]*covs[i, 5] + beta4[covs[i,2]]*covs[i, 6] + alpha[site[i]]
  #'       #log(mu[i]) <- alpha0 + beta[1]*covs[i, 5] + beta[2]*covs[i, 6] + beta1[covs[i,1]] + beta2[covs[i,2]] + alpha[site[i]]
  #'       #log(mu[i]) <- alpha0 + beta*covs[i, 4] + beta1[covs[i,1]] + beta2[covs[i,2]] + alpha[site[i]]
  #'       #log(mu[i]) <- alpha0 + beta*covs[i, 3] + beta1[covs[i,1]] + beta2[covs[i,2]] + beta3[covs[i,2]]*covs[i, 5] + alpha[site[i]]
  #'     }
  #'     
  #'     #'  Derived parameters
  #'     #'  ------------------
  #'     #'  mu.mu = mean number of minutes between events
  #'     mu.mu <- mean(mu[])
  #' 
  #'     }
  #'     ")
  #'     
  #'     
  #' #'  Define and bundle data
  #' tbd_pp <- tbd_md_short
  #' #'  Number of time-btwn-detection observations
  #' ntbd <- nrow(tbd_pp)
  #' #'  Number of unique camera locations
  #' ncams <- length(unique(tbd_pp$CameraLocation))
  #' #'  Format covariate data
  #' tbd_dat <- dplyr::select(tbd_pp, c(tbd_min, tbd_hour, tbd_day, CameraLocation,
  #'                                    Season, PredatorID, HuntingMode, TrophicLevel,
  #'                                    Complexity_index1, backgroundRisk_HCI,
  #'                                    backgroundRisk_TRI, backgroundRisk_For,
  #'                                    TRI, TRI_250m, PercForest, Species, spp_pair)) %>%
  #'   mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
  #'          Season = as.numeric(factor(Season, levels = c("Summer", "Fall", "Winter", "Spring"))), # levels must be 1-4 (not 0-3) for nested indexing
  #'          PredatorID = as.numeric(factor(PredatorID, levels = c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))), # levels must be 1-5 for nested indexing
  #'          HCI_level = ifelse(backgroundRisk_HCI == "Low", 0, 1),
  #'          Complexity_index1 = scale(Complexity_index1),
  #'          TRI = scale(TRI),
  #'          PercForest = scale(PercForest))
  #' 
  #' summary(tbd_dat)
  #' head(tbd_dat)
  #' #'  Covariate matrix for JAGS
  #' covs <- matrix(NA, ncol = 6, nrow = ntbd)
  #' covs[,1] <- tbd_dat$Season
  #' covs[,2] <- tbd_dat$PredatorID
  #' covs[,3] <- tbd_dat$HCI_level
  #' covs[,4] <- tbd_dat$Complexity_index1
  #' covs[,5] <- tbd_dat$TRI
  #' covs[,6] <- tbd_dat$PercForest
  #' head(covs)
  #' #'  Number of covariates
  #' ncovs <- ncol(covs)
  #' 
  #' #'  Time-between-detections
  #' tbd <- tbd_pp$tbd_min
  #' # tbd <- tbd_pp$tbd_hour
  #' # tbd <- tbd_pp$tbd_day
  #' summary(tbd)
  #' hist(tbd)
  #' 
  #' bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
  #'                 site = tbd_dat$cams)
  #' 
  #' #####  Initial values, monitor parameters and specify MCMC settings  ####
  #' #'  -----------------------------------------------------------------
  #' #'  Set up initial values
  #' alpha.init <- log(aggregate(tbd, list(tbd_dat$cams), FUN = mean)[,2])
  #' inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))}  #beta = runif(ncovs,-1,1)
  #' #why does it not need inits for alpha0, beta1, beta2, etc.???
  #' 
  #' #'  Parameters to be monitored
  #' params <- c("mu.mu", "alpha0", "beta", "beta1", "beta2", "beta3", "beta4", "sigma") # 
  #' 
  #' #'  MCMC settings
  #' nc <- 3; ni <- 50000; nb <- 10000; nt <- 5; na <- 2000
  #' 
  #' #'  Run model in JAGS
  #' #'  -----------------
  #' #'  Set up initial values
  #' alpha.init <- log(aggregate(tbd, list(tbd_dat$cams), FUN = mean)[,2])
  #' inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))}  #beta = runif(ncovs,-1,1)
  #' 
  #' #'  Parameters to be monitored
  #' params <- c("mu.mu", "alpha0", "beta", "beta1", "beta2", "beta3", "beta4", "sigma") # 
  #' 
  #' #'  Source JAGS model
  #' source("./Scripts/JAGS_models/JAGS_tbdpredprey_season_predID_X_habitat.R")
  #' 
  #' #'  Run model
  #' start.time <- Sys.time()
  #' tbd.mod <- jags(bundled, params, './Outputs/TimeBtwnDetections/tbd_season_predID_X_habitat.txt',
  #'                 inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
  #'                 n.adapt = na, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(tbd.mod)
  #' mcmcplot(tbd.mod$samples)
  #' # tbd.mod$mean
  #' # tbd.mod$summary
  #' # which(tbd.mod$summary[,"Rhat"] > 1.1)
  #' # save(tbd.mod, file="./Outputs/TimeBtwnDetections/tbd_global.Rdata")

  
  

  
  