  #'  ---------------------------------------------
  #'  Predation risk effects on latency of site use
  #'  WA Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Sept. 2022
  #'  ---------------------------------------------
  #'  Read in time-between-detection data (times between detection of a predator
  #'  followed by a prey species) and build generalized linear model with an 
  #'  exponential distribution in JAGS to estimate effect of predation risk on 
  #'  latency of prey use of camera sites following detection of a predator. 
  #'  Predation risk defined as "risky habitat" (e.g., habitat complexity) associated
  #'  with different predator hunting modes, presence of apex vs mesopredators,
  #'  and presence of ambush vs cursorial predators at a camera site.
  #'  
  #'  Time-between-detections data generated in TimeBetweenDetections.R script
  #'  ---------------------------------------------
  
  #'  Clear workspace & load libraries
  rm(list = ls())
  
  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  
  #'  Read in data
  load("./Outputs/tbd_pred.prey_2022-09-07.RData") 
  #'  Remove observations involving lynx due to too few lynx detections
  tbd_pred.prey <- filter(tbd_pred.prey, PredatorID != "Lynx")
  
  #'  Set up model in BUGS language
  #'  -----------------------------
  cat(file = './Outputs/TimeBtwnDetections/global.txt', "
      model{
      
      #'  Define priors
      #'  -------------
      #'  Prior for random intercept for each camera location
      for(j in 1:ncams){
        alpha[j] ~ dnorm(mu.alpha, tau.alpha)
      }
      
      #'  Prior for intercept
      mu ~ dunif(-10, 10)
      
      #'  Priors for beta coefficients
      for(k in 1:ncovs){
        beta[k] ~ dunif(-10, 10)
        # beta[k] ~ dnorm(0, 1e-04)
      }
      
      #'  Hyperpriors for random intercept
      #'  --------------------------------
      #'  Mean and precision = 1/variance
      mu.alpha ~ dnorm(0, 0.001)
      sigma ~ dunif(0, 10)
      # sigma ~ dgamma(0.001, 0.001) #gamma(shape, rate)
      tau.alpha <- pow(sigma, -2) 
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(lambda[i])
        lambda[i] <- log(tbd[i])
        tbd[i] <- mu + alpha[site[i]] + #beta[1]*covs[i, 1] + beta[2]*covs[i, 2] +
                  #beta[3]*covs[i, 3] + beta[4]*covs[i, 4] + beta[5]*covs[i, 5] +
                  beta[6]*covs[i, 6]
      }
      
      }
      ")
  
  #'  Define and bundle data for JAGS
  #'  -------------------------------
  #'  Number of time-btwn-deteciton observations
  ntbd <- nrow(tbd_pred.prey)
  #'  Number of unique camera locations
  ncams <- length(unique(tbd_pred.prey$CameraLocation))
  #'  Format covariate data
  tbd_dat <- dplyr::select(tbd_pred.prey, c(TimeSinceLastDet, CameraLocation, Season, 
                                            PredatorID, HuntingMode, TrophicLevel, 
                                            Complexity_index1, backgroundRisk, Species, 
                                            spp_pair)) %>%
    mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
           Season = as.numeric(factor(Season), levels = c("Summer", "Fall", "Winter", "Spring")) -1, # -1 changes levels from 1-4 to 0-3
           PredatorID = as.numeric(factor(PredatorID), levels = c("Cougar", "Black Bear", "Wolf", "Coyote", "Bobcat")) -1,
           HuntingMode = ifelse(HuntingMode == "Ambush", 0, 1), 
           TrophicLevel = ifelse(TrophicLevel == "Apex", 0, 1),
           backgroundRisk = ifelse(backgroundRisk == "Low", 0, 1),
           Complexity_index1 = scale(Complexity_index1),
           Species = as.numeric(factor(Species), levels = c("White-tailed Deer", "Mule Deer", "Moose", "Elk")) -1,
           spp_pair = as.numeric(factor(spp_pair), levels = spp_pair) -1)
  summary(tbd_dat)
  head(tbd_dat)
  #'  Covariate matrix for JAGS
  covs <- matrix(NA, ncol = 6, nrow = ntbd)
  covs[,1] <- tbd_dat$Season
  covs[,2] <- tbd_dat$PredatorID
  covs[,3] <- tbd_dat$HuntingMode
  covs[,4] <- tbd_dat$TrophicLevel
  covs[,5] <- tbd_dat$backgroundRisk
  covs[,6] <- tbd_dat$Complexity_index1
  head(covs)
  #'  Number of covariates
  ncovs <- ncol(covs)
  
  #'  Time-between-detections
  tbd <- tbd_pred.prey$TimeSinceLastDet
  summary(tbd)
  
  bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                  site = tbd_dat$cams)
  
  
  #'  Initial values, monitor parameters and specify MCMC settings
  #'  ------------------------------------------------------------
  #'  Set up initial values
  inits <- function(){list(alpha = runif(ncams,-1,1), beta = runif(ncovs,-1,1))}
  
  #'  Parameters to be monitored
  params <- c("lambda", "alpha", "beta", "sigma")
  
  #'  MCMC settings
  nc <- 3; ni <- 6000; nb <- 3000; nt <- 1; na <- 1000
  
  #'  Run model in JAGS
  #'  -----------------
  start.time <- Sys.time()
  tbd.mod <- jags(bundled, params, './Outputs/TimeBtwnDetections/global.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, 
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(tbd.mod$samples)
  which(tbd.mod$summary[,"Rhat"] > 1.1) 
  save(tbd.mod, file="./Outputs/TimeBtwnDetections/tbd_global.Rdata")
  
  
  
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT  ####
  
  