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
      # for(j in 1:ncams){
      #   alpha[j] ~ dnorm(0, tau)
      # }
      for(i in 1:ncams){
        alpha[cam_locs[i]] ~ dnorm(0, tau)
      }
      #'  Priors for beta coefficients
      for(i in 2:ncovs){
        beta[i] ~ dnorm(0, 1e-04)
      }
      
      #'  Hyperpriors for random intercept
      #'  --------------------------------
      # mu ~ dnorm(0, 0.01)
      sigma ~ dgamma(0.001, 0.001) #gamma(shape, rate)
      tau <- pow(sigma, -2) 
      
      #'  Define likelihood
      #'  -----------------
      # for(i in 1:ntbd){
      #   for(j in 1:ncams){
      #     y[i] ~ dexp(lambda[i, j])
      #     lambda[i, j] <- log(tbd[i, j])
      #     tbd[i, j] <- alpha[j]*covs[i, 1] + beta[1]*covs[i, 2] + beta[2]*covs[i, 3] + 
      #               beta[3]*covs[i, 4] + beta[4]*covs[i, 5] + beta[5]*covs[i, 6] + 
      #               beta[6]*covs[i, 7]
      #   }
      # }
      for(i in 1:ntbd){
          y[i] ~ dexp(lambda[i])
          lambda[i] <- log(tbd[i])
          tbd[i] <- alpha[cam_locs[i]] + beta[1]*covs[i, 2] + beta[2]*covs[i, 3] + 
                    beta[3]*covs[i, 4] + beta[4]*covs[i, 5] + beta[5]*covs[i, 6] + 
                    beta[6]*covs[i, 7]
      }
      
      }
      ")
  
  #'  Define and bundle data for JAGS
  #'  -------------------------------
  #'  Number of time-btwn-deteciton observations
  ntbd <- nrow(tbd_pred.prey)
  #'  Number of unique camera locations
  ncams <- length(unique(tbd_pred.prey$CameraLocation))
  #'  Format covariate data ---- NOT SURE IF THIS IS THE CORRECT WAY TO HANDLE CATEGORICAL VARIABLES
  tbd_dat <- dplyr::select(tbd_pred.prey, c(TimeSinceLastDet, CameraLocation, Season, PredatorID, HuntingMode, TrophicLevel, 
                                         Complexity_index1, backgroundRisk, Species, spp_pair)) %>%
    mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation),
           Season = as.numeric(factor(Season), levels = Season) -1,
           PredatorID = as.numeric(factor(PredatorID), levels = PredatorID) -1,
           HuntingMode = ifelse(HuntingMode == "Ambush", 0, 1),
           TrophicLevel = ifelse(TrophicLevel == "Apex", 0, 1),
           backgroundRisk = ifelse(backgroundRisk == "Low", 0, 1),
           Species = as.numeric(factor(Species), levels = Species) -1,
           spp_pair = as.numeric(factor(spp_pair), levels = spp_pair) -1)
  summary(tbd_dat)
  #'  Covariate matrix for JAGS
  covs <- matrix(NA, ncol = 7, nrow = ntbd)
  covs[,1] <- 1
  covs[,2] <- tbd_dat$Season
  covs[,3] <- tbd_dat$PredatorID
  covs[,4] <- tbd_dat$HuntingMode
  covs[,5] <- tbd_dat$TrophicLevel
  covs[,6] <- tbd_dat$backgroundRisk
  covs[,7] <- tbd_dat$Complexity_index1
  head(covs)
  
  #'  Time-between-detections
  tbd <- tbd_pred.prey$TimeSinceLastDet
  summary(tbd)
  
  bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncol(covs), ntbd = ntbd,
                  cam_locs = tbd_dat$cams)
  
  
  #'  Initial values, monitor parameters and specify MCMC settings
  #'  ------------------------------------------------------------
  #'  Set up initial values
  inits <- function(){list(alpha = runif(cam_locs,-1,1), beta = runif(ncol(covs),-1,1))}
  
  #'  Parameters to be monitored
  params <- c("lambda", "alpha", "beta", "sigma")
  
  #'  MCMC settings
  nc <- 3; ni <- 500; nb <- 300; nt <- 1; na <- 100
  
  #'  Run model in JAGS
  #'  -----------------
  start.time <- Sys.time()
  tbd.mod <- jags(bundled, params, './Outputs/TimeBtwnDetections/global.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, 
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(tbd.mod$samples)
  which(tbd.mod$summary[,"Rhat"] > 1.1) 
  # print(tbd.mod$summary[1:36, -c(4:6)], 3)
  save(tbd.mod, file="./Outputs/TimeBtwnDetections/tbd_global.Rdata")
  
  
  
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT  ####
  
  