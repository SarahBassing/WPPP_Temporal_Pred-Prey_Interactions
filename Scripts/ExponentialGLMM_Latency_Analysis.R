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
  #'  Handy website for Poisson and exponential lambda interpretation:
  #'  https://towardsdatascience.com/what-is-exponential-distribution-7bdd08590e2a#:~:text=The%20definition%20of%20exponential%20distribution,Poisson%20(X%3D0)
  #'  The definition of exponential distribution is the probability distribution 
  #'  of the time *between* the events in a Poisson process.
  #'  Poisson distribution assumes that events occur independent of one another.
  #'  y ~ Exp(lambda) --> lambda is the reciprocal (1/λ) of the rate (λ) in a Poisson distribution.
  #'  lambda is often expressed in the amount of time before an event occurs.
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
  tbd_pred.prey <- filter(tbd_pred.prey, PredatorID != "Lynx") %>%
    mutate(tbd_round = round(TimeSinceLastDet, 0),
           tbd_min = TimeSinceLastDet,
           tbd_hour = tbd_min/60,
           tbd_day = tbd_hour/24)
  
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
      
      #'  Priors for beta coefficients
      #'  Season
      beta1[1] <- 0
      for(hh in 2:4){
        beta1[hh] ~ dnorm(0, 0.01)  # dunif(-10,10)
      }
      #'  Predator species ID
      beta2[1] <- 0
      for(jj in 2:5){
        beta2[jj] ~ dnorm(0, 0.01)  # dunif(-10,10)
      }
      #'  Prey species ID
      beta3[1] <- 0
      for(jj in 2:4){
        beta3[jj] ~ dnorm(0, 0.01)  # dunif(-10,10)
      }
      #'  Binary predator variables and habitat complexity
      for(k in 1:2){  #ncovs
        beta[k] ~ dnorm(0, 0.01)  # dunif(-10, 10)
      }
      
      #'  Hyperpriors for random intercept
      #'  --------------------------------
      #'  Mean and precision = 1/variance
      mu.alpha ~ dnorm(0, 0.001)
      sigma ~ dunif(0, 10)
      # sigma ~ dgamma(0.001, 0.001)  # gamma(shape, rate)
      tau.alpha <- pow(sigma, -2) 
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(lambda[i])
        lambda[i] <- log(tbd[i])
        #' HuntingMode + Complexity
        tbd[i] <- alpha[site[i]] + beta[1]*covs[i, 4] + beta[2]*covs[i, 7]
        #' #'  Global
        #' tbd[i] <- alpha[site[i]] + beta1[covs[i,1]] + beta2[covs[i,2]] + beta3[covs[i,3]] +
        #'           beta[1]*covs[i, 4] + beta[2]*covs[i, 5] + beta[3]*covs[i, 6] +
        #'           beta[4]*covs[i, 7]
      }
      
      #'  Derived parameters
      #'  ------------------
      mu.lambda <- mean(lambda[])
      
      }
      ")
  
  #'  Define and bundle data for JAGS
  #'  -------------------------------
  #'  Filter to a single ungulate species
  tbd_md <- tbd_pred.prey[tbd_pred.prey$Species == "Mule Deer",]
  #'  Number of time-btwn-detection observations
  ntbd <- nrow(tbd_md)
  #'  Number of unique camera locations
  ncams <- length(unique(tbd_md$CameraLocation))
  #'  Format covariate data
  tbd_dat <- dplyr::select(tbd_md, c(TimeSinceLastDet, CameraLocation, Season,
                                      PredatorID, HuntingMode, TrophicLevel,
                                      Complexity_index1, backgroundRisk, Species,
                                      spp_pair)) %>%
    mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
           Season = as.numeric(factor(Season), levels = c("Fall", "Spring", "Summer", "Winter")), # levels must be 1-4 (not 0-3) for nested indexing
           PredatorID = as.numeric(factor(PredatorID), levels = c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf")), # levels must be 1-5 for nested indexing
           Species = as.numeric(factor(Species), levels = c("Mule Deer")), # levels must be 1-4 for nested indexing
           HuntingMode = ifelse(HuntingMode == "Ambush", 0, 1),  # 0's OK here b/c not nested indexing and binary variable
           TrophicLevel = ifelse(TrophicLevel == "Apex", 0, 1),
           backgroundRisk = ifelse(backgroundRisk == "Low", 0, 1),
           Complexity_index1 = scale(Complexity_index1),
           spp_pair = as.numeric(factor(spp_pair), levels = spp_pair))
  summary(tbd_dat)
  head(tbd_dat)
  #'  Covariate matrix for JAGS
  covs <- matrix(NA, ncol = 7, nrow = ntbd)
  covs[,1] <- tbd_dat$Season
  covs[,2] <- tbd_dat$PredatorID
  covs[,3] <- tbd_dat$Species
  covs[,4] <- tbd_dat$HuntingMode
  covs[,5] <- tbd_dat$TrophicLevel
  covs[,6] <- tbd_dat$backgroundRisk
  covs[,7] <- tbd_dat$Complexity_index1
  head(covs)
  #'  Number of covariates
  ncovs <- ncol(covs)

  #'  Time-between-detections
  tbd <- tbd_md$TimeSinceLastDet
  tbd <- tbd_md$tbd_hour
  tbd <- tbd_md$tbd_day
  summary(tbd)

  bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                  site = tbd_dat$cams)

  #'  Initial values, monitor parameters and specify MCMC settings
  #'  ------------------------------------------------------------
  #'  Set up initial values
  alpha.init <- log(aggregate(tbd, list(tbd_dat$cams), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))}  #beta = runif(ncovs,-1,1)
  
  #'  Parameters to be monitored
  params <- c("alpha", "beta", "sigma", "lambda", "mu.lambda") # "beta1", "beta2", "beta3"
  
  #'  MCMC settings
  nc <- 3; ni <- 1000; nb <- 750; nt <- 5; na <- 100

  #'  Run model in JAGS
  #'  -----------------
  start.time <- Sys.time()
  tbd.mod <- jags(bundled, params, './Outputs/TimeBtwnDetections/global.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(tbd.mod$samples)
  print(tbd.mod)
  # tbd.mod$mean
  # tbd.mod$summary
  # which(tbd.mod$summary[,"Rhat"] > 1.1)
  # save(tbd.mod, file="./Outputs/TimeBtwnDetections/tbd_global.Rdata")
  
  
  
  #' #'  Define and bundle data for JAGS
  #' #'  -------------------------------
  #' #'  Number of time-btwn-detection observations
  #' ntbd <- nrow(tbd_pred.prey)
  #' #'  Number of unique camera locations
  #' ncams <- length(unique(tbd_pred.prey$CameraLocation))
  #' #'  Format covariate data
  #' tbd_dat <- dplyr::select(tbd_pred.prey, c(TimeSinceLastDet, CameraLocation, Season,
  #'                                           PredatorID, HuntingMode, TrophicLevel,
  #'                                           Complexity_index1, backgroundRisk, Species,
  #'                                           spp_pair)) %>%
  #'   mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
  #'          Season = as.numeric(factor(Season), levels = c("Fall", "Spring", "Summer", "Winter")), # levels must be 1-4 (not 0-3) for nested indexing
  #'          PredatorID = as.numeric(factor(PredatorID), levels = c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf")), # levels must be 1-5 for nested indexing
  #'          Species = as.numeric(factor(Species), levels = c("Elk", "Moose", "Mule Deer", "White-tailed Deer")), # levels must be 1-4 for nested indexing
  #'          HuntingMode = ifelse(HuntingMode == "Ambush", 0, 1),  # 0's OK here b/c not nested indexing and binary variable
  #'          TrophicLevel = ifelse(TrophicLevel == "Apex", 0, 1),
  #'          backgroundRisk = ifelse(backgroundRisk == "Low", 0, 1),
  #'          Complexity_index1 = scale(Complexity_index1),
  #'          spp_pair = as.numeric(factor(spp_pair), levels = spp_pair))
  #' summary(tbd_dat)
  #' head(tbd_dat)
  #' #'  Covariate matrix for JAGS
  #' covs <- matrix(NA, ncol = 7, nrow = ntbd)
  #' covs[,1] <- tbd_dat$Season
  #' covs[,2] <- tbd_dat$PredatorID
  #' covs[,3] <- tbd_dat$Species
  #' covs[,4] <- tbd_dat$HuntingMode
  #' covs[,5] <- tbd_dat$TrophicLevel
  #' covs[,6] <- tbd_dat$backgroundRisk
  #' covs[,7] <- tbd_dat$Complexity_index1
  #' head(covs)
  #' #'  Number of covariates
  #' ncovs <- ncol(covs)
  #' 
  #' #'  Time-between-detections
  #' tbd <- tbd_pred.prey$TimeSinceLastDet
  #' summary(tbd)
  #' 
  #' bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
  #'                 site = tbd_dat$cams)
  #' 
  #' 
  #' #'  Initial values, monitor parameters and specify MCMC settings
  #' #'  ------------------------------------------------------------
  #' #'  Set up initial values
  #' # inits <- function(){list(alpha = runif(ncams,-1,1), beta = runif(ncovs,-1,1))}
  #' alpha.init <- log(aggregate(tbd, list(tbd_dat$cams), FUN = mean)[,2])
  #' inits <- function(){list(alpha = alpha.init, beta = runif(4,-1,1))}  #beta = runif(ncovs,-1,1)
  #' 
  #' #'  Parameters to be monitored
  #' params <- c("alpha", "beta", "sigma", "mu.lambda") # "beta1", "beta2", "beta3"
  #' 
  #' #'  MCMC settings
  #' nc <- 3; ni <- 100000; nb <- 50000; nt <- 100; na <- 10000
  #' 
  #' #'  Run model in JAGS
  #' #'  -----------------
  #' start.time <- Sys.time()
  #' tbd.mod <- jags(bundled, params, './Outputs/TimeBtwnDetections/global.txt',
  #'                 inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
  #'                 n.adapt = na, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' mcmcplot(tbd.mod$samples)
  #' print(tbd.mod)
  #' tbd.mod$mean
  #' tbd.mod$summary
  #' which(tbd.mod$summary[,"Rhat"] > 1.1)
  #' save(tbd.mod, file="./Outputs/TimeBtwnDetections/tbd_global.Rdata")
  
  
  
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT  ####
  
  