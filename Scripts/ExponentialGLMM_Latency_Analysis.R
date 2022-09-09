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
      for(k in 1:4){  #ncovs
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
        tbd[i] <- alpha[site[i]] + beta[1]*covs[i, 4] + beta[4]*covs[i, 7]
        #' #'  Global
        #' tbd[i] <- alpha[site[i]] + beta1[covs[i,1]] + beta2[covs[i,2]] + beta3[covs[i,3]] +
        #'           beta[1]*covs[i, 4] + beta[2]*covs[i, 5] + beta[3]*covs[i, 6] +
        #'           beta[4]*covs[i, 7]
      }
      
      #'  Derived parameters
      #'  ------------------
      mu.lambda <- mean(lambda[])
      #' #'  Predict lambda under predation risk from ambush vs coursing predators at mean habitat complexity
      #' for(i in 1:ntbd){
      #'   for(hm in 1:2){
      #'     lam.pred[i, hm] <- exp(alpha[site[i]] + beta[1]*(hm - 1) + beta[4]*0)
      #'   }
      #' }
  
      # mu.lambda.pred <- mean(lam.pred[])
      
      }
      ")
  
  #'  Define and bundle data for JAGS
  #'  -------------------------------
  #'  Filter to a single ungulate species
  tbd_wtd <- tbd_pred.prey[tbd_pred.prey$Species == "White-tailed Deer",]
  #'  Number of time-btwn-deteciton observations
  ntbd <- nrow(tbd_wtd)
  #'  Number of unique camera locations
  ncams <- length(unique(tbd_wtd$CameraLocation))
  #'  Format covariate data
  tbd_dat <- dplyr::select(tbd_wtd, c(TimeSinceLastDet, CameraLocation, Season,
                                      PredatorID, HuntingMode, TrophicLevel,
                                      Complexity_index1, backgroundRisk, Species,
                                      spp_pair)) %>%
    mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
           Season = as.numeric(factor(Season), levels = c("Fall", "Spring", "Summer", "Winter")), # levels must be 1-4 (not 0-3) for nested indexing
           PredatorID = as.numeric(factor(PredatorID), levels = c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf")), # levels must be 1-5 for nested indexing
           Species = as.numeric(factor(Species), levels = c("White-tailed Deer")), # levels must be 1-4 for nested indexing
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
  tbd <- tbd_wtd$TimeSinceLastDet
  summary(tbd)

  bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                  site = tbd_dat$cams)

  #'  Initial values, monitor parameters and specify MCMC settings
  #'  ------------------------------------------------------------
  #'  Set up initial values
  # inits <- function(){list(alpha = runif(ncams,-1,1), beta = runif(ncovs,-1,1))}
  alpha.init <- log(aggregate(tbd, list(tbd_dat$cams), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(4,-1,1))}  #beta = runif(ncovs,-1,1)

  #'  Parameters to be monitored
  params <- c("alpha", "beta", "sigma", "mu.lambda") #"lam.pred", "lambda", "mu.lambda.pred"

  #'  MCMC settings
  nc <- 3; ni <- 75000; nb <- 50000; nt <- 25; na <- 10000

  #'  Run model in JAGS
  #'  -----------------
  start.time <- Sys.time()
  tbd.mod <- jags(bundled, params, './Outputs/TimeBtwnDetections/global.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(tbd.mod$samples)
  print(tbd.mod)
  tbd.mod$mean
  tbd.mod$summary
  which(tbd.mod$summary[,"Rhat"] > 1.1)
  save(tbd.mod, file="./Outputs/TimeBtwnDetections/tbd_global.Rdata")
  
  
  
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
           Season = as.numeric(factor(Season), levels = c("Fall", "Spring", "Summer", "Winter")), # levels must be 1-4 (not 0-3) for nested indexing
           PredatorID = as.numeric(factor(PredatorID), levels = c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf")), # levels must be 1-5 for nested indexing
           Species = as.numeric(factor(Species), levels = c("Elk", "Moose", "Mule Deer", "White-tailed Deer")), # levels must be 1-4 for nested indexing
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
  tbd <- tbd_pred.prey$TimeSinceLastDet
  summary(tbd)

  bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                  site = tbd_dat$cams)


  #'  Initial values, monitor parameters and specify MCMC settings
  #'  ------------------------------------------------------------
  #'  Set up initial values
  # inits <- function(){list(alpha = runif(ncams,-1,1), beta = runif(ncovs,-1,1))}
  alpha.init <- log(aggregate(tbd, list(tbd_dat$cams), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(4,-1,1))}  #beta = runif(ncovs,-1,1)

  #'  Parameters to be monitored
  params <- c("alpha", "beta", "sigma", "mu.lambda") #"lam.pred", "lambda", "mu.lambda.pred"

  #'  MCMC settings
  nc <- 3; ni <- 100000; nb <- 50000; nt <- 100; na <- 10000

  #'  Run model in JAGS
  #'  -----------------
  start.time <- Sys.time()
  tbd.mod <- jags(bundled, params, './Outputs/TimeBtwnDetections/global.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(tbd.mod$samples)
  print(tbd.mod)
  tbd.mod$mean
  tbd.mod$summary
  which(tbd.mod$summary[,"Rhat"] > 1.1)
  save(tbd.mod, file="./Outputs/TimeBtwnDetections/tbd_global.Rdata")
  
  
  
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT  ####
  
  