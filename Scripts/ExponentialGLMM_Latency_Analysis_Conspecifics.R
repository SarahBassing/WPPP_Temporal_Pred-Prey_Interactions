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
  load("./Outputs/tbd_conspif_2022-09-20.RData") #2022-09-07 2022-09-13
  #'  Remove observations involving lynx due to too few lynx detections
  tbd_conspif <- tbd_conspif %>%
    #'  Change units of time
    mutate(tbd_min_round = round(TimeSinceLastDet, 0),
           tbd_min = TimeSinceLastDet,
           tbd_hour = tbd_min/60,
           tbd_day = tbd_hour/24) %>%
    dplyr::select(-TimeSinceLastDet)
  
  #'  Filter images by ungulate species
  tbd_md <- filter(tbd_conspif, Species == "Mule Deer")
  tbd_elk <- filter(tbd_conspif, Species == "Elk")
  tbd_moose <- filter(tbd_conspif, Species == "Moose")
  tbd_wtd <- filter(tbd_conspif, Species == "White-tailed Deer")
  
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
  tbd_md_short <- tbd_summary(tbd_md, spp = "mule deer", quant = 0.97)
  tbd_elk_short <- tbd_summary(tbd_elk, spp = "elk", quant = 0.97)
  tbd_moose_short <- tbd_summary(tbd_moose, spp = "moose", quant = 0.95)
  tbd_wtd_short <- tbd_summary(tbd_wtd, spp = "white-tailed deer", quant = 0.97)
  
  tbd_all_short <- rbind(tbd_md_short, tbd_elk_short, tbd_moose_short, tbd_wtd_short)
  
  
  #####  Set up model in BUGS language  ####
  #'  -----------------------------------
  cat(file = './Outputs/TimeBtwnDetections/tbd_conspecific.txt', "
        model{
        
        #'  Define priors
        #'  -------------
        #'  Prior for intercept
        alpha0 ~ dnorm(0, 0.001)
        
        #'  Priors for categorical beta coefficients
        #'  Season
        beta1[1] <- 0
        for(hh in 2:4){
          beta1[hh] ~ dnorm(0, 0.01)  # dunif(-10,10)
        }
        
        #'  Prior for random effect for each camera location
        for(j in 1:ncams){
          alpha[j] ~ dnorm(0, tau.alpha) # mu.alpha for mean if no intercept (alpha0)
        } 
        
        #'  Priors for TRI and PercForest
        for(k in 1:2){  #ncovs
          beta[k] ~ dnorm(0, 0.0001)
        }
        
        #'  Hyperpriors for random effect
        #'  -----------------------------
        #'  Precision = 1/variance
        sigma ~ dunif(0, 10)  
        tau.alpha <- pow(sigma, -2) 
        
        #'  Define likelihood
        #'  -----------------
        for(i in 1:ntbd){
          y[i] ~ dexp(lambda[i])   
  
          lambda[i] <- 1/mu[i]
        
          log(mu[i]) <- alpha0 + beta[1]*covs[i, 2] + beta[2]*covs[i, 3] + beta1[covs[i,1]] + alpha[site[i]]
        }
        
        #'  Derived parameters
        #'  ------------------
        #'  mu.mu = mean number of minutes between events
        mu.mu <- mean(mu[])
    
        }
        ")
  
  #####  Define and bundle data for JAGS  ####
  #'  -------------------------------------
  #'  Filter to a single ungulate species
  tbd_pp <- tbd_md_short
  #'  Number of time-btwn-detection observations
  ntbd <- nrow(tbd_pp)
  #'  Number of unique camera locations
  ncams <- length(unique(tbd_pp$CameraLocation))
  #'  Format covariate data
  tbd_dat <- dplyr::select(tbd_pp, c(tbd_min, tbd_hour, tbd_day, CameraLocation, 
                                     Season, Complexity_index1, TRI, TRI_250m,
                                     PercForest, Species)) %>%
    mutate(cams = as.numeric(factor(CameraLocation), levels = CameraLocation), # must be 1 - 313 (not 0 - 312) if using nested indexing for random effect
           Season = as.numeric(factor(Season, levels = c("Summer", "Fall", "Winter", "Spring"))), # levels must be 1-4 (not 0-3) for nested indexing
           Complexity_index1 = scale(Complexity_index1),
           TRI = scale(TRI),
           PercForest = scale(PercForest))
  
  summary(tbd_dat)
  head(tbd_dat)
  #'  Covariate matrix for JAGS
  covs <- matrix(NA, ncol = 4, nrow = ntbd)
  covs[,1] <- tbd_dat$Season
  covs[,2] <- tbd_dat$TRI
  covs[,3] <- tbd_dat$PercForest
  covs[,4] <- tbd_dat$Complexity_index1
  head(covs)
  #'  Number of covariates
  ncovs <- ncol(covs)
  
  #'  Time-between-detections
  tbd <- tbd_pp$tbd_min
  # tbd <- tbd_pp$tbd_hour
  # tbd <- tbd_pp$tbd_day
  summary(tbd)
  hist(tbd)
  
  bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                  site = tbd_dat$cams)
  
  #####  Initial values, monitor parameters and specify MCMC settings  ####
  #'  -----------------------------------------------------------------
  #'  Set up initial values
  alpha.init <- log(aggregate(tbd, list(tbd_dat$cams), FUN = mean)[,2])
  inits <- function(){list(alpha = alpha.init, beta = runif(2,-1,1))}  #beta = runif(ncovs,-1,1)
  #why does it not need inits for alpha0, beta1, beta2, etc.???
  
  #'  Parameters to be monitored
  params <- c("mu.mu", "alpha0", "beta", "beta1", "sigma") # 
  
  #'  MCMC settings
  nc <- 3; ni <- 50000; nb <- 10000; nt <- 5; na <- 2000
  
  #'  Run model in JAGS
  #'  -----------------
  start.time <- Sys.time()
  tbd.mod <- jags(bundled, params, './Outputs/TimeBtwnDetections/tbd_conspecific.txt',
                  inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.mod)
  mcmcplot(tbd.mod$samples)
  # tbd.mod$mean
  # tbd.mod$summary
  # which(tbd.mod$summary[,"Rhat"] > 1.1)
  # save(tbd.mod, file="./Outputs/TimeBtwnDetections/tbd_global.Rdata")
  
  
  
  ####  EVENTUALLY DO SOME ASSESSMENT OF GOODNESS OF FIT  ####
  
