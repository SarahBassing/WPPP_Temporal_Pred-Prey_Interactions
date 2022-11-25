  #'  ================================
  #'  Prey-Wolf Latency Model 2
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Nov. 2022
  #'  ================================
  #'  Model sourced in ExponentialGLMM_Latency_Analysis_Prey-Predator.R
  #'  
  #'  Test the effects of season, ungulate ID, and terrain ruggedness, on the latency
  #'  between detections of ungulates followed by predator at camera sites.
  #'  Excludes spring category from season covariate and elk & moose from ungualteID
  #'  covariate in the model owing to so observations within these categories.
  #'  tbd ~ season + ungulateID + TRI 
  #'  ================================
  
  cat(file = './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat_noSprgElkMoose.txt', "
              model{
              
              #'  Define priors
              #'  -------------
              #'  Prior for intercept
              alpha0 ~ dnorm(0, 0.001)
              
              #'  Prior for TRI
              beta ~ dnorm(0, 0.0001)
              
              #'  Priors for categorical covariates 
              #'  Season
              beta1[1] <- 0
              for(hh in 2:3){
                beta1[hh] ~ dnorm(0, 0.01)  
              }
            
              #'  Ungulate species ID
              beta2[1] <- 0
              for(jj in 2:2){
                beta2[jj] ~ dnorm(0, 0.01)  
              }
              
              #'  Prior for random effect for each camera location
              for(j in 1:ncams){
                alpha[j] ~ dnorm(0, tau.alpha) 
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
              
                log(mu[i]) <- alpha0 + beta1[covs[i,1]] + beta2[covs[i,2]] + 
                              beta*covs[i, 4] + alpha[site[i]]
              }
              
              #'  Derived parameters
              #'  ------------------
              #'  Mean tbd per season & ungulate at mean TRI
              for(hh in 1:3){
                for(jj in 1:2){
                  tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta*0)
                }
              }  
            
            
              #'  Mean tbd per season 
              for(hh in 1:3){
                season.tbd[hh] <- mean(tbd[hh,])
              }
            
              #'  Mean tbd per ungulate species
              for(jj in 1:2){
                prey.tbd[jj] <- mean(tbd[,jj])
              }
            
              #' Mean number of minutes between events
              mu.tbd <- mean(tbd[,])
              
            
              }
              ")
  
  
  
