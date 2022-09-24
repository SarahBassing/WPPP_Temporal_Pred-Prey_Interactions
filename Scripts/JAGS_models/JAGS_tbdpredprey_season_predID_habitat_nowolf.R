  #'  ================================
  #'  Predator-Prey Latency Model 2
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Sept. 2022
  #'  ================================
  #'  Model sourced in ExponentialGLMM_Latency_Analysis.R
  #'  
  #'  Test the effects of season, predator ID, terrain ruggedness, percent forest,
  #'  and interactions between terrain, forest and predator ID on the latency
  #'  between detections of predators followed by ungulate prey at camera sites.
  #'  tbd ~ season + predatorID + TRI + PercForest + (predatorID*TRI) + (predatorID*PercForest)
  #'  Note: two predators removed from list of predatorID, beta2, beta3, & beta4 
  #'  only have 3 categories each now. Dropping interactions as well owing to 
  #'  relatively small sample sizes.
  #'  ================================
  
  cat(file = './Outputs/TimeBtwnDetections/tbd_season_predID_habitat_nowolf.txt', "
          model{
          
          #'  Define priors
          #'  -------------
          #'  Prior for intercept
          alpha0 ~ dnorm(0, 0.001)
          
          #'  Priors for TRI and PercForest
          for(k in 1:2){  
            beta[k] ~ dnorm(0, 0.0001)
          }
          
          #'  Priors for categorical covariates and interactions
          #'  Season
          beta1[1] <- 0
          for(hh in 2:4){
            beta1[hh] ~ dnorm(0, 0.01)  
          }
        
          #'  Predator species ID (note: one less predator in this model)
          beta2[1] <- 0
          for(jj in 2:4){
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
                          beta[1]*covs[i, 4] + beta[2]*covs[i, 5] + alpha[site[i]]
          }
          
          #'  Derived parameters
          #'  ------------------
          #'  Mean tbd per season & predator at mean TRI & PercForest
            for(hh in 1:4){
              for(jj in 1:4){
                tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta[1]*0 + beta[2]*0)
              }
            } 
             
      
            #'  Mean tbd per season 
            for(hh in 1:4){
              season.tbd[hh] <- mean(tbd[hh,])
            }
      
            #'  Mean tbd per predator 
            for(jj in 1:4){
              pred.tbd[jj] <- mean(tbd[,jj])
            }
      
            #' Mean number of minutes between events
            mu.tbd <- mean(tbd[,])
            mu.mu <- mean(mu[])
        
      
          }
          ")