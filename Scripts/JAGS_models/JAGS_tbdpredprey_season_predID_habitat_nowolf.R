  #'  ================================
  #'  Predator-Prey Latency Model 2
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Sept. 2022
  #'  ================================
  #'  Model sourced in ExponentialGLMM_Latency_Analysis.R
  #'  
  #'  Test the effects of season, predator ID, terrain ruggedness & percent forest 
  #'  on the latency between detections of predators followed by ungulate prey 
  #'  at camera sites. Same as JAGS_tbdpredprey_season_predID_habitat model but
  #'  wolves (beta2[5]) excluded from analysis owing to too few wolf-ungulate 
  #'  observations so only 4 predatorID categories in this parameterization.
  #'  tbd ~ season + predatorID + TRI + PercForest
  #'  ================================
  
  cat(file = './Outputs/TimeBtwnDetections/tbd_season_predID_habitat_nowolf.txt', "
          model{
          
          #'  Define priors
          #'  -------------
          #'  Prior for intercept
          alpha0 ~ dnorm(0, 0.001)
          
          #' #'  Priors for TRI, PercForest, & Study_Area
          #' for(k in 1:2){  #3
          #'   beta[k] ~ dnorm(0, 0.0001)
          #' }
      
          #'  Prior for TRI
          beta ~ dnorm(0, 0.0001)
          
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
          
            # log(mu[i]) <- alpha0 + beta1[covs[i,1]] + beta2[covs[i,2]] + 
            #               beta[1]*covs[i, 4] + beta[2]*covs[i, 5] + #beta[3]*covs[i, 9] + 
            #               alpha[site[i]]
            log(mu[i]) <- alpha0 + beta1[covs[i,1]] + beta2[covs[i,2]] + 
                          beta*covs[i, 4] + alpha[site[i]]
          }
          
      
          #'  Derived parameters
          #'  ------------------
          #'  Mean tbd per season & predator at mean TRI & PercForest
          for(hh in 1:4){
            for(jj in 1:4){
              # tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta[1]*0 + beta[2]*0)
              tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta*0)
            }
          }
      
          #' for(sa in 1:2){
          #'   for(hh in 1:4){
          #'     for(jj in 1:4){
          #'       tbd[sa, hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta[1]*0 + beta[2]*0 + beta[3]*(sa-1))
          #'     }
          #'   }
          #' }
          #' 
          #' #'  Mean tbd per season 
          #' for(sa in 1:2){
          #'   sa.tbd[sa] <- mean(tbd[sa,,])
          #' }
             
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
          # mu.mu <- mean(mu[])
      
          #' #'  Mean tbd per season & predator across range of TRI values
          #' for(i in 1:100){
          #'   for(hh in 1:4){
          #'     for(jj in 1:4){
          #'       tri.tbd[i, hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + 
          #'                                 beta[1]*newcovs[i,1] + beta[2]*0)
          #'     }
          #'   }
          #' }
          #' 
          #' #'  Mean tbd per predator across range of TRI values
          #' for(i in 1:100){
          #'   for(jj in 1:4){
          #'     pred.tbd.tri[i,jj] <- mean(tri.tbd[i,,jj])
          #'   }
          #' }
          #' 
          #' #'  Mean tbd per season & predator across range of % forest values
          #' for(i in 1:100){
          #'   for(hh in 1:4){
          #'     for(jj in 1:4){
          #'       for.tbd[i, hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + 
          #'                                 beta[1]*0 + beta[2]*newcovs[i,2])
          #'       } 
          #'     }
          #'   }
          #' 
          #' #'  Mean tbd per predator across range of % forest values
          #' for(i in 1:100){
          #'   for(jj in 1:4){
          #'     pred.tbd.for[i,jj] <- mean(for.tbd[i,,jj])
          #'   }
          #' }
      
          }
          ")