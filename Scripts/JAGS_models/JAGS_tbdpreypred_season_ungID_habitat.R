  #'  ================================
  #'  Prey-Predator Latency Model 2
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Nov. 2022
  #'  ================================
  #'  Model sourced in ExponentialGLMM_Latency_Analysis_Prey-Predator.R
  #'  
  #'  Test the effects of season, ungulate ID, and terrain ruggedness,(no
  #'  interactions) on the latency between detections of ungulates followed by 
  #'  predator at camera sites.
  #'  tbd ~ season + ungulateID + TRI
  #'  ================================
  
  cat(file = './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat.txt', "
            model{
            
            #'  Define priors
            #'  -------------
            #'  Prior for intercept
            alpha0 ~ dnorm(0, 0.001)
            
            #'  Prior for TRI
            beta ~ dnorm(0, 0.0001)
            
            #'  Priors for categorical covariates and interactions
            #'  Season
            beta1[1] <- 0
            for(hh in 2:4){
              beta1[hh] ~ dnorm(0, 0.01)  
            }
          
            #'  Ungulate species ID
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
                            beta*covs[i, 4] + alpha[site[i]]
            }
            
            #'  Derived parameters
            #'  ------------------
            #'  Mean tbd per season & ungulate at mean TRI
            for(hh in 1:4){
              for(jj in 1:4){
                tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta*0)
              }
            }  
          
          
            #'  Mean tbd per season 
            for(hh in 1:4){
              season.tbd[hh] <- mean(tbd[hh,])
            }
          
            #'  Mean tbd per ungulate species
            for(jj in 1:4){
              prey.tbd[jj] <- mean(tbd[,jj])
            }
          
            #' Mean number of minutes between events
            mu.tbd <- mean(tbd[,])
            
          
            }
            ")
  
  
  
