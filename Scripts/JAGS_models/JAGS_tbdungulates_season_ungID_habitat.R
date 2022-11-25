  #'  ================================
  #'  Ungulate-ungulate Latency Model 
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Nov. 2022
  #'  ================================
  #'  Model sourced in ExponentialGLMM_Latency_Analysis_Ungulate.R
  #'  
  #'  Test the effects of season, ugulate species, terrain ruggedness, and percent 
  #'  forest on the latency between detections of an ungulate followed by a 
  #'  different ungulate species at camera sites.
  #'  tbd ~ season + previousUngulate + TRI + PercForest
  #'  ================================
  
  cat(file = './Outputs/TimeBtwnDetections/tbd_season_ungID_habitat.txt', "
            model{
            
            #'  Define priors
            #'  -------------
            #'  Prior for intercept
            alpha0 ~ dnorm(0, 0.001)
            
            #' #'  Priors for TRI and PercForest
            #' for(k in 1:2){  
            #'   beta[k] ~ dnorm(0, 0.0001)
            #' }
      
            #'  Prior for TRI
            beta ~ dnorm(0, 0.0001)
            
            #'  Priors for categorical covariates 
            #'  Season
            beta1[1] <- 0
            for(hh in 2:4){
              beta1[hh] ~ dnorm(0, 0.01)  
            }
          
            #'  Ungulate species ID 
            beta2[1] <- 0
            for(jj in 2:3){
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
                            beta*covs[i,3] + alpha[site[i]]
            }
            
            #'  Derived parameters
            #'  ------------------
            #'  Mean tbd per season & ungulate species at mean TRI #& PercForest
            for(hh in 1:4){
              for(jj in 1:3){
                tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta*0)
                # tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta[1]*0 + beta[2]*0)
              }
            } 
      
            #'  Mean tbd per season 
            for(hh in 1:4){
              season.tbd[hh] <- mean(tbd[hh,])
            }
      
            #'  Mean tbd per ungulate species 
            for(jj in 1:3){
              ung.tbd[jj] <- mean(tbd[,jj])
            }
      
            #' Mean number of minutes between events
            mu.tbd <- mean(tbd[,])
            mu.mu <- mean(mu[])
      
            #' #'  Mean tbd per season & ungulate species across range of TRI values
            #' for(i in 1:100){
            #'   for(hh in 1:4){
            #'     for(jj in 1:5){
            #'       tri.tbd[i, hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] +
            #'                                 beta[1]*newcovs[i,1] + beta[2]*0)
            #'     }
            #'   }
            #' }
            #' 
            #' #'  Mean tbd per ungulate species across range of TRI values
            #' for(i in 1:100){
            #'   for(jj in 1:5){
            #'     pred.tbd.tri[i,jj] <- mean(tri.tbd[i,,jj])
            #'   }
            #' }
            #' 
            #' #'  Mean tbd per season & ungulate species across range of % forest values
            #' for(i in 1:100){
            #'   for(hh in 1:4){
            #'     for(jj in 1:5){
            #'       for.tbd[i, hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] +
            #'                                 beta[1]*0 + beta[2]*newcovs[i,2])
            #'       }
            #'     }
            #'   }
            #' 
            #' #'  Mean tbd per ungulate species across range of % forest values
            #' for(i in 1:100){
            #'   for(jj in 1:5){
            #'     pred.tbd.for[i,jj] <- mean(for.tbd[i,,jj])
            #'   }
            #' }
      
        
            }
            ")
