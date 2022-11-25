  #'  ================================
  #'  Ungulate-ungulate Latency Model 2
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Nov. 2022
  #'  ================================
  #'  Model sourced in ExponentialGLMM_Latency_Analysis_Ungulate.R
  #'  
  #'  Test the effects of season, ungulate species, terrain ruggedness, percent 
  #'  forest, and interactions between terrain, forest and ungulate species on 
  #'  the latency between detections of an ungulate followed by a different ungulate.
  #'  tbd ~ season + previousUngulate + TRI + PercForest + (previousUngulate*TRI) + (previousUngulate*PercForest)
  #'  ================================
  
  cat(file = './Outputs/TimeBtwnDetections/tbd_season_ungID_X_habitat.txt', "
        model{
        
        #'  Define priors
        #'  -------------
        #'  Prior for intercept
        alpha0 ~ dnorm(0, 0.001)
        
        #' #'  Priors for TRI, PercForest, and Study_Area
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
      
        #'  Ungulate species ID
        beta2[1] <- 0
        for(jj in 2:3){
          beta2[jj] ~ dnorm(0, 0.01)  
        }
        
        #'  Interaction for TRI * ungulate species ID
        beta3[1] <- 0
        for(jj in 2:3){
          beta3[jj] ~ dnorm(0, 0.01)  
        }
        
        #' #'  Interaction for PercForest * ungulate species ID
        #' beta4[1] <- 0
        #' for(jj in 2:3){
        #'   beta4[jj] ~ dnorm(0, 0.01)  
        #' }
        
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
        
          log(mu[i]) <- alpha0 + beta1[covs[i,1]] + beta2[covs[i,2]] + beta*covs[i, 3] +
                        beta3[covs[i,2]]*covs[i, 3] + alpha[site[i]]
          # log(mu[i]) <- alpha0 + beta1[covs[i,1]] + beta2[covs[i,2]] + 
          #               beta[1]*covs[i, 3] + beta[2]*covs[i, 4] + #beta[3]*covs[i, 9] +
          #               beta3[covs[i,2]]*covs[i, 3] + beta4[covs[i,2]]*covs[i, 4] + 
          #               alpha[site[i]]
        }
        
        #'  Derived parameters
        #'  ------------------
        #'  Mean tbd per season & ungulate at mean TRI #& PercForest
        for(hh in 1:4){
          for(jj in 1:3){
            tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta*0 + beta3[jj]*0)
            #tbd[hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta[1]*0 + beta[2]*0 + beta3[jj]*0 + beta4[jj]*0) 
          }
        }  
      
        #' for(sa in 1:2){
        #'   for(hh in 1:4){
        #'       for(jj in 1:3){
        #'         tbd[sa, hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta[1]*0 + beta[2]*0 +
        #'                        beta[3]*(sa-1) + beta3[jj]*0 + beta4[jj]*0)
        #'     }
        #'   }
        #' }
        #' 
        #' #'  Mean tbd per study area
        #' for(sa in 1:2){
        #'   sa.tbd[sa] <- mean(tbd[sa,,])
        #' }
      
        #'  Mean tbd per season 
        for(hh in 1:4){
          season.tbd[hh] <- mean(tbd[hh,])
        }
      
        #'  Mean tbd per ungulate 
        for(jj in 1:3){
          ung.tbd[jj] <- mean(tbd[,jj])
        }
      
        #' Mean number of minutes between events
        mu.tbd <- mean(tbd[,])
        # mu.mu <- mean(mu[])
    
      
        #' #'  Mean tbd per season & ungulate across range of TRI values
        #' for(i in 1:100){
        #'   for(hh in 1:4){
        #'     for(jj in 1:3){
        #'       tri.tbd[i, hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta[1]*newcovs[i,1] +
        #'                          beta[2]*0 + beta[3]*1 + beta3[jj]*newcovs[i,1] + beta4[jj]*0)
        #'     }
        #'   }
        #' }
        #' 
        #' #'  Mean tbd per ungulate across range of TRI values
        #' for(i in 1:100){
        #'   for(jj in 1:3){
        #'     pred.tbd.tri[i,jj] <- mean(tri.tbd[i,,jj])
        #'   }
        #' }
        #' 
        #' #'  Mean tbd per season & ungulate across range of % forest values
        #' for(i in 1:100){
        #'   for(hh in 1:4){
        #'     for(jj in 1:3){
        #'       for.tbd[i, hh, jj] <- exp(alpha0 + beta1[hh] + beta2[jj] + beta[1]*0 +
        #'                                 beta[2]*newcovs[i,2] + beta[3]*1 + beta3[jj]*0 +
        #'                                 beta4[jj]*beta[2]*newcovs[i,2])
        #'     }
        #'   }
        #' }
        #' 
        #' #'  Mean tbd per ungulate across range of % forest values
        #' for(i in 1:100){
        #'   for(jj in 1:3){
        #'     pred.tbd.for[i,jj] <- mean(for.tbd[i,,jj])
        #'   }
        #' }
        
      
        }
        ")



