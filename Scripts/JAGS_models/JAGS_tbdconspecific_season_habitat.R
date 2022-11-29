  #'  ================================
  #'  Conspecific Latency Model
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Sept. 2022
  #'  ================================
  #'  Model sourced in ExponentialGLMM_Latency_Analysis_Conspecifics.R
  #'  
  #'  Test effects of season, terrain ruggedness, & percent forest on the latency
  #'  between detections of conspecific ungulates at camera sites.
  #'  tbd ~ season + TRI + PercForest
  #'  ================================
  
  cat(file = './Outputs/TimeBtwnDetections/tbd_season_habitat.txt', "
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
              
                # log(mu[i]) <- alpha0 + beta1[covs[i,1]] + beta[1]*covs[i, 2] + beta[2]*covs[i, 3] + 
                #               alpha[site[i]]  #beta[3]*covs[i, 5] + 
                log(mu[i]) <- alpha0 + beta1[covs[i,1]] + beta*covs[i, 2] + alpha[site[i]]
              }
              
              #'  Derived parameters
              #'  ------------------
              #'  Mean tbd per season at mean TRI & PercForest
              # for(sa in 1:2){
              #   for(hh in 1:4){
              #     season.sa.tbd[sa, hh] <- exp(alpha0 + beta1[hh] + beta[1]*0 + beta[2]*0 + beta[3]*(sa-1))
              #   } 
              # }
              # 
              # for(sa in 1:2){
              #   sa.tbd[sa] <- mean(season.sa.tbd[sa,])
              # }
              # 
              # for(hh in 1:4){
              #   season.tbd[hh] <- mean(season.sa.tbd[,hh])
              # }
                
              for(hh in 1:4){
                  # season.tbd[hh] <- exp(alpha0 + beta1[hh] + beta[1]*0 + beta[2]*0)
                  season.tbd[hh] <- exp(alpha0 + beta1[hh] + beta*0)
              }
      
              #' Mean number of minutes between events 
              #' Consider multiplying each season by total number of observations, 
              #' then average so it accounts for differences in sample size across seasons
              mu.tbd <- mean(season.tbd[])
      
              #' #'  Mean tbd per season across range of TRI values
              #' for(i in 1:100){
              #'   for(hh in 1:4){
              #'     tri.tbd[i, hh] <- exp(alpha0 + beta1[hh] + beta[1]*newcovs[i,1] + beta[2]*0 + beta[3]*1)
              #'   }
              #' }
              #' 
              #' #'  Mean tbd per season across range of % Forest values
              #' for(i in 1:100){
              #'   for(hh in 1:4){
              #'     for.tbd[i, hh] <- exp(alpha0 + beta1[hh] + beta[1]*0 + beta[2]*newcovs[i,2] + beta[3]*1)
              #'   }
              #' }
              #' 
              #' #' Mean number of minutes between events across range of TRI values
              #' for(i in 1:100){
              #'   con.tbd.tri[i] <- mean(tri.tbd[i,])
              #' }
              #' 
              #' #' Mean number of minutes between events across range of % Forest values
              #' for(i in 1:100){
              #'   con.tbd.for[i] <- mean(for.tbd[i,])
              #' }
                
              }
              ")