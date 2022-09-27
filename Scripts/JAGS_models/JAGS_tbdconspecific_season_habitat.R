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
              
                log(mu[i]) <- alpha0 + beta1[covs[i,1]] + beta[1]*covs[i, 2] + beta[2]*covs[i, 3] + 
                              alpha[site[i]]
              }
              
              #'  Derived parameters
              #'  ------------------
              #'  Mean tbd per season at mean TRI & PercForest
              for(hh in 1:4){
                  season.tbd[hh] <- exp(alpha0 + beta1[hh] + beta[1]*0 + beta[2]*0)
              } 
      
              #' Mean number of minutes between events
              mu.tbd <- mean(season.tbd[])
      
              #'  Mean tbd per season across range of TRI values
              for(i in 1:100){
                for(hh in 1:4){
                  tri.tbd[i, hh] <- exp(alpha0 + beta1[hh] + beta[1]*newcovs[i,1] + beta[2]*0)
                }
              }

              #'  Mean tbd per season across range of % Forest values
              for(i in 1:100){
                for(hh in 1:4){
                  for.tbd[i, hh] <- exp(alpha0 + beta1[hh] + beta[1]*0 + beta[2]*newcovs[i,2])
                }
              }

              #' Mean number of minutes between events across range of TRI values
              for(i in 1:100){
                con.tbd.tri[i] <- mean(tri.tbd[i,])
              }

              #' Mean number of minutes between events across range of % Forest values
              for(i in 1:100){
                con.tbd.for[i] <- mean(for.tbd[i,])
              }
                
              }
              ")