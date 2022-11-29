  #'  ================================
  #'  Conspecific Latency Model - intercept only
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Sept. 2022
  #'  ================================
  #'  Model sourced in ExponentialGLMM_Latency_Analysis_Conspecifics.R
  #'  
  #'  Test effects of season, terrain ruggedness, & percent forest on the latency
  #'  between detections of conspecific ungulates at camera sites.
  #'  ================================
  
  cat(file = './Outputs/TimeBtwnDetections/tbd_intercept_only.txt', "
                model{
                
                #'  Define priors
                #'  -------------
                #'  Prior for intercept
                alpha0 ~ dnorm(0, 0.001)
                
                #' #'  Prior for random effect for each camera location
                #' for(j in 1:ncams){
                #'   alpha[j] ~ dnorm(0, tau.alpha) 
                #' } 
                #' 
                #' #'  Hyperpriors for random effect
                #' #'  -----------------------------
                #' #'  Precision = 1/variance
                #' sigma ~ dunif(0, 10)   
                #' tau.alpha <- pow(sigma, -2) 
                
                #'  Define likelihood
                #'  -----------------
                for(i in 1:ntbd){
                  y[i] ~ dexp(lambda[i])   
          
                  lambda[i] <- 1/mu[i]
                
                  log(mu[i]) <- alpha0 #+ alpha[site[i]]
                }
      
                  
                #'  Derived parameters
                #'  ------------------
                mu.tbd <- exp(alpha0)
                  
                }
                ")