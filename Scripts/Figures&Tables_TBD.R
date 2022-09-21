  #'  ================================
  #'  TBD figures & tables
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  Sept. 2022
  #'  ================================
  #'  Create result tables and figures based on results from latency analyses.
  #'  ================================
  
  #'  Clear workspace & load libraries
  rm(list = ls())
  
  library(ggplot2)
  library(tidyverse)
  
  #'  Load model results
  load("./Outputs/TimeBtwnDetections/tbd.pred.md-season_predID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.pred.elk-season_predID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.pred.moose-season_predID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.pred.wtd-season_predID_habitat.RData")
  
  load("./Outputs/TimeBtwnDetections/tbd.md-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.elk-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.moose-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.wtd-season_habitat.RData")

  #'  Pull out coefficient estimates
  coefs <- function(mod_out, spp) {
    Species <- spp
    Estimate <- round(unlist(mod_out$mean), 2)
    lci <- round(unlist(mod_out$q2.5), 2)
    uci <- round(unlist(mod_out$q97.5), 2)
    CI <- paste(" ",lci, "-", uci)
    out <- as.data.frame(cbind(Species, Estimate, CI))
    out <- tibble::rownames_to_column(out, "row_names") %>%
      relocate(row_names, .after = Species)
    colnames(out) <- c("Species", "Parameter", "Estimate", "95% CI")
    return(out)
  }
  pred.md.out <- coefs(tbd.pred.md, spp = "Mule deer") %>%
    mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
           Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
           Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
           Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
           Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
           Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
           Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter),
           Parameter = ifelse(Parameter == "beta21", "Predator: Black bear", Parameter),
           Parameter = ifelse(Parameter == "beta22", "Predator: Bobcat", Parameter),
           Parameter = ifelse(Parameter == "beta23", "Predator: Cougar", Parameter),
           Parameter = ifelse(Parameter == "beta24", "Predator: Coyote", Parameter),
           Parameter = ifelse(Parameter == "beta25", "Predator: Wolf", Parameter)) %>%
    filter(Estimate != 0)
  pred.md.out <- pred.md.out[1:10,]
  pred.wtd.out <- coefs(tbd.pred.wtd, spp = "White-tailed deer") %>%
    mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
           Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
           Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
           Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
           Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
           Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
           Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter),
           Parameter = ifelse(Parameter == "beta21", "Predator: Black bear", Parameter),
           Parameter = ifelse(Parameter == "beta22", "Predator: Bobcat", Parameter),
           Parameter = ifelse(Parameter == "beta23", "Predator: Cougar", Parameter),
           Parameter = ifelse(Parameter == "beta24", "Predator: Coyote", Parameter),
           Parameter = ifelse(Parameter == "beta25", "Predator: Wolf", Parameter)) %>%
    filter(Estimate != 0)
  pred.wtd.out <- pred.wtd.out[1:10,]
  #'  For Elk model: NO bobcat or wolf observations in this model 
  #'  PredID 1 = black bear, PredID2 = cougar, PredID3 = coyote
  pred.elk.out <- coefs(tbd.pred.elk, spp = "Elk") %>%
    mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
           Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
           Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
           Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
           Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
           Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
           Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter),
           Parameter = ifelse(Parameter == "beta21", "Predator: Black bear", Parameter),
           Parameter = ifelse(Parameter == "beta22", "Predator: Cougar", Parameter),
           Parameter = ifelse(Parameter == "beta23", "Predator: Coyote", Parameter)) %>%
    filter(Estimate != 0)
  pred.elk.out <- pred.elk.out[1:8,]
  #'  For moose model: NO bobcat or coyote observations in this model
  #'  PredID 1 = black bear, PredID2 = cougar, PredID3 = wolf
  pred.moose.out <- coefs(tbd.pred.moose, spp = "Moose") %>%
    mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
           Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
           Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
           Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
           Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
           Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
           Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter),
           Parameter = ifelse(Parameter == "beta21", "Predator: Black bear", Parameter),
           Parameter = ifelse(Parameter == "beta22", "Predator: Cougar", Parameter),
           Parameter = ifelse(Parameter == "beta23", "Predator: Wolf", Parameter)) %>%
    filter(Estimate != 0)
  pred.moose.out <- pred.moose.out[1:8,]
  
  pred.prey.out <- rbind(pred.elk.out, pred.moose.out, pred.md.out, pred.wtd.out)
  # write.csv(pred.prey.out, "./Outputs/TimeBtwnDetections/tbd.pred.prey_results_table.csv")
  
  md.out <- coefs(tbd.md, spp = "Mule deer") %>%
    mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
           Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
           Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
           Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
           Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
           Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
           Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter)) %>%
    filter(Estimate != 0)
  md.out <- md.out[1:10,]
  elk.out <- coefs(tbd.elk, spp = "Elk") %>%
    mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
           Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
           Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
           Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
           Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
           Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
           Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter)) %>%
    filter(Estimate != 0)
  elk.out <- elk.out[1:10,]
  moose.out <- coefs(tbd.moose, spp = "Moose") %>%
    mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
           Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
           Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
           Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
           Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
           Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
           Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter)) %>%
    filter(Estimate != 0)
  moose.out <- moose.out[1:10,]
  wtd.out <- coefs(tbd.wtd, spp = "White-tailed deer") %>%
    mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
           Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
           Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
           Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
           Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
           Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
           Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter)) %>%
    filter(Estimate != 0)
  wtd.out <- wtd.out[1:10,]
  
  conspecific.out <- rbind(elk.out, moose.out, md.out, wtd.out)
  # write.csv(conspecific.out, "./Outputs/TimeBtwnDetections/tbd.conspecific_results_table.csv")
  
  
  
  # md <- mutate(pred.md.out, 
  #              Parameter = ifelse(Parameter == "season.tbd1", "Mean TBD: Summer", Parameter),
  #              Parameter = ifelse(Parameter == "season.tbd2", "Mean TBD: Fall", Parameter),
  #              Parameter = ifelse(Parameter == "season.tbd3", "Mean TBD: Winter", Parameter),
  #              Parameter = ifelse(Parameter == "season.tbd4", "Mean TBD: Spring", Parameter),
  #              Parameter = ifelse(Parameter == "pred.tbd1", "Mean TBD: Black bear", Parameter),
  #              Parameter = ifelse(Parameter == "pred.tbd2", "Mean TBD: Bobcat", Parameter),
  #              Parameter = ifelse(Parameter == "pred.tbd3", "Mean TBD: Cougar", Parameter),
  #              Parameter = ifelse(Parameter == "pred.tbd4", "Mean TBD: Coyote", Parameter),
  #              Parameter = ifelse(Parameter == "pred.tbd5", "Mean TBD: Wolf", Parameter))
  # 
  # wtd <- mutate(pred.wtd.out,
  #               Parameter = ifelse(Parameter == "season.tbd1", "Mean TBD: Summer", Parameter),
  #               Parameter = ifelse(Parameter == "season.tbd2", "Mean TBD: Fall", Parameter),
  #               Parameter = ifelse(Parameter == "season.tbd3", "Mean TBD: Winter", Parameter),
  #               Parameter = ifelse(Parameter == "season.tbd4", "Mean TBD: Spring", Parameter),
  #               Parameter = ifelse(Parameter == "pred.tbd1", "Mean TBD: Black bear", Parameter),
  #               Parameter = ifelse(Parameter == "pred.tbd2", "Mean TBD: Bobcat", Parameter),
  #               Parameter = ifelse(Parameter == "pred.tbd3", "Mean TBD: Cougar", Parameter),
  #               Parameter = ifelse(Parameter == "pred.tbd4", "Mean TBD: Coyote", Parameter),
  #               Parameter = ifelse(Parameter == "pred.tbd5", "Mean TBD: Wolf", Parameter))
  # 
  # elk <- mutate(pred.elk.out,
  #               Parameter = ifelse(Parameter == "season.tbd1", "Mean TBD: Summer", Parameter),
  #               Parameter = ifelse(Parameter == "season.tbd2", "Mean TBD: Fall", Parameter),
  #               Parameter = ifelse(Parameter == "season.tbd3", "Mean TBD: Winter", Parameter),
  #               Parameter = ifelse(Parameter == "season.tbd4", "Mean TBD: Spring", Parameter),
  #               Parameter = ifelse(Parameter == "pred.tbd1", "Mean TBD: Black bear", Parameter),
  #               Parameter = ifelse(Parameter == "pred.tbd2", "Mean TBD: Cougar", Parameter),
  #               Parameter = ifelse(Parameter == "pred.tbd3", "Mean TBD: Coyote", Parameter))
  # 
  # moose <- mutate(pred.moose.out,
  #                 Parameter = ifelse(Parameter == "season.tbd1", "Mean TBD: Summer", Parameter),
  #                 Parameter = ifelse(Parameter == "season.tbd2", "Mean TBD: Fall", Parameter),
  #                 Parameter = ifelse(Parameter == "season.tbd3", "Mean TBD: Winter", Parameter),
  #                 Parameter = ifelse(Parameter == "season.tbd4", "Mean TBD: Spring", Parameter),
  #                 Parameter = ifelse(Parameter == "pred.tbd1", "Mean TBD: Black bear", Parameter),
  #                 Parameter = ifelse(Parameter == "pred.tbd2", "Mean TBD: Cougar", Parameter),
  #                 Parameter = ifelse(Parameter == "pred.tbd3", "Mean TBD: Wolf", Parameter))
  