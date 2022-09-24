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

  #'  Pull out and rename coefficient estimates 
  coefs <- function(mod_out, spp) {
    Species <- spp
    Estimate <- round(unlist(mod_out$mean), 2)
    lci <- round(unlist(mod_out$q2.5), 2)
    uci <- round(unlist(mod_out$q97.5), 2)
    CI <- paste(" ",lci, "-", uci)
    out <- as.data.frame(cbind(Species, Estimate, CI, lci, uci))
    out <- tibble::rownames_to_column(out, "row_names") %>%
      relocate(row_names, .after = Species)
    colnames(out) <- c("Species", "Parameter", "Estimate", "95% CI", "lci", "uci")
    renamed <- out %>%
      mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
             Parameter = ifelse(Parameter == "beta1", "Terrain ruggendess", Parameter),
             Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
             Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
             Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
             Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
             Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter),
             Parameter = ifelse(Parameter == "beta21", "Predator: Bobcat", Parameter),
             Parameter = ifelse(Parameter == "beta22", "Predator: Coyote", Parameter),
             Parameter = ifelse(Parameter == "beta23", "Predator: Black bear", Parameter),
             Parameter = ifelse(Parameter == "beta24", "Predator: Cougar", Parameter),
             Parameter = ifelse(Parameter == "beta25", "Predator: Wolf", Parameter),
             Parameter = ifelse(Parameter == "mu.tbd", "Mean TBD", Parameter),
             Parameter = ifelse(Parameter == "season.tbd1", "Mean TBD: Summer", Parameter),
             Parameter = ifelse(Parameter == "season.tbd2", "Mean TBD: Fall", Parameter),
             Parameter = ifelse(Parameter == "season.tbd3", "Mean TBD: Winter", Parameter),
             Parameter = ifelse(Parameter == "season.tbd4", "Mean TBD: Spring", Parameter),
             Parameter = ifelse(Parameter == "pred.tbd1", "Mean TBD: Bobcat", Parameter),
             Parameter = ifelse(Parameter == "pred.tbd2", "Mean TBD: Coyote", Parameter),
             Parameter = ifelse(Parameter == "pred.tbd3", "Mean TBD: Black bear", Parameter),
             Parameter = ifelse(Parameter == "pred.tbd4", "Mean TBD: Cougar", Parameter),
             Parameter = ifelse(Parameter == "pred.tbd5", "Mean TBD: Wolf", Parameter)) %>%
      filter(Estimate != 0)
    return(renamed)
  }
  pred.md.out <- coefs(tbd.pred.md, spp = "Mule deer") 
  pred.wtd.out <- coefs(tbd.pred.wtd, spp = "White-tailed deer") 
  pred.moose.out <- coefs(tbd.pred.moose, spp = "Moose")
  pred.elk.out <- coefs(tbd.pred.elk, spp = "Elk")
  
  con.md.out <- coefs(tbd.md, spp = "Mule deer") 
  con.wtd.out <- coefs(tbd.wtd, spp = "White-tailed deer") 
  con.moose.out <- coefs(tbd.moose, spp = "Moose")
  con.elk.out <- coefs(tbd.elk, spp = "Elk") 
  
  #'  Save coefficient estimates only
  pred.md.coef.out <- pred.md.out[1:10,1:4]
  pred.wtd.coef.out <- pred.wtd.out[1:10,1:4]
  pred.moose.coef.out <- pred.moose.out[1:10,1:4]
  pred.elk.coef.out <- pred.elk.out[1:9,1:4]
  
  con.md.coef.out <- con.md.out[1:6,1:4]
  con.wtd.coef.out <- con.wtd.out[1:6,1:4]
  con.moose.coef.out <- con.moose.out[1:6,1:4]
  con.elk.coef.out <- con.elk.out[1:6,1:4]
  
  pred.prey.coef.out <- rbind(pred.elk.coef.out, pred.moose.coef.out, pred.md.coef.out, pred.wtd.coef.out)
  write.csv(pred.prey.coef.out, "./Outputs/TimeBtwnDetections/tbd.pred.prey_coef_table.csv")
  
  conspif.coef.out <- rbind(con.elk.coef.out, con.moose.coef.out, con.md.coef.out, con.wtd.coef.out)
  write.csv(conspif.coef.out, "./Outputs/TimeBtwnDetections/tbd.conspecific_coef_table.csv")
  
  #'  Save mean time-between-detection results only
  pred.md.mutbd.out <- pred.md.out[12:21,1:4]
  pred.wtd.mutbd.out <- pred.wtd.out[12:21,1:4]
  pred.moose.mutbd.out <- pred.moose.out[12:21,1:4]
  pred.elk.mutbd.out <- pred.elk.out[11:19,1:4]
  
  con.md.mutbd.out <- con.md.out[8:12,1:4]
  con.wtd.mutbd.out <- con.wtd.out[8:12,1:4]
  con.moose.mutbd.out <- con.moose.out[8:12,1:4]
  con.elk.mutbd.out <- con.elk.out[8:12,1:4]
  
  pred.prey.mu.tbd.out <- rbind(pred.elk.mutbd.out, pred.moose.mutbd.out, pred.md.mutbd.out, pred.wtd.mutbd.out)
  write.csv(pred.prey.mu.tbd.out, "./Outputs/TimeBtwnDetections/tbd.pred.prey_meanTBD_table.csv")
  
  conspif.mu.tbd.out <- rbind(con.elk.mutbd.out, con.moose.mutbd.out, con.md.mutbd.out, con.wtd.mutbd.out)
  write.csv(conspif.mu.tbd.out, "./Outputs/TimeBtwnDetections/tbd.conspecific_meanTBD_table.csv")
  
  
  ####  Plot mean TBD  ####
  
  