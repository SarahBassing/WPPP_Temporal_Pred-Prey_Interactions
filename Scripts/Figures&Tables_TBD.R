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
  library(stringr)
  library(tidyverse)
  library(khroma)
  
  #'  Load model results
  #'  Predator-Prey Analyses
  load("./Outputs/TimeBtwnDetections/tbd.pred.md-season_predID_X_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.pred.elk-season_predID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.pred.moose-season_predID_X_habitat.RData") 
  load("./Outputs/TimeBtwnDetections/tbd.pred.wtd-season_predID_X_habitat.RData")
  
  #'  Ungulate-Ungulate Analyses
  load("./Outputs/TimeBtwnDetections/tbd.ung.md-season_ungID_X_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.ung.elk-season_ungID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.ung.moose-season_ungID_X_habitat.RData") 
  load("./Outputs/TimeBtwnDetections/tbd.ung.wtd-season_ungID_X_habitat.RData")
  
  #'  Prey-Predator Analyses
  load("./Outputs/TimeBtwnDetections/tbd.prey.bear-season_ungID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.prey.bob-season_ungID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.prey.coug-season_ungID_habitat.RData") 
  load("./Outputs/TimeBtwnDetections/tbd.prey.coy-season_ungID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.prey.wolf-season_ungID_habitat.RData")
  
  #'  Conspecific Analyses based on 30-min interval for unique detection events
  load("./Outputs/TimeBtwnDetections/tbd.md-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.elk-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.moose-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.wtd-season_habitat.RData")
  
  #'  Conspecific Analyses based on 5-min interval for unique detection events
  load("./Outputs/TimeBtwnDetections/tbd.md-season_habitat_5min.RData")
  load("./Outputs/TimeBtwnDetections/tbd.elk-season_habitat_5min.RData")
  load("./Outputs/TimeBtwnDetections/tbd.moose-season_habitat_5min.RData")
  load("./Outputs/TimeBtwnDetections/tbd.wtd-season_habitat_5min.RData")

  #'  Pull out and rename coefficient estimates for pred.prey & conspecific models
  coefs <- function(mod_out, spp) {
    Species <- spp
    Estimate <- unlist(mod_out$mean)
    Median <- unlist(mod_out$q50)
    lci <- round(unlist(mod_out$q2.5), 3)
    uci <- round(unlist(mod_out$q97.5), 3)
    CI <- paste(" ",lci, "-", uci)
    out <- as.data.frame(cbind(Species, Estimate, Median, CI, lci, uci))
    out <- tibble::rownames_to_column(out, "row_names") %>%
      relocate(row_names, .after = Species)
    colnames(out) <- c("Species", "Parameter", "Estimate", "Median", "95% CI", "lci", "uci")
    renamed <- out %>%
      mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
             Parameter = ifelse(Parameter == "beta", "Terrain ruggendess", Parameter),
             # Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
             Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
             Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
             Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
             Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter),
             Parameter = ifelse(Parameter == "beta21", "Predator: Bobcat", Parameter),
             Parameter = ifelse(Parameter == "beta22", "Predator: Coyote", Parameter),
             Parameter = ifelse(Parameter == "beta23", "Predator: Black bear", Parameter),
             Parameter = ifelse(Parameter == "beta24", "Predator: Cougar", Parameter),
             Parameter = ifelse(Parameter == "beta25", "Predator: Wolf", Parameter),
             Parameter = ifelse(Parameter == "beta31", "Predator: Bobcat TRI", Parameter),
             Parameter = ifelse(Parameter == "beta32", "Predator: Coyote TRI", Parameter),
             Parameter = ifelse(Parameter == "beta33", "Predator: Black bear TRI", Parameter),
             Parameter = ifelse(Parameter == "beta34", "Predator: Cougar TRI", Parameter),
             Parameter = ifelse(Parameter == "beta35", "Predator: Wolf TRI", Parameter),
             # Parameter = ifelse(Parameter == "beta41", "Predator: Bobcat Forest", Parameter),
             # Parameter = ifelse(Parameter == "beta42", "Predator: Coyote Forest", Parameter),
             # Parameter = ifelse(Parameter == "beta43", "Predator: Black bear Forest", Parameter),
             # Parameter = ifelse(Parameter == "beta44", "Predator: Cougar Forest", Parameter),
             # Parameter = ifelse(Parameter == "beta45", "Predator: Wolf Forest", Parameter),
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
      filter(Estimate != 0) %>%
      mutate(Estimate = round(as.numeric(Estimate), 2),
             Median = round(as.numeric(Median), 2),
             lci = as.numeric(lci),
             uci = as.numeric(uci))
    return(renamed)
  }
  pred.md.out <- coefs(tbd.pred.md, spp = "Mule deer")
  pred.elk.out <- coefs(tbd.pred.elk, spp = "Elk")
  pred.moose.out <- coefs(tbd.pred.moose, spp = "Moose")
  pred.wtd.out <- coefs(tbd.pred.wtd, spp = "White-tailed deer") 

  con.md.out <- coefs(tbd.md, spp = "Mule deer") 
  con.wtd.out <- coefs(tbd.wtd, spp = "White-tailed deer") 
  con.moose.out <- coefs(tbd.moose, spp = "Moose")
  con.elk.out <- coefs(tbd.elk, spp = "Elk") 
  
  #'  Pull out and rename coefficient estimates for prey.pred models
  coefs <- function(mod_out, spp) {
    Species <- spp
    Estimate <- unlist(mod_out$mean)
    Median <- unlist(mod_out$q50)
    lci <- round(unlist(mod_out$q2.5), 3)
    uci <- round(unlist(mod_out$q97.5), 3)
    CI <- paste(" ",lci, "-", uci)
    out <- as.data.frame(cbind(Species, Estimate, Median, CI, lci, uci))
    out <- tibble::rownames_to_column(out, "row_names") %>%
      relocate(row_names, .after = Species)
    colnames(out) <- c("Species", "Parameter", "Estimate", "Median", "95% CI", "lci", "uci")
    renamed <- out %>%
      mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
             Parameter = ifelse(Parameter == "beta", "Terrain ruggendess", Parameter),
             # Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
             Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
             Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
             Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
             Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter),
             Parameter = ifelse(Parameter == "beta21", "Ungulate: Elk", Parameter),
             Parameter = ifelse(Parameter == "beta22", "Ungulate: Moose", Parameter),
             Parameter = ifelse(Parameter == "beta23", "Ungulate: Mule deer", Parameter),
             Parameter = ifelse(Parameter == "beta24", "Ungulate: White-tailed deer", Parameter),
             # Parameter = ifelse(Parameter == "beta31", "Ungulate: Elk TRI", Parameter),
             # Parameter = ifelse(Parameter == "beta32", "Ungulate: Moose TRI", Parameter),
             # Parameter = ifelse(Parameter == "beta33", "Ungulate: Mule deer TRI", Parameter),
             # Parameter = ifelse(Parameter == "beta34", "Ungulate: White-tailed deer TRI", Parameter),
             Parameter = ifelse(Parameter == "mu.tbd", "Mean TBD", Parameter),
             Parameter = ifelse(Parameter == "season.tbd1", "Mean TBD: Summer", Parameter),
             Parameter = ifelse(Parameter == "season.tbd2", "Mean TBD: Fall", Parameter),
             Parameter = ifelse(Parameter == "season.tbd3", "Mean TBD: Winter", Parameter),
             Parameter = ifelse(Parameter == "season.tbd4", "Mean TBD: Spring", Parameter),
             Parameter = ifelse(Parameter == "prey.tbd1", "Mean TBD: Elk", Parameter),
             Parameter = ifelse(Parameter == "prey.tbd2", "Mean TBD: Moose", Parameter),
             Parameter = ifelse(Parameter == "prey.tbd3", "Mean TBD: Mule deer", Parameter),
             Parameter = ifelse(Parameter == "prey.tbd4", "Mean TBD: White-tailed deer", Parameter)) %>%
      filter(Estimate != 0) %>%
      mutate(Estimate = round(as.numeric(Estimate), 2),
             Median = round(as.numeric(Median), 2),
             lci = as.numeric(lci),
             uci = as.numeric(uci))
    return(renamed)
  }
  prey.bear.out <- coefs(tbd.prey.bear, spp = "Black bear") %>%
    mutate(Parameter = ifelse(Parameter == "Season: Winter", "Season: Spring", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: Winter", "Mean TBD: Spring", Parameter))
  prey.bob.out <- coefs(tbd.prey.bob, spp = "Bobcat")
  prey.coug.out <- coefs(tbd.prey.coug, spp = "Cougar")
  prey.coy.out <- coefs(tbd.prey.coy, spp = "Coyote") 
  prey.wolf.out <- coefs(tbd.prey.wolf, spp = "Wolf") %>%
    mutate(Parameter = ifelse(Parameter == "Ungulate: Moose", "Ungulate: White-tailed deer", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: Elk", "Mean TBD: Mule deer", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: Moose", "Mean TBD: White-tailed deer", Parameter))
  
  
  #'  Pull out and rename coefficient estimates for ungulate-ungulate models
  coefs <- function(mod_out, spp) {
    Species <- spp
    Estimate <- unlist(mod_out$mean)
    Median <- unlist(mod_out$q50)
    lci <- round(unlist(mod_out$q2.5), 3)
    uci <- round(unlist(mod_out$q97.5), 3)
    CI <- paste(" ",lci, "-", uci)
    out <- as.data.frame(cbind(Species, Estimate, Median, CI, lci, uci))
    out <- tibble::rownames_to_column(out, "row_names") %>%
      relocate(row_names, .after = Species)
    colnames(out) <- c("Species", "Parameter", "Estimate", "Median", "95% CI", "lci", "uci")
    renamed <- out %>%
      mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
             Parameter = ifelse(Parameter == "beta", "Terrain ruggendess", Parameter),
             # Parameter = ifelse(Parameter == "beta2", "Percent forest", Parameter),
             Parameter = ifelse(Parameter == "beta11", "Season: Summer", Parameter),
             Parameter = ifelse(Parameter == "beta12", "Season: Fall", Parameter),
             Parameter = ifelse(Parameter == "beta13", "Season: Winter", Parameter),
             Parameter = ifelse(Parameter == "beta14", "Season: Spring", Parameter),
             Parameter = ifelse(Parameter == "beta21", "Ungulate: 1", Parameter),
             Parameter = ifelse(Parameter == "beta22", "Ungulate: 2", Parameter),
             Parameter = ifelse(Parameter == "beta23", "Ungulate: 3", Parameter),
             Parameter = ifelse(Parameter == "beta31", "Ungulate: 1 TRI", Parameter),
             Parameter = ifelse(Parameter == "beta32", "Ungulate: 2 TRI", Parameter),
             Parameter = ifelse(Parameter == "beta33", "Ungulate: 3 TRI", Parameter),
             # Parameter = ifelse(Parameter == "beta41", "Ungulate: 1 Forest", Parameter),
             # Parameter = ifelse(Parameter == "beta42", "Ungulate: 2 Forest", Parameter),
             # Parameter = ifelse(Parameter == "beta43", "Ungulate: 3 Forest", Parameter),
             Parameter = ifelse(Parameter == "mu.tbd", "Mean TBD", Parameter),
             Parameter = ifelse(Parameter == "season.tbd1", "Mean TBD: Summer", Parameter),
             Parameter = ifelse(Parameter == "season.tbd2", "Mean TBD: Fall", Parameter),
             Parameter = ifelse(Parameter == "season.tbd3", "Mean TBD: Winter", Parameter),
             Parameter = ifelse(Parameter == "season.tbd4", "Mean TBD: Spring", Parameter),
             Parameter = ifelse(Parameter == "ung.tbd1", "Mean TBD: 1", Parameter),
             Parameter = ifelse(Parameter == "ung.tbd2", "Mean TBD: 2", Parameter),
             Parameter = ifelse(Parameter == "ung.tbd3", "Mean TBD: 3", Parameter)) %>%
      filter(Estimate != 0) %>%
      mutate(Estimate = round(as.numeric(Estimate), 2),
             Median = round(as.numeric(Median), 2),
             lci = as.numeric(lci),
             uci = as.numeric(uci))
    return(renamed)
  }
  ung.md.out <- coefs(tbd.ung.md, spp = "Mule deer") %>%
    mutate(Parameter = ifelse(Parameter == "Ungulate: 1", "Ungulate: White-tailed Deer", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 2", "Ungulate: Elk", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 3", "Ungulate: Moose", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 1 TRI", "Ungulate: White-tailed Deer TRI", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 2 TRI", "Ungulate: Elk TRI", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 3 TRI", "Ungulate: Moose TRI", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 1", "Mean TBD: White-tailed Deer", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 2", "Mean TBD: Elk", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 3", "Mean TBD: Moose", Parameter))
  ung.elk.out <- coefs(tbd.ung.elk, spp = "Elk") %>%
    mutate(Parameter = ifelse(Parameter == "Ungulate: 1", "Ungulate: Mule Deer", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 2", "Ungulate: White-tailed Deer", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 3", "Ungulate: Moose", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 1 TRI", "Ungulate: Mule Deer TRI", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 2 TRI", "Ungulate: White-tailed Deer TRI", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 3 TRI", "Ungulate: Moose TRI", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 1", "Mean TBD: Mule Deer", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 2", "Mean TBD: White-tailed Deer", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 3", "Mean TBD: Moose", Parameter))
  ung.moose.out <- coefs(tbd.ung.moose, spp = "Moose") %>%
    mutate(Parameter = ifelse(Parameter == "Ungulate: 1", "Ungulate: Mule Deer", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 2", "Ungulate: White-tailed Deer", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 3", "Ungulate: Elk", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 1 TRI", "Ungulate: Mule Deer TRI", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 2 TRI", "Ungulate: White-tailed Deer TRI", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 3 TRI", "Ungulate: Elk TRI", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 1", "Mean TBD: Mule Deer", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 2", "Mean TBD: White-tailed Deer", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 3", "Mean TBD: Elk", Parameter))
  ung.wtd.out <- coefs(tbd.ung.wtd, spp = "White-tailed deer") %>%
    mutate(Parameter = ifelse(Parameter == "Ungulate: 1", "Ungulate: Mule Deer", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 2", "Ungulate: Elk", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 3", "Ungulate: Moose", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 1 TRI", "Ungulate: Mule Deer TRI", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 2 TRI", "Ungulate: Elk TRI", Parameter),
           Parameter = ifelse(Parameter == "Ungulate: 3 TRI", "Ungulate: Moose TRI", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 1", "Mean TBD: Mule Deer", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 2", "Mean TBD: Elk", Parameter),
           Parameter = ifelse(Parameter == "Mean TBD: 3", "Mean TBD: Moose", Parameter))
  
  
  #'  Save coefficient estimates only
  pred.md.coef.out <- pred.md.out[1:13,1:5] 
  pred.elk.coef.out <- pred.elk.out[1:8,1:5] 
  pred.moose.coef.out <- pred.moose.out[1:13,1:5] 
  pred.wtd.coef.out <- pred.wtd.out[1:13,1:5] 
  
  ung.md.coef.out <- ung.md.out[1:9,1:5] 
  ung.elk.coef.out <- ung.elk.out[1:7,1:5]
  ung.moose.coef.out <- ung.moose.out[1:9,1:5] 
  ung.wtd.coef.out <- ung.wtd.out[1:9,1:5] 
  
  prey.bear.coef.out <- prey.bear.out[1:7,1:5] 
  prey.bob.coef.out <- prey.bob.out[1:8,1:5]
  prey.coug.coef.out <- prey.coug.out[1:8,1:5] 
  prey.coy.coef.out <- prey.coy.out[1:8,1:5]
  prey.wolf.coef.out <- prey.wolf.out[1:5,1:5]
  
  con.md.coef.out <- con.md.out[1:5,1:5] 
  con.elk.coef.out <- con.elk.out[1:5,1:5] 
  con.moose.coef.out <- con.moose.out[1:5,1:5] 
  con.wtd.coef.out <- con.wtd.out[1:5,1:5] 
  
  pred.prey.coef.out <- rbind(pred.elk.coef.out, pred.moose.coef.out, pred.md.coef.out, pred.wtd.coef.out)
  # write.csv(pred.prey.coef.out, "./Outputs/TimeBtwnDetections/Tables/tbd.pred.prey_coef_table.csv")
  
  ungulate.coef.out <- rbind(ung.elk.coef.out, ung.moose.coef.out, ung.md.coef.out, ung.wtd.coef.out)
  # write.csv(ungulate.coef.out, "./Outputs/TimeBtwnDetections/Tables/tbd.ungulate_coef_table.csv")
  
  prey.pred.coef.out <- rbind(prey.bear.coef.out, prey.bob.coef.out, prey.coug.coef.out, prey.coy.coef.out, prey.wolf.coef.out)
  # write.csv(prey.pred.coef.out, "./Outputs/TimeBtwnDetections/Tables/tbd.prey.pred_coef_table.csv")
  
  conspif.coef.out <- rbind(con.elk.coef.out, con.moose.coef.out, con.md.coef.out, con.wtd.coef.out)
  # write.csv(conspif.coef.out, "./Outputs/TimeBtwnDetections/Tables/tbd.conspecific_coef_table.csv")
  
  #'  Save mean time-between-detection results only
  #'  But first change 95% CRI to 2 digits
  short_CI <- function(out) {
    out <- out %>%
      mutate(lci = round(lci, 2),
             uci = round(uci, 2), 
             CRI = paste(" ",lci, "-", uci)) %>%
      dplyr::select(-c(`95% CI`)) %>%
      relocate(CRI, .before = lci) %>%
      rename(`95% CI` = CRI)
    return(out)
  }
  pred.md.out2 <- short_CI(pred.md.out)
  pred.elk.out2 <- short_CI(pred.elk.out)
  pred.moose.out2 <- short_CI(pred.moose.out)
  pred.wtd.out2 <- short_CI(pred.wtd.out)
  
  pred.md.mutbd.out <- pred.md.out2[15:24,1:7] #pred.md.out[20:29,1:6]
  pred.elk.mutbd.out <- pred.elk.out2[10:18,1:7] #pred.elk.out[11:19,1:6]
  pred.moose.mutbd.out <- pred.moose.out2[15:24,1:7] #pred.moose.out[20:29,1:6]
  pred.wtd.mutbd.out <- pred.wtd.out2[15:24,1:7] #pred.wtd.out[20:29,1:6]
  
  ung.md.out2 <- short_CI(ung.md.out)
  ung.elk.out2 <- short_CI(ung.elk.out)
  ung.moose.out2 <- short_CI(ung.moose.out)
  ung.wtd.out2 <- short_CI(ung.wtd.out)
  
  ung.md.mutbd.out <- ung.md.out2[11:18,1:7] #ung.md.out[11:18,1:6]
  ung.elk.mutbd.out <- ung.elk.out2[9:16,1:7] #ung.elk.out[9:13,1:6]
  ung.moose.mutbd.out <- ung.moose.out2[11:18,1:7] #ung.moose.out[11:18,1:6]
  ung.wtd.mutbd.out <- ung.wtd.out2[11:18,1:7] #ung.wtd.out[11:18,1:6]
  
  prey.bear.out2 <- short_CI(prey.bear.out)
  prey.bob.out2 <- short_CI(prey.bob.out)
  prey.coug.out2 <- short_CI(prey.coug.out)
  prey.coy.out2 <- short_CI(prey.coy.out)
  prey.wolf.out2 <- short_CI(prey.wolf.out)
  
  prey.bear.mutbd.out <- prey.bear.out2[9:16,1:7] #prey.md.out[20:29,1:6]
  prey.bob.mutbd.out <- prey.bob.out2[10:18,1:7] #prey.elk.out[11:19,1:6]
  prey.coug.mutbd.out <- prey.coug.out2[10:18,1:7] #prey.moose.out[20:29,1:6]
  prey.coy.mutbd.out <- prey.coy.out2[10:18,1:7] #prey.wtd.out[20:29,1:6]
  prey.wolf.mutbd.out <- prey.wolf.out2[7:12,1:7] #prey.wtd.out[20:29,1:6]
  
  con.md.out2 <- short_CI(con.md.out)
  con.wtd.out2 <- short_CI(con.wtd.out)
  con.moose.out2 <- short_CI(con.moose.out)
  con.elk.out2 <- short_CI(con.elk.out)
  
  con.md.mutbd.out <- con.md.out2[7:11,1:7] #con.md.out[8:12,1:6]
  con.wtd.mutbd.out <- con.wtd.out2[7:11,1:7] #con.wtd.out[8:12,1:6]
  con.moose.mutbd.out <- con.moose.out2[7:11,1:7] #con.moose.out[8:12,1:6]
  con.elk.mutbd.out <- con.elk.out2[7:11,1:7] #con.elk.out[8:12,1:6]
  
  pred.prey.mu.tbd.out <- rbind(pred.elk.mutbd.out, pred.moose.mutbd.out, pred.md.mutbd.out, pred.wtd.mutbd.out) 
  pred.prey.mu.tbd.tbl <- dplyr::select(pred.prey.mu.tbd.out, -c("lci", "uci"))
  # write.csv(pred.prey.mu.tbd.tbl, "./Outputs/TimeBtwnDetections/Tables/tbd.pred.prey_meanTBD_table.csv")
  
  ungulate.mu.tbd.out <- rbind(ung.elk.mutbd.out, ung.moose.mutbd.out, ung.md.mutbd.out, ung.wtd.mutbd.out) 
  ungulate.mu.tbd.tbl <- dplyr::select(ungulate.mu.tbd.out, -c("lci", "uci"))
  # write.csv(ungulate.mu.tbd.tbl, "./Outputs/TimeBtwnDetections/Tables/tbd.ungulate_meanTBD_table.csv")
  
  prey.pred.mu.tbd.out <- rbind(prey.bear.mutbd.out, prey.bob.mutbd.out, prey.coug.mutbd.out, prey.coy.mutbd.out, prey.wolf.mutbd.out) 
  prey.pred.mu.tbd.tbl <- dplyr::select(prey.pred.mu.tbd.out, -c("lci", "uci"))
  # write.csv(prey.pred.mu.tbd.tbl, "./Outputs/TimeBtwnDetections/Tables/tbd.prey.pred_meanTBD_table.csv")
  
  conspif.mu.tbd.out <- rbind(con.elk.mutbd.out, con.moose.mutbd.out, con.md.mutbd.out, con.wtd.mutbd.out)
  conspif.mu.tbd.tbl <- dplyr::select(conspif.mu.tbd.out, -c("lci", "uci"))
  # write.csv(conspif.mu.tbd.tbl, "./Outputs/TimeBtwnDetections/Tables/tbd.conspecific_meanTBD_table.csv")
  
  
  ####  Plot mean TBD  ####
  #'  ------------------
  #'  Effect of predator, ungulate, and conspecific on tbd
  #'  Filter to means of interest
  mean_conspif <- conspif.mu.tbd.out[conspif.mu.tbd.out$Parameter == "Mean TBD",] %>%
    dplyr::select(-c("95% CI")) %>%
    #'  Rename Parameters so they are species-specific
    # mutate(Parameter = paste0(Parameter, ": ", Species))
    # mutate(Parameter = paste0(Parameter, ": Conspecific"))
    mutate(Parameter = "Conspecific")
  mean_bob <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Bobcat",]
  mean_coy <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Coyote",]
  mean_bear <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Black bear",]
  mean_coug <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Cougar",]
  mean_wolf <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Wolf",]
  mean_pred <- rbind(mean_bob, mean_coy, mean_bear, mean_coug, mean_wolf) %>%
    dplyr::select(-c("95% CI")) %>%
    #'  Reduce parameter name to just the interacting species
    mutate(Parameter = gsub("Mean TBD: *", "", Parameter))
  mean_ungulate <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD",] %>%
    mutate(Parameter = "Mean Ungulate")
  mean_elk <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD: Elk",]
  mean_moose <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD: Moose",]
  mean_md <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD: Mule Deer",]
  mean_wtd <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD: White-tailed Deer",]
  mean_ung <- rbind(mean_ungulate, mean_elk, mean_moose, mean_md, mean_wtd) %>%  
    dplyr::select(-c("95% CI")) %>%
    #'  Reduce parameter name to just the interacting species
    mutate(Parameter = gsub("Mean TBD: *", "", Parameter))
  mean_elkprey <- prey.pred.mu.tbd.out[prey.pred.mu.tbd.out$Parameter == "Mean TBD: Elk",]
  mean_mooseprey <- prey.pred.mu.tbd.out[prey.pred.mu.tbd.out$Parameter == "Mean TBD: Moose",]
  mean_mdprey <- prey.pred.mu.tbd.out[prey.pred.mu.tbd.out$Parameter == "Mean TBD: Mule deer",]
  mean_wtdprey <- prey.pred.mu.tbd.out[prey.pred.mu.tbd.out$Parameter == "Mean TBD: White-tailed deer",]
  mean_prey <- rbind(mean_elkprey, mean_mooseprey, mean_mdprey, mean_wtdprey) %>%
    dplyr::select(-c("95% CI")) %>%
    #'  Reduce parameter name to just the interacting species
    mutate(Parameter = gsub("Mean TBD: *", "", Parameter))
  
  mean_TBD <- rbind(mean_conspif, mean_ung, mean_pred) %>%
    arrange(Species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Elk", "Moose", "Mule deer", "White-tailed deer")),
           Parameter = factor(Parameter, levels = c("Conspecific", "Mean Ungulate", "Elk", "Moose", "Mule Deer", "White-tailed Deer", "Bobcat", "Coyote", "Black bear", "Cougar", "Wolf")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  mean_TBD_v2 <- rbind(mean_pred, mean_prey) %>%
    arrange(Species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Elk", "Moose", "Mule deer", "White-tailed deer", "Bobcat", "Coyote", "Black bear", "Cougar", "Wolf")),
           Parameter = factor(Parameter, levels = c("Elk", "Moose", "Mule deer", "White-tailed deer", "Bobcat", "Coyote", "Black bear", "Cougar", "Wolf")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci))

  #'  Choose colorblind-friendly scheme
  bright <- colour("bright")
  bright(6)
  
  #'  Effect of predator, conspecific, and ungulate species on latency
  mean_TBD_skinny <- dplyr::filter(mean_TBD, mean_TBD$Parameter != "Mean Ungulate")
  mean_tbd_plot <- ggplot(mean_TBD_skinny, aes(x = Parameter, y = Estimate, group = Species)) +
    geom_errorbar(aes(ymin = lci, ymax = uci, color = Parameter), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = Parameter), size = 2.5, position = position_dodge(width = 0.4)) +
    # scale_colour_bright() +
    theme_bw() +
    ylim(0,10000) +
    facet_wrap(~Species, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Species detected prior to ungulate detection") +
    ylab("Mean number of minutes between detections") +
    ggtitle("Mean time between detections of interacting species")
  mean_tbd_plot
  # ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_SpeciesID_plot.tiff", mean_tbd_plot,
  #        units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Effect of ungulates (averaged over all ungulate species), conspecifics, and
  #'  predators on latency of site use
  mean_TBD_skinnier <- mean_TBD %>%
    filter(Parameter != "Elk") %>%
    filter(Parameter != "Moose") %>%
    filter(Parameter != "Mule Deer") %>%
    filter(Parameter != "White-tailed Deer")
  mean_tbd_plot2 <- ggplot(mean_TBD_skinnier, aes(x = Parameter, y = Estimate, group = Species)) +
    geom_errorbar(aes(ymin = lci, ymax = uci, color = Parameter), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = Parameter), size = 2.5, position = position_dodge(width = 0.4)) +
    # scale_colour_bright() +
    theme_bw() +
    ylim(0,10000) +
    facet_wrap(~Species, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Species detected prior to ungulate detection") +
    ylab("Mean number of minutes between detections") +
    ggtitle("Mean time between detections of interacting species")
  mean_tbd_plot2
  # ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_SpeciesID_plot2.tiff", mean_tbd_plot2,
  #        units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Effect of predator on ungulate latency vs ungulate prey on predator latency
  species_tbd_plot <- ggplot(mean_TBD_v2, aes(x = Parameter, y = Estimate, group = Species)) +
    geom_errorbar(aes(ymin = lci, ymax = uci, color = Parameter), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = Parameter), size = 2.5, position = position_dodge(width = 0.4)) +
    # scale_colour_bright() +
    theme_bw() +
    ylim(0,10000) +
    facet_wrap(~Species, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Species detected prior to ungulate detection") +
    ylab("Mean number of minutes between detections") +
    ggtitle("Mean time between detections of interacting species: \nPredator-prey & prey-predator")
  species_tbd_plot
  # ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_pred&prey_plot.tiff", species_tbd_plot,
  #        units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Effect of predator species vs conspecific on tbd
  #'  Filter to means of interest
  season_conspif <- conspif.mu.tbd.out[conspif.mu.tbd.out$Parameter != "Mean TBD",] %>%
    dplyr::select(-c("95% CI")) %>%
    #'  Rename Parameters so they are species-specific
    mutate(Season = gsub("Mean TBD: *", "", Parameter), 
           Parameter = "Conspecific") %>%
    rename(Interacting_spp = Parameter)
  mean_smr <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Summer",]
  mean_fall <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Fall",]
  mean_wtr <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Winter",]
  mean_sprg <- pred.prey.mu.tbd.out[pred.prey.mu.tbd.out$Parameter == "Mean TBD: Spring",]
  season_pred <- rbind(mean_smr, mean_fall, mean_wtr, mean_sprg) %>%
    dplyr::select(-c("95% CI")) %>%
    #'  Reduce parameter name to just the interacting species
    mutate(Season = gsub("Mean TBD: *", "", Parameter),
           Parameter = "Predator") %>%
    rename(Interacting_spp = Parameter) 
  mean_ung_smr <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD: Summer",]
  mean_ung_fall <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD: Fall",]
  mean_ung_wtr <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD: Winter",]
  mean_ung_sprg <- ungulate.mu.tbd.out[ungulate.mu.tbd.out$Parameter == "Mean TBD: Spring",]
  season_ung <- rbind(mean_ung_smr, mean_ung_fall, mean_ung_wtr, mean_ung_sprg) %>%
    dplyr::select(-c("95% CI")) %>%
    #'  Reduce parameter name to just the interacting species
    mutate(Season = gsub("Mean TBD: *", "", Parameter),
           Parameter = "Ungulate") %>%
    rename(Interacting_spp = Parameter) 
  mean_prey_smr <- prey.pred.mu.tbd.out[prey.pred.mu.tbd.out$Parameter == "Mean TBD: Summer",]
  mean_prey_fall <- prey.pred.mu.tbd.out[prey.pred.mu.tbd.out$Parameter == "Mean TBD: Fall",]
  mean_prey_wtr <- prey.pred.mu.tbd.out[prey.pred.mu.tbd.out$Parameter == "Mean TBD: Winter",]
  mean_prey_sprg <- prey.pred.mu.tbd.out[prey.pred.mu.tbd.out$Parameter == "Mean TBD: Spring",]
  season_prey <- rbind(mean_prey_smr, mean_prey_fall, mean_prey_wtr, mean_prey_sprg) %>%
    dplyr::select(-c("95% CI")) %>%
    #'  Reduce parameter name to just the interacting species
    mutate(Season = gsub("Mean TBD: *", "", Parameter),
           Parameter = "Prey") %>%
    rename(Interacting_spp = Parameter) 
  season_TBD <- rbind(season_conspif, season_ung, season_pred) %>%
    arrange(Species, Season, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Elk", "Moose", "Mule deer", "White-tailed deer")),
           Interacting_spp = factor(Interacting_spp, levels = c("Conspecific", "Ungulate", "Predator")),
           Season = factor(Season, levels = c("Summer", "Fall", "Winter", "Spring")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  season_TBD2 <- rbind(season_pred, season_prey) %>%
    arrange(Species, Season, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Elk", "Moose", "Mule deer", "White-tailed deer", "Bobcat", "Coyote", "Black bear", "Cougar", "Wolf")),
           Interacting_spp = factor(Interacting_spp, levels = c("Predator", "Prey")),
           Season = factor(Season, levels = c("Summer", "Fall", "Winter", "Spring")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
 
  season_tbd_plot <- ggplot(season_TBD, aes(x = Season, y = Estimate, group = Interacting_spp)) +
    geom_errorbar(aes(ymin = lci, ymax = uci, color = Interacting_spp), width = 0.2, position = position_dodge(0.4)) +
    geom_point(stat = 'identity', aes(col = Interacting_spp), size = 2.5, position = position_dodge(0.4)) +
    scale_colour_bright() +
    theme_bw() +
    ylim(0,12500) +
    facet_wrap(~Species, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(color = "Interacting species") +
    # xlab("Species detected prior to ungulate detection") +
    ylab("Mean number of minutes between detections") +
    ggtitle("Seasonal mean time between detections of interacting species")
  season_tbd_plot
  # ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_Season_plot.tiff", season_tbd_plot,
  #        units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  season_predprey_tbd_plot <- ggplot(season_TBD2, aes(x = Season, y = Estimate, group = Interacting_spp)) +
    geom_errorbar(aes(ymin = lci, ymax = uci, color = Interacting_spp), width = 0.2, position = position_dodge(0.4)) +
    geom_point(stat = 'identity', aes(col = Interacting_spp), size = 2.5, position = position_dodge(0.4)) +
    scale_colour_bright() +
    theme_bw() +
    ylim(0,12500) +
    facet_wrap(~Species, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(color = "Interacting species") +
    xlab("Species detected prior to focal species") +
    ylab("Mean number of minutes between detections") +
    ggtitle("Seasonal mean time between detections of interacting species: \nPredator-prey & prey-predator")
  season_predprey_tbd_plot
  # ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_Season_pred&prey_plot.tiff", season_predprey_tbd_plot,
  #        units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')

  
  #' #'  Load tbd data with covariates
  #' load("./Outputs/tbd_pred.prey_2022-09-23.RData")
  #' #'  Generate range of TRI and % Forest values to plot across
  #' #'  Generate range of TRI & % Forest covariate values to predict across based 
  #' #'  on range of values across all camera sites
  #' scaledTRI <- scale(tbd_pred.prey$TRI); scaledFor <- scale(tbd_pred.prey$PercForest)
  #' print(minmax_tri <- range(scaledTRI)); print(minmax_for <- range(scaledFor))
  #' newTRI <- seq(from = minmax_tri[1], to = minmax_tri[2], length.out = 100)
  #' newFor <- seq(from = minmax_for[1], to = minmax_for[2], length.out = 100)
  #' newcovs <- as.data.frame(cbind(newTRI, newFor)) %>%
  #'   mutate(backTRI = (newTRI * attr(scaledTRI, 'scaled:scale')) + attr(scaledTRI, 'scaled:center'),
  #'          backFor = (newFor * attr(scaledFor, 'scaled:scale')) + attr(scaledFor, 'scaled:center'))

  #'  Grab scaled TRI and % Forest from bundled data so specific to each analysis
  load("./Data/pred_md_bundled.RData"); pred.md_newcovs <- pred_md_bundled[[7]]
  load("./Data/pred_elk_bundled.RData"); pred.elk_newcovs <- pred_elk_bundled[[7]]
  load("./Data/pred_moose_bundled.RData"); pred.moose_newcovs <- pred_moose_bundled[[7]]
  load("./Data/pred_wtd_bundled.RData"); pred.wtd_newcovs <- pred_wtd_bundled[[7]]
  
  load("./Data/con_md_bundled.RData"); con.md_newcovs <- con_md_bundled[[7]]
  load("./Data/con_elk_bundled.RData"); con.elk_newcovs <- con_elk_bundled[[7]]
  load("./Data/con_moose_bundled.RData"); con.moose_newcovs <- con_moose_bundled[[7]]
  load("./Data/con_wtd_bundled.RData"); con.wtd_newcovs <- con_wtd_bundled[[7]]

  #'  Append scaled covariate data to each predicted TBD value - PREDATOR-PREY ANALYSIS
  predicted_dat <- function(newcovs, pred.out) {
    #'  Scaled TRI & % Forest
    newTRI <- newcovs[,1]
    newFor <- newcovs[,2]
    
    #'  Predator-specific mean TBD values predicted across scaled TRI values
    tbd_bob_tri <- cbind(pred.out[30:129,], newTRI) %>% mutate(Predator = "Bobcat")
    tbd_coy_tri <- cbind(pred.out[130:229,], newTRI) %>% mutate(Predator = "Coyote")
    tbd_bear_tri <- cbind(pred.out[230:329,], newTRI) %>% mutate(Predator = "Black bear")
    tbd_coug_tri <- cbind(pred.out[330:429,], newTRI) %>% mutate(Predator = "Cougar")
    tbd_wolf_tri <- cbind(pred.out[430:529,], newTRI) %>% mutate(Predator = "Wolf")
    tbd_tri <- rbind(tbd_bob_tri, tbd_coy_tri, tbd_bear_tri, tbd_coug_tri, tbd_wolf_tri)
    
    #'  Predator-specific mean TBD predicted across scaled % forest values
    tbd_bob_for <- cbind(pred.out[530:629,], newFor) %>% mutate(Predator = "Bobcat")
    tbd_coy_for <- cbind(pred.out[630:729,], newFor) %>% mutate(Predator = "Coyote")
    tbd_bear_for <- cbind(pred.out[730:829,], newFor) %>% mutate(Predator = "Black bear")
    tbd_coug_for <- cbind(pred.out[830:929,], newFor) %>% mutate(Predator = "Cougar")
    tbd_wolf_for <- cbind(pred.out[930:1029,], newFor) %>% mutate(Predator = "Wolf")
    tbd_for <- rbind(tbd_bob_for, tbd_coy_for, tbd_bear_for, tbd_coug_for, tbd_wolf_for)
    
    tbd_list <- list(tbd_tri, tbd_for)
    return(tbd_list)
  }
  predicted_pred.md <- predicted_dat(pred.md_newcovs, pred.md.out)
  predicted_pred.moose <- predicted_dat(pred.moose_newcovs, pred.moose.out)
  predicted_pred.wtd <- predicted_dat(pred.wtd_newcovs, pred.wtd.out)
  
  #'  Append scaled covariate data to predicted ELK TBD values - PREDATOR-PREY ANALYSIS
  #'  Note: different b/c no wolves or interactions in this model
  #'  Scaled TRI & % Forest
  predicted_elk <- function(newcovs, pred.out) {
    newTRI <- newcovs[,1]
    newFor <- newcovs[,2]
    
    #'  Predator-specific mean TBD values predicted across scaled TRI values
    tbd_bob_tri <- cbind(pred.out[20:119,], newTRI) %>% mutate(Predator = "Bobcat")
    tbd_coy_tri <- cbind(pred.out[120:219,], newTRI) %>% mutate(Predator = "Coyote")
    tbd_bear_tri <- cbind(pred.out[220:319,], newTRI) %>% mutate(Predator = "Black bear")
    tbd_coug_tri <- cbind(pred.out[320:419,], newTRI) %>% mutate(Predator = "Cougar")
    tbd_tri <- rbind(tbd_bob_tri, tbd_coy_tri, tbd_bear_tri, tbd_coug_tri)
    
    #'  Predator-specific mean TBD predicted across scaled % forest values
    tbd_bob_for <- cbind(pred.out[420:519,], newFor) %>% mutate(Predator = "Bobcat")
    tbd_coy_for <- cbind(pred.out[520:619,], newFor) %>% mutate(Predator = "Coyote")
    tbd_bear_for <- cbind(pred.out[620:719,], newFor) %>% mutate(Predator = "Black bear")
    tbd_coug_for <- cbind(pred.out[720:819,], newFor) %>% mutate(Predator = "Cougar")
    tbd_for <- rbind(tbd_bob_for, tbd_coy_for, tbd_bear_for, tbd_coug_for)
    
    tbd_list <- list(tbd_tri, tbd_for)
    return(tbd_list)
  }
  predicted_pred.elk <- predicted_elk(pred.elk_newcovs, pred.elk.out)
  
  #'  Append scaled covariate data to each predicted TBD value - CONSPECIFIC ANALYSIS
  predicted_dat_con <- function(newcovs, con.out, conspif) {
    #'  Scaled TRI & % Forest
    newTRI <- newcovs[,1]
    newFor <- newcovs[,2]
    
    #'  Predator-specific mean TBD values predicted across scaled TRI values
    tbd_tri <- cbind(con.out[13:112,], newTRI) %>% mutate(Predator = conspif)
    
    #'  Predator-specific mean TBD predicted across scaled % forest values
    tbd_for <- cbind(con.out[113:212,], newFor) %>% mutate(Predator = conspif)
    
    tbd_list <- list(tbd_tri, tbd_for)
    return(tbd_list)
  }
  predicted_con.md <- predicted_dat_con(con.md_newcovs, con.md.out, conspif = "Mule deer")
  predicted_con.elk <- predicted_dat_con(con.elk_newcovs, con.elk.out, conspif = "Elk")
  predicted_con.moose <- predicted_dat_con(con.moose_newcovs, con.moose.out, conspif = "Moose")
  predicted_con.wtd <- predicted_dat_con(con.wtd_newcovs, con.wtd.out, conspif = "White-tailed deer")
  
  #' #'  Merge all together
  #' predicted_md <- mapply(rbind, predicted_pred.md, predicted_con.md, SIMPLIFY = FALSE) 
  #' predicted_elk <- mapply(rbind, predicted_pred.elk, predicted_con.elk, SIMPLIFY = FALSE) 
  #' predicted_moose <- mapply(rbind, predicted_pred.moose, predicted_con.moose, SIMPLIFY = FALSE) 
  #' predicted_wtd <- mapply(rbind, predicted_pred.wtd, predicted_con.wtd, SIMPLIFY = FALSE) 
  #' 
  # predicted_tri <- rbind(predicted_pred.md[[1]], predicted_pred.elk[[1]], predicted_pred.moose[[1]], predicted_pred.wtd[[1]],
  #                             predicted_con.md[[1]], predicted_con.elk[[1]], predicted_con.moose[[1]], predicted_con.wtd[[1]]) %>%
  #   mutate(Predator = factor(Predator, levels = c("Bobcat", "Coyote", "Black bear", "Cougar", "Wolf", "Mule deer", "Elk", "Moose", "White-tailed deer")))
  # predicted_for <- rbind(predicted_pred.md[[2]], predicted_pred.elk[[2]], predicted_pred.moose[[2]], predicted_pred.wtd[[2]],
  #                             predicted_con.md[[2]], predicted_con.elk[[2]], predicted_con.moose[[2]], predicted_con.wtd[[2]]) %>%
  #   mutate(Predator = factor(Predator, levels = c("Bobcat", "Coyote", "Black bear", "Cougar", "Wolf", "Mule deer", "Elk", "Moose", "White-tailed deer")))
  # predicted_pred <- list(predicted_tri, predicted_for)
  
  #'  ---------------------------------------------------
  ####  Plot interaction between predators and habitat  ####
  #'  ---------------------------------------------------
  #'  Mule deer
  md_TRI_plot <- ggplot(NULL, aes(x = newTRI, y = Estimate, colour = Predator)) + 
    geom_line(data  = predicted_pred.md[[1]], size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(data = predicted_pred.md[[1]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Add conspecific data to plot
    geom_line(data = predicted_con.md[[1]], size = 0.75, lty = "dashed", col = "black") +
    geom_ribbon(data = predicted_con.md[[1]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA, fill = "black") +
    # scale_color_manual(values=c("red", "blue", "green", "yellow", "orange", "black")) + 
    # scale_fill_manual(values=c("red", "blue", "green", "yellow", "orange", "black")) + 
    #'  Get rid of lines and gray backgroundtheme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    #theme(legend.position="bottom") +
    coord_cartesian(ylim = c(0, 10000)) +
    xlab("Scaled terrain ruggedness") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of terrain ruggedness and predators on mule deer latency")
  md_TRI_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_md_TRI_plot.tiff", md_TRI_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  md_For_plot <- ggplot(NULL, aes(x = newFor, y = Estimate, colour = Predator)) + 
    geom_line(data = predicted_pred.md[[2]], size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(data = predicted_pred.md[[2]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Add conspecific data to plot
    geom_line(data = predicted_con.md[[2]], size = 0.75, lty = "dashed", col = "black") +
    geom_ribbon(data = predicted_con.md[[2]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA, fill = "black") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    coord_cartesian(ylim = c(0, 10000), xlim = c(-.5, 2.5)) +
    xlab("Scaled percent forested habitat") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of forested habitat and predators on mule deer latency")
  md_For_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_md_PercForest_plot.tiff", md_For_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  Elk
  elk_TRI_plot <- ggplot(NULL, aes(x = newTRI, y = Estimate, colour = Predator)) + 
    geom_line(data = predicted_pred.elk[[1]], size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(data = predicted_pred.elk[[1]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Add conspecific data to plot
    geom_line(data = predicted_con.elk[[1]], size = 0.75, lty = "dashed", col = "black") +
    geom_ribbon(data = predicted_con.elk[[1]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA, fill = "black") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    coord_cartesian(ylim = c(0, 25000), xlim = c(-1.5, 3.5)) +
    xlab("Scaled terrain ruggedness") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of terrain ruggedness and predators on elk latency")
  elk_TRI_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_elk_TRI_plot.tiff", elk_TRI_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  elk_For_plot <- ggplot(NULL, aes(x = newFor, y = Estimate, colour = Predator)) + 
    geom_line(data = predicted_pred.elk[[2]], size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(data = predicted_pred.elk[[2]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Add conspecific data to plot
    geom_line(data = predicted_con.elk[[2]], size = 0.75, lty = "dashed", col = "black") +
    geom_ribbon(data = predicted_con.elk[[2]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA, fill = "black") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    coord_cartesian(ylim = c(0, 10000), xlim = c(-.5, 2.5)) +
    xlab("Scaled percent forest habitat") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of forested habitat and predators on elk latency")
  elk_For_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_elk_PercForest_plot.tiff", elk_For_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Moose
  moose_TRI_plot <- ggplot(NULL, aes(x = newTRI, y = Estimate, colour = Predator)) + 
    geom_line(data = predicted_pred.moose[[1]], size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(data = predicted_pred.moose[[1]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Add conspecific data to plot
    geom_line(data = predicted_con.moose[[1]], size = 0.75, lty = "dashed", col = "black") +
    geom_ribbon(data = predicted_con.moose[[1]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA, fill = "black") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    coord_cartesian(ylim = c(0, 25000), xlim = c(-1.5, 2.5)) +
    xlab("Scaled terrain ruggedness") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of terrain ruggedness and predators on moose latency")
  moose_TRI_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_moose_TRI_plot.tiff", moose_TRI_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  moose_For_plot <- ggplot(NULL, aes(x = newFor, y = Estimate, colour = Predator)) + 
    geom_line(data = predicted_pred.moose[[2]], size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(data = predicted_pred.moose[[2]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Add conspecific data to plot
    geom_line(data = predicted_con.moose[[2]], size = 0.75, lty = "dashed", col = "black") +
    geom_ribbon(data = predicted_con.moose[[2]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA, fill = "black") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    coord_cartesian(ylim = c(0, 30000), xlim = c(-.5, 2.5)) +
    xlab("Scaled percent forest habitat") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of forested habitat and predators on moose latency")
  moose_For_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_moose_PercForest_plot.tiff", moose_For_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  White-tailed deer
  wtd_TRI_plot <- ggplot(NULL, aes(x = newTRI, y = Estimate, colour = Predator)) + 
    geom_line(data = predicted_pred.wtd[[1]], size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(data = predicted_pred.wtd[[1]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Add conspecific data to plot
    geom_line(data = predicted_con.wtd[[1]], size = 0.75, lty = "dashed", col = "black") +
    geom_ribbon(data = predicted_con.wtd[[1]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA, fill = "black") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    coord_cartesian(ylim = c(0, 15000), xlim = c(-1.5, 3)) +
    xlab("Scaled terrain ruggedness") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of terrain ruggedness and predators on white-tailed deer latency")
  wtd_TRI_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_wtd_TRI_plot.tiff", wtd_TRI_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  wtd_For_plot <- ggplot(NULL, aes(x = newFor, y = Estimate, colour = Predator)) + 
    geom_line(data = predicted_pred.wtd[[2]], size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(data = predicted_pred.wtd[[2]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Add conspecific data to plot
    geom_line(data = predicted_con.wtd[[2]], size = 0.75, lty = "dashed", col = "black") +
    geom_ribbon(data = predicted_con.wtd[[2]], aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA, fill = "black") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    coord_cartesian(xlim = c(-0.5, 2.5)) +
    xlab("Scaled percent forest habitat") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of forested habitat and predators on white-tailed deer latency")
  wtd_For_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_wtd_PercForest_plot.tiff", wtd_For_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
 
  
  