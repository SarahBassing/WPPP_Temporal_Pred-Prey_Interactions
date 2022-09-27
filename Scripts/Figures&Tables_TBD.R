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
  load("./Outputs/TimeBtwnDetections/tbd.pred.md-season_predID_X_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.pred.elk-season_predID_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.pred.moose-season_predID_X_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.pred.wtd-season_predID_X_habitat.RData")
  
  load("./Outputs/TimeBtwnDetections/tbd.md-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.elk-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.moose-season_habitat.RData")
  load("./Outputs/TimeBtwnDetections/tbd.wtd-season_habitat.RData")

  #'  Pull out and rename coefficient estimates 
  coefs <- function(mod_out, spp) {
    Species <- spp
    Estimate <- unlist(mod_out$mean)
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
             Parameter = ifelse(Parameter == "beta31", "Predator: Bobcat TRI", Parameter),
             Parameter = ifelse(Parameter == "beta32", "Predator: Coyote TRI", Parameter),
             Parameter = ifelse(Parameter == "beta33", "Predator: Black bear TRI", Parameter),
             Parameter = ifelse(Parameter == "beta34", "Predator: Cougar TRI", Parameter),
             Parameter = ifelse(Parameter == "beta35", "Predator: Wolf TRI", Parameter),
             Parameter = ifelse(Parameter == "beta41", "Predator: Bobcat Forest", Parameter),
             Parameter = ifelse(Parameter == "beta42", "Predator: Coyote Forest", Parameter),
             Parameter = ifelse(Parameter == "beta43", "Predator: Black bear Forest", Parameter),
             Parameter = ifelse(Parameter == "beta44", "Predator: Cougar Forest", Parameter),
             Parameter = ifelse(Parameter == "beta45", "Predator: Wolf Forest", Parameter),
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
  
  #'  Save coefficient estimates only
  pred.md.coef.out <- pred.md.out[1:18,1:4]
  pred.elk.coef.out <- pred.elk.out[1:9,1:4]
  pred.moose.coef.out <- pred.moose.out[1:18,1:4]
  pred.wtd.coef.out <- pred.wtd.out[1:18,1:4]
  
  con.md.coef.out <- con.md.out[1:6,1:4]
  con.elk.coef.out <- con.elk.out[1:6,1:4]  
  con.moose.coef.out <- con.moose.out[1:6,1:4]
  con.wtd.coef.out <- con.wtd.out[1:6,1:4]
  
  pred.prey.coef.out <- rbind(pred.elk.coef.out, pred.moose.coef.out, pred.md.coef.out, pred.wtd.coef.out)
  # write.csv(pred.prey.coef.out, "./Outputs/TimeBtwnDetections/Tables/tbd.pred.prey_coef_table.csv")
  
  conspif.coef.out <- rbind(con.elk.coef.out, con.moose.coef.out, con.md.coef.out, con.wtd.coef.out)
  # write.csv(conspif.coef.out, "./Outputs/TimeBtwnDetections/Tables/tbd.conspecific_coef_table.csv")
  
  #'  Save mean time-between-detection results only
  pred.md.mutbd.out <- pred.md.out[20:29,1:6]
  pred.elk.mutbd.out <- pred.elk.out[11:19,1:6]
  pred.moose.mutbd.out <- pred.moose.out[20:29,1:6]
  pred.wtd.mutbd.out <- pred.wtd.out[20:29,1:6]
  
  con.md.mutbd.out <- con.md.out[8:12,1:6]
  con.wtd.mutbd.out <- con.wtd.out[8:12,1:6]
  con.moose.mutbd.out <- con.moose.out[8:12,1:6]
  con.elk.mutbd.out <- con.elk.out[8:12,1:6]
  
  pred.prey.mu.tbd.out <- rbind(pred.elk.mutbd.out, pred.moose.mutbd.out, pred.md.mutbd.out, pred.wtd.mutbd.out) 
  pred.prey.mu.tbd.tbl <- dplyr::select(pred.prey.mu.tbd.out, -c("lci", "uci"))
  # write.csv(pred.prey.mu.tbd.tbl, "./Outputs/TimeBtwnDetections/Tables/tbd.pred.prey_meanTBD_table.csv")
  
  conspif.mu.tbd.out <- rbind(con.elk.mutbd.out, con.moose.mutbd.out, con.md.mutbd.out, con.wtd.mutbd.out)
  conspif.mu.tbd.tbl <- dplyr::select(conspif.mu.tbd.out, -c("lci", "uci"))
  # write.csv(conspif.mu.tbd.tbl, "./Outputs/TimeBtwnDetections/Tables/tbd.conspecific_meanTBD_table.csv")
  
  
  ####  Plot mean TBD  ####
  #'  ------------------
  #'  Effect of predator species vs conspecific on tbd
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
  mean_TBD <- rbind(mean_conspif, mean_pred) %>%
    arrange(Species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Elk", "Moose", "Mule deer", "White-tailed deer")),
           Parameter = factor(Parameter, levels = c("Conspecific", "Bobcat", "Coyote", "Black bear", "Cougar", "Wolf")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci))

  #'  Choose colorblind-friendly scheme
  bright <- colour("bright")
  bright(6)
  
  mean_tbd_plot <- ggplot(mean_TBD, aes(x = Parameter, y = Estimate, group = Species)) +
    geom_errorbar(aes(ymin = lci, ymax = uci, color = Parameter), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = Parameter), size = 2.5, position = position_dodge(width = 0.4)) +
    scale_colour_bright() +
    theme_bw() +
    ylim(0,10000) +
    facet_wrap(~Species, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Species detected prior to ungulate detection") +
    ylab("Mean number of minutes between detections") +
    ggtitle("Mean time between detections of interacting species")
  mean_tbd_plot
  # ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_PredatorID_plot.tiff", mean_tbd_plot, 
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
  season_TBD <- rbind(season_conspif, season_pred) %>%
    arrange(Species, Season, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Elk", "Moose", "Mule deer", "White-tailed deer")),
           Interacting_spp = factor(Interacting_spp, levels = c("Conspecific", "Predator")),
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

  #'  Function to append scaled covariate data to each predicted TBD value
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
  
  #'  Append scaled covariate data to predicted ELK TBD values (different b/c no wolves or interactions)
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
  
  #'  Merge all together
  predicted_pred.tri <- rbind(predicted_pred.md[[1]], predicted_pred.elk[[1]], predicted_pred.moose[[1]], predicted_pred.wtd[[1]]) %>%
    mutate(Predator = factor(Predator, levels = c("Bobcat", "Coyote", "Black bear", "Cougar", "Wolf")))
  predicted_pred.for <- rbind(predicted_pred.md[[2]], predicted_pred.elk[[2]], predicted_pred.moose[[2]], predicted_pred.wtd[[2]]) %>%
    mutate(Predator = factor(Predator, levels = c("Bobcat", "Coyote", "Black bear", "Cougar", "Wolf")))
  predicted_pred <- list(predicted_pred.tri, predicted_pred.for)
  
  #'  ---------------------------------------------------
  ####  Plot interaction between predators and habitat  ####
  #'  ---------------------------------------------------
  #'  Mule deer
  md_TRI_plot <- ggplot(predicted_pred.md[[1]], aes(x = newTRI, y = Estimate, colour = Predator)) + 
    geom_line(size = 0.75) +
    # scale_color_manual(values=c("#40B0A6", "#E66100")) +  #, "#5D3A9B"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    # scale_fill_manual(values=c("#40B0A6", "#E66100")) + #, "#5D3A9B"
    #'  Get rid of lines and gray background
    theme_bw() +
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
  
  md_For_plot <- ggplot(predicted_pred.md[[2]], aes(x = newFor, y = Estimate, colour = Predator)) + 
    geom_line(size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    # coord_cartesian(ylim = c(0, 10000)) +
    xlab("Scaled percent forested habitat") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of forested habitat and predators on mule deer latency")
  md_For_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_md_PercForest_plot.tiff", md_For_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  Elk
  elk_TRI_plot <- ggplot(predicted_pred.elk[[1]], aes(x = newTRI, y = Estimate, colour = Predator)) + 
    geom_line(size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    # coord_cartesian(ylim = c(0, 10000)) +
    xlab("Scaled terrain ruggedness") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of terrain ruggedness and predators on elk latency")
  elk_TRI_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_elk_TRI_plot.tiff", elk_TRI_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  elk_For_plot <- ggplot(predicted_pred.elk[[2]], aes(x = newFor, y = Estimate, colour = Predator)) + 
    geom_line(size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    # coord_cartesian(ylim = c(0, 10000)) +
    xlab("Scaled percent forest habitat") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of forested habitat and predators on elk latency")
  elk_For_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_elk_PercForest_plot.tiff", elk_For_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Moose
  moose_TRI_plot <- ggplot(predicted_pred.moose[[1]], aes(x = newTRI, y = Estimate, colour = Predator)) + 
    geom_line(size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    # coord_cartesian(ylim = c(0, 10000)) +
    xlab("Scaled terrain ruggedness") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of terrain ruggedness and predators on moose latency")
  moose_TRI_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_moose_TRI_plot.tiff", moose_TRI_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  moose_For_plot <- ggplot(predicted_pred.moose[[2]], aes(x = newFor, y = Estimate, colour = Predator)) + 
    geom_line(size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    # coord_cartesian(ylim = c(0, 10000)) +
    xlab("Scaled percent forest habitat") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of forested habitat and predators on moose latency")
  moose_For_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_moose_PercForest_plot.tiff", moose_For_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  White-tailed deer
  wtd_TRI_plot <- ggplot(predicted_pred.wtd[[1]], aes(x = newTRI, y = Estimate, colour = Predator)) + 
    geom_line(size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    coord_cartesian(ylim = c(0, 10000)) +
    xlab("Scaled terrain ruggedness") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of terrain ruggedness and predators on white-tailed deer latency")
  wtd_TRI_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_wtd_TRI_plot.tiff", wtd_TRI_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  wtd_For_plot <- ggplot(predicted_pred.wtd[[2]], aes(x = newFor, y = Estimate, colour = Predator)) + 
    geom_line(size = 0.75) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Predator), alpha = 0.3, colour = NA) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    xlab("Scaled percent forest habitat") +
    ylab("Mean time-between-detections") +
    ggtitle("Effect of forested habitat and predators on white-tailed deer latency")
  wtd_For_plot
  ggsave("./Outputs/TimeBtwnDetections/Figures/TBD_wtd_PercForest_plot.tiff", wtd_For_plot,
         units = "in", width = 6, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  
  
  
  
  
  
  
  
  
  