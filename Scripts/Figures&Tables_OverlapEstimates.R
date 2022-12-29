  #'  =====================================================
  #'  Temporal overlap analysis: results tables and figures
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  August 2022
  #'  =====================================================
  #'  Script pulls in outputs from temporal overlap analyses conducted with the
  #'  Temporal_overlap_PredPrey_HabitatComplexity.R (script that calculates overlap 
  #'  of predator-prey activity under different levels of habtiat compleixty) & 
  #'  Temporal_overlap_PredPrey_HabitatComplexity&Predators.R (script that calculates 
  #'  prey activity under varying levels of risk and predator activity under varying 
  #'  levels of habitat complexity) scripts, then makes tables & figures for publication.
  #'  =====================================================
  
  #'  Libraries
  library(ggplot2)
  library(khroma)
  library(patchwork)
  library(sp)
  library(raster)
  library(tidyverse)  

  #'  Load output from temporal overlap analysis
  load("./Outputs/Temporal Overlap/PredPrey_TRI_Overlap_2022-09-19.RData")
  load("./Outputs/Temporal Overlap/PredPrey_PercForest_Overlap_2022-09-19.RData")
  load("./Outputs/Temporal Overlap/PreyOnly_TRI_Forest_Pred_Overlap_2022-09-26.RData")
  load("./Outputs/Temporal Overlap/PredOnly_TRI_Forest_Overlap_2022-10-15.RData")
  
  
  ####  Results Tables  ####
  #'  Format data into tables for easier viewing and to create figures
  #'  ----------------------------
  #'  Predator-prey results tables
  #'  ----------------------------
  #'  Create tables from predator-prey overlap estimates in low vs high risk habitats
  #'  Indexing [2,i] gives norm0 CIs, [4,i] gives basic0 CIs
  predprey_table <- function(overlap_out, spp1, spp2, risk_type, season) {
    dhat_low <- round(overlap_out$dhat_lowrisk, 2)
    low_lci <- round(overlap_out$lowrisk_CI[4,1], 2)
    low_uci <- round(overlap_out$lowrisk_CI[4,2], 2)
    dhat_high <- round(overlap_out$dhat_highrisk, 2)
    high_lci <- round(overlap_out$highrisk_CI[4,1], 2)
    high_uci <- round(overlap_out$highrisk_CI[4,2], 2)
    pair <- paste0(spp1, "-", spp2)
    spp <- c(pair, pair)
    predator <- c(spp1, spp1)
    prey <- c(spp2, spp2)
    risk <- c("Low", "High")
    season <- c(season, season)
    Dhat <- c(dhat_low, dhat_high)
    l95 <- c(low_lci, high_lci)
    u95 <- c(low_uci, high_uci)
    ndet_predator <- c(overlap_out$ndet_spp1_lowrisk, overlap_out$ndet_spp1_highrisk)
    ndet_prey <- c(overlap_out$ndet_spp2_lowrisk, overlap_out$ndet_spp2_highrisk)
    df <- as.data.frame(cbind(spp, predator, prey, season, risk, Dhat, l95, u95, 
                              ndet_predator, ndet_prey))
    rownames(df) <- NULL
    names(df)[names(df) == "spp"] <- "Species.pair"
    names(df)[names(df) == "risk"] <- paste0(risk_type, "_level")
    df <- mutate(df, Dhat = as.numeric(Dhat),
                 l95 = as.numeric(l95),
                 u95 = as.numeric(u95))
    return(df)
  }
  #'  When using two different overlap estimators for same species pairing (e.g.,
  #'  when low vs high risk sample sizes requires different overlap estimators), 
  #'  need to filter to the appropriate estimator given sample size- 
  #'  dhat1 for when n < 50 observations, dhat4 for n > 50 observations
  #'  Note: overlap estimate when covariate value is low is always first row [1,],  
  #'  overlap estimate when covariate value is high is always second row [2,]
  ####  Cougar - Prey TRI Output  ####
  #'  Cougar-mule deer
  coug_md_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[1]][[1]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "TRI", season = "Summer")
  coug_md_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[1]][[2]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "TRI", season = "Fall")
  coug_md_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[1]][[3]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "TRI", season = "Winter")
  coug_md_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[1]][[4]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "TRI", season = "Spring")
  #'  Cougar-elk
  #' Low risk >50 cougars & elk; High risk <50 cougars & elk
  coug_elk_tri_smr_out1 <- predprey_table(pred_prey_tri_overlap[[1]][[5]], spp1 = "Cougar", spp2 = "Elk", risk_type = "TRI", season = "Summer")
  coug_elk_tri_smr_out4 <- predprey_table(pred_prey_tri_overlap[[1]][[6]], spp1 = "Cougar", spp2 = "Elk", risk_type = "TRI", season = "Summer")
  coug_elk_tri_smr_out <- rbind(coug_elk_tri_smr_out4[1,], coug_elk_tri_smr_out1[2,])
  coug_elk_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[1]][[7]], spp1 = "Cougar", spp2 = "Elk", risk_type = "TRI", season = "Fall")
  coug_elk_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[1]][[8]], spp1 = "Cougar", spp2 = "Elk", risk_type = "TRI", season = "Winter")
  coug_elk_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[1]][[9]], spp1 = "Cougar", spp2 = "Elk", risk_type = "TRI", season = "Spring")
  #'  Cougar-moose
  coug_moose_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[1]][[10]], spp1 = "Cougar", spp2 = "Moose", risk_type = "TRI", season = "Summer")
  coug_moose_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[1]][[11]], spp1 = "Cougar", spp2 = "Moose", risk_type = "TRI", season = "Fall")
  coug_moose_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[1]][[12]], spp1 = "Cougar", spp2 = "Moose", risk_type = "TRI", season = "Winter")
  coug_moose_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[1]][[13]], spp1 = "Cougar", spp2 = "Moose", risk_type = "TRI", season = "Spring")
  #'  Cougar-white-tailed deer
  coug_wtd_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[1]][[14]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Summer")
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_wtd_tri_fall_out1 <- predprey_table(pred_prey_tri_overlap[[1]][[15]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Fall")
  coug_wtd_tri_fall_out4 <- predprey_table(pred_prey_tri_overlap[[1]][[16]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Fall")
  coug_wtd_tri_fall_out <- rbind(coug_wtd_tri_fall_out4[1,], coug_wtd_tri_fall_out1[2,])
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_wtd_tri_wtr_out1 <- predprey_table(pred_prey_tri_overlap[[1]][[17]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Winter")
  coug_wtd_tri_wtr_out4 <- predprey_table(pred_prey_tri_overlap[[1]][[18]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Winter")
  coug_wtd_tri_wtr_out <- rbind(coug_wtd_tri_wtr_out4[1,], coug_wtd_tri_wtr_out1[2,])
  #'  Low risk >50 cougars; High risk < 50 cougars
  coug_wtd_tri_sprg_out1 <- predprey_table(pred_prey_tri_overlap[[1]][[19]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Spring")
  coug_wtd_tri_sprg_out4 <- predprey_table(pred_prey_tri_overlap[[1]][[20]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Spring")
  coug_wtd_tri_sprg_out <- rbind(coug_wtd_tri_sprg_out4[1,], coug_wtd_tri_sprg_out1[2,])
  #'  Cougar-prey results table
  coug_prey_tri_out <- rbind(coug_md_tri_smr_out, coug_md_tri_fall_out, coug_md_tri_wtr_out, coug_md_tri_sprg_out,
                         coug_elk_tri_smr_out, coug_elk_tri_fall_out, coug_elk_tri_wtr_out, coug_elk_tri_sprg_out,
                         coug_moose_tri_smr_out, coug_moose_tri_fall_out, coug_moose_tri_wtr_out, coug_moose_tri_sprg_out,
                         coug_wtd_tri_smr_out, coug_wtd_tri_fall_out, coug_wtd_tri_wtr_out, coug_wtd_tri_sprg_out)
  
  ####  Wolf - Prey TRI Output  ####
  #'  Wolf-mule deer
  wolf_md_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[2]][[1]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "TRI", season = "Summer")
  wolf_md_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[2]][[2]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "TRI", season = "Fall")
  wolf_md_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[2]][[3]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "TRI", season = "Winter")
  wolf_md_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[2]][[4]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "TRI", season = "Spring")
  #'  Wolf-elk
  wolf_elk_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[2]][[5]], spp1 = "Wolf", spp2 = "Elk", risk_type = "TRI", season = "Summer")
  #'  Wolf-moose
  wolf_moose_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[2]][[6]], spp1 = "Wolf", spp2 = "Moose", risk_type = "TRI", season = "Summer")
  wolf_moose_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[2]][[7]], spp1 = "Wolf", spp2 = "Moose", risk_type = "TRI", season = "Fall")
  wolf_moose_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[2]][[8]], spp1 = "Wolf", spp2 = "Moose", risk_type = "TRI", season = "Winter")
  wolf_moose_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[2]][[9]], spp1 = "Wolf", spp2 = "Moose", risk_type = "TRI", season = "Spring")
  #'  Wolf-white-tailed deer
  wolf_wtd_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[2]][[10]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Summer")
  wolf_wtd_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[2]][[11]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Fall")
  wolf_wtd_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[2]][[12]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Winter")
  wolf_wtd_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[2]][[13]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Spring")
  #'  Wolf-prey results table
  wolf_prey_tri_out <- rbind(wolf_md_tri_smr_out, wolf_md_tri_fall_out, wolf_md_tri_wtr_out, wolf_md_tri_sprg_out, wolf_elk_tri_smr_out,
                         wolf_moose_tri_smr_out, wolf_moose_tri_fall_out, wolf_moose_tri_wtr_out,wolf_moose_tri_sprg_out,
                         wolf_wtd_tri_smr_out, wolf_wtd_tri_fall_out, wolf_wtd_tri_wtr_out, wolf_wtd_tri_sprg_out)
  
  ####  Black bear - Prey TRI Output  ####
  #'  Black bear-mule deer
  bear_md_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[3]][[1]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "TRI", season = "Summer")
  #'  Low risk 50 black bears; High risk <50 black bears
  bear_md_tri_fall_out1 <- predprey_table(pred_prey_tri_overlap[[3]][[2]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "TRI", season = "Fall")
  bear_md_tri_fall_out4 <- predprey_table(pred_prey_tri_overlap[[3]][[3]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "TRI", season = "Fall")
  bear_md_tri_fall_out <- rbind(bear_md_tri_fall_out4[1,], bear_md_tri_fall_out1[2,])
  bear_md_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[3]][[4]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "TRI", season = "Spring")
  #'  Black bear-elk
  bear_elk_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[3]][[5]], spp1 = "Black bear", spp2 = "Elk", risk_type = "TRI", season = "Summer")
  bear_elk_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[3]][[6]], spp1 = "Black bear", spp2 = "Elk", risk_type = "TRI", season = "Fall")
  bear_elk_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[3]][[7]], spp1 = "Black bear", spp2 = "Elk", risk_type = "TRI", season = "Spring")
  #'  Black bear-moose
  bear_moose_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[3]][[8]], spp1 = "Black bear", spp2 = "Moose", risk_type = "TRI", season = "Summer")
  bear_moose_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[3]][[9]], spp1 = "Black bear", spp2 = "Moose", risk_type = "TRI", season = "Fall")
  bear_moose_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[3]][[10]], spp1 = "Black bear", spp2 = "Moose", risk_type = "TRI", season = "Spring")
  #'  Black bear-white-tailed deer
  bear_wtd_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[3]][[11]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Summer")
  bear_wtd_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[3]][[12]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Fall")
  #'  Low risk >50 black bears; High risk <50 black bears
  bear_wtd_tri_sprg_out1 <- predprey_table(pred_prey_tri_overlap[[3]][[13]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Spring")
  bear_wtd_tri_sprg_out4 <- predprey_table(pred_prey_tri_overlap[[3]][[14]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Spring")
  bear_wtd_tri_sprg_out <- rbind(bear_wtd_tri_sprg_out4[1,], bear_wtd_tri_sprg_out1[2,])
  #'  Black bear-prey results table
  bear_prey_tri_out <- rbind(bear_md_tri_smr_out, bear_md_tri_fall_out, bear_md_tri_sprg_out,
                         bear_elk_tri_smr_out, bear_elk_tri_fall_out, bear_elk_tri_sprg_out,
                         bear_moose_tri_smr_out, bear_moose_tri_fall_out, bear_moose_tri_sprg_out,
                         bear_wtd_tri_smr_out, bear_wtd_tri_fall_out, bear_wtd_tri_sprg_out)
  
  ####  Bobcat - Prey TRI Output  ####
  #'  Bobcat-mule deer
  bob_md_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[4]][[1]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "TRI", season = "Summer")
  bob_md_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[4]][[2]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "TRI", season = "Fall")
  bob_md_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[4]][[3]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "TRI", season = "Winter")
  bob_md_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[4]][[4]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "TRI", season = "Spring")
  #'  Bobcat-white-tailed deer
  bob_wtd_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[4]][[5]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Summer")
  #'  Low risk >50 bobcat; High risk <50 bobcat
  bob_wtd_tri_fall_out1 <- predprey_table(pred_prey_tri_overlap[[4]][[6]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Fall")
  bob_wtd_tri_fall_out4 <- predprey_table(pred_prey_tri_overlap[[4]][[7]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Fall")
  bob_wtd_tri_fall_out <- rbind(bob_wtd_tri_fall_out4[1,], bob_wtd_tri_fall_out1[2,])
  #'  Low risk >50 bobcat; High risk <50 bobcat
  bob_wtd_tri_wtr_out1 <- predprey_table(pred_prey_tri_overlap[[4]][[8]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Winter")
  bob_wtd_tri_wtr_out4 <- predprey_table(pred_prey_tri_overlap[[4]][[9]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Winter")
  bob_wtd_tri_wtr_out <- rbind(bob_wtd_tri_wtr_out4[1,], bob_wtd_tri_wtr_out1[2,])
  bob_wtd_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[4]][[10]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Spring")
  #'  Bobcat-prey results table
  bob_prey_tri_out <- rbind(bob_md_tri_smr_out, bob_md_tri_fall_out, bob_md_tri_wtr_out, bob_md_tri_sprg_out,
                            bob_wtd_tri_smr_out, bob_wtd_tri_fall_out, bob_wtd_tri_wtr_out, bob_wtd_tri_sprg_out)
  
  ####  Coyote - Prey TRI Output  ####
  #'  Coyote-mule deer
  coy_md_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[5]][[1]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "TRI", season = "Summer")
  coy_md_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[5]][[2]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "TRI", season = "Fall")
  coy_md_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[5]][[3]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "TRI", season = "Winter")
  coy_md_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[5]][[4]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "TRI", season = "Spring")
  #'  Coyote-white-tailed deer
  coy_wtd_tri_smr_out <- predprey_table(pred_prey_tri_overlap[[5]][[5]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Summer")
  coy_wtd_tri_fall_out <- predprey_table(pred_prey_tri_overlap[[5]][[6]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Fall")
  coy_wtd_tri_wtr_out <- predprey_table(pred_prey_tri_overlap[[5]][[7]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Winter")
  coy_wtd_tri_sprg_out <- predprey_table(pred_prey_tri_overlap[[5]][[8]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "TRI", season = "Spring")
  #'  Coyote-prey results table
  coy_prey_tri_out <- rbind(coy_md_tri_smr_out, coy_md_tri_fall_out, coy_md_tri_wtr_out, coy_md_tri_sprg_out,
                        coy_wtd_tri_smr_out, coy_wtd_tri_fall_out, coy_wtd_tri_wtr_out, coy_wtd_tri_sprg_out)
  
  pred_prey_tri_out_tbl <- rbind(coug_prey_tri_out, wolf_prey_tri_out, bear_prey_tri_out, bob_prey_tri_out, coy_prey_tri_out) %>%
    mutate(ndet_predator = as.numeric(ndet_predator),
           ndet_prey = as.numeric(ndet_prey))
  # write.csv(pred_prey_tri_out_tbl, file = paste0("./Outputs/Temporal Overlap/pred-prey_TRI_overlap_tbl_", Sys.Date(), ".csv"))

  #'  Drop results that are based on small sample sizes (<20 observations for one spp; Niedballa et al. 2019)
  pred_prey_tri_skinny <- pred_prey_tri_out_tbl %>%
    filter(ndet_predator >= 20) %>%  
    filter(ndet_prey >= 20) 
  # write.csv(pred_prey_tri_skinny, file = paste0("./Outputs/Temporal Overlap/pred-prey_TRI_overlap_tbl_skinny_", Sys.Date(), ".csv"))
  #'  Drop results that no longer had something for comparison in actual csv file - too annoying to do it in R
  
  
  ####  Cougar - Prey % Forest Output  ####
  #'  Cougar-mule deer
  coug_md_for_smr_out <- predprey_table(pred_prey_for_overlap[[1]][[1]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Summer")
  coug_md_for_fall_out <- predprey_table(pred_prey_for_overlap[[1]][[2]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Fall")
  coug_md_for_wtr_out <- predprey_table(pred_prey_for_overlap[[1]][[3]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Winter")
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_md_for_sprg_out1 <- predprey_table(pred_prey_for_overlap[[1]][[4]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Spring")
  coug_md_for_sprg_out4 <- predprey_table(pred_prey_for_overlap[[1]][[5]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Spring")
  coug_md_for_sprg_out <- rbind(coug_md_for_sprg_out4[1,], coug_md_for_sprg_out1[2,])
  #'  Cougar-elk
  #' Low risk <50 cougars; High risk >50 cougars
  coug_elk_for_smr_out1 <- predprey_table(pred_prey_for_overlap[[1]][[6]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Perc_Forest", season = "Summer")
  coug_elk_for_smr_out4 <- predprey_table(pred_prey_for_overlap[[1]][[7]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Perc_Forest", season = "Summer")
  coug_elk_for_smr_out <- rbind(coug_elk_for_smr_out1[1,], coug_elk_for_smr_out4[2,])
  coug_elk_for_fall_out <- predprey_table(pred_prey_for_overlap[[1]][[8]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Perc_Forest", season = "Fall")
  coug_elk_for_wtr_out <- predprey_table(pred_prey_for_overlap[[1]][[9]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Perc_Forest", season = "Winter")
  coug_elk_for_sprg_out <- predprey_table(pred_prey_for_overlap[[1]][[10]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Perc_Forest", season = "Spring")
  #'  Cougar-moose
  coug_moose_for_smr_out <- predprey_table(pred_prey_for_overlap[[1]][[11]], spp1 = "Cougar", spp2 = "Moose", risk_type = "Perc_Forest", season = "Summer")
  coug_moose_for_fall_out <- predprey_table(pred_prey_for_overlap[[1]][[12]], spp1 = "Cougar", spp2 = "Moose", risk_type = "Perc_Forest", season = "Fall")
  coug_moose_for_wtr_out <- predprey_table(pred_prey_for_overlap[[1]][[13]], spp1 = "Cougar", spp2 = "Moose", risk_type = "Perc_Forest", season = "Winter")
  coug_moose_for_sprg_out <- predprey_table(pred_prey_for_overlap[[1]][[14]], spp1 = "Cougar", spp2 = "Moose", risk_type = "Perc_Forest", season = "Spring")
  #'  Cougar-white-tailed deer
  coug_wtd_for_smr_out <- predprey_table(pred_prey_for_overlap[[1]][[15]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Summer")
  #'  Low risk <50 cougars; High risk >50 cougars
  coug_wtd_for_fall_out1 <- predprey_table(pred_prey_for_overlap[[1]][[16]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Fall")
  coug_wtd_for_fall_out4 <- predprey_table(pred_prey_for_overlap[[1]][[17]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Fall")
  coug_wtd_for_fall_out <- rbind(coug_wtd_for_fall_out1[1,], coug_wtd_for_fall_out4[2,])
  coug_wtd_for_wtr_out <- predprey_table(pred_prey_for_overlap[[1]][[18]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Winter")
  coug_wtd_for_sprg_out <- predprey_table(pred_prey_for_overlap[[1]][[19]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Spring")
    #'  Cougar-prey results table
  coug_prey_for_out <- rbind(coug_md_for_smr_out, coug_md_for_fall_out, coug_md_for_wtr_out, coug_md_for_sprg_out,
                             coug_elk_for_smr_out, coug_elk_for_fall_out, coug_elk_for_wtr_out, coug_elk_for_sprg_out,
                             coug_moose_for_smr_out, coug_moose_for_fall_out, coug_moose_for_wtr_out, coug_moose_for_sprg_out,
                             coug_wtd_for_smr_out, coug_wtd_for_fall_out, coug_wtd_for_wtr_out, coug_wtd_for_sprg_out)
  
  ####  Wolf - Prey % Forest Output  ####
  #'  Wolf-mule deer
  wolf_md_for_smr_out <- predprey_table(pred_prey_for_overlap[[2]][[1]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Summer")
  wolf_md_for_fall_out <- predprey_table(pred_prey_for_overlap[[2]][[2]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Fall")
  wolf_md_for_sprg_out <- predprey_table(pred_prey_for_overlap[[2]][[3]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Spring")
  #'  Wolf-elk
  wolf_elk_for_smr_out <- predprey_table(pred_prey_for_overlap[[2]][[4]], spp1 = "Wolf", spp2 = "Elk", risk_type = "Perc_Forest", season = "Summer")
  #'  Wolf-moose
  wolf_moose_for_smr_out <- predprey_table(pred_prey_for_overlap[[2]][[5]], spp1 = "Wolf", spp2 = "Moose", risk_type = "Perc_Forest", season = "Summer")
  wolf_moose_for_fall_out <- predprey_table(pred_prey_for_overlap[[2]][[6]], spp1 = "Wolf", spp2 = "Moose", risk_type = "Perc_Forest", season = "Fall")
  wolf_moose_for_wtr_out <- predprey_table(pred_prey_for_overlap[[2]][[7]], spp1 = "Wolf", spp2 = "Moose", risk_type = "Perc_Forest", season = "Winter")
  #'  Wolf-white-tailed deer
  wolf_wtd_for_smr_out <- predprey_table(pred_prey_for_overlap[[2]][[8]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Summer")
  wolf_wtd_for_fall_out <- predprey_table(pred_prey_for_overlap[[2]][[9]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Fall")
  wolf_wtd_for_wtr_out <- predprey_table(pred_prey_for_overlap[[2]][[10]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Winter")
  #'  Wolf-prey results table
  wolf_prey_for_out <- rbind(wolf_md_for_smr_out, wolf_md_for_fall_out, wolf_md_for_sprg_out, wolf_elk_for_smr_out,
                             wolf_moose_for_smr_out, wolf_moose_for_fall_out, wolf_moose_for_wtr_out,
                             wolf_wtd_for_smr_out, wolf_wtd_for_fall_out, wolf_wtd_for_wtr_out)
  
  ####  Black bear - Prey % Forest Output  ####
  #'  Black bear-mule deer
  bear_md_for_smr_out <- predprey_table(pred_prey_for_overlap[[3]][[1]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Summer")
  #'  Low risk >50 black bears; High risk <50 black bears
  bear_md_for_fall_out1 <- predprey_table(pred_prey_for_overlap[[3]][[2]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Fall")
  bear_md_for_fall_out4 <- predprey_table(pred_prey_for_overlap[[3]][[3]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Fall")
  bear_md_for_fall_out <- rbind(bear_md_for_fall_out4[1,], bear_md_for_fall_out1[2,])
  bear_md_for_sprg_out <- predprey_table(pred_prey_for_overlap[[3]][[4]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Spring")
  #'  Black bear-elk
  bear_elk_for_smr_out <- predprey_table(pred_prey_for_overlap[[3]][[5]], spp1 = "Black bear", spp2 = "Elk", risk_type = "Perc_Forest", season = "Summer")
  bear_elk_for_fall_out <- predprey_table(pred_prey_for_overlap[[3]][[6]], spp1 = "Black bear", spp2 = "Elk", risk_type = "Perc_Forest", season = "Fall")
  bear_elk_for_sprg_out <- predprey_table(pred_prey_for_overlap[[3]][[7]], spp1 = "Black bear", spp2 = "Elk", risk_type = "Perc_Forest", season = "Spring")
  #'  Black bear-moose
  bear_moose_for_smr_out <- predprey_table(pred_prey_for_overlap[[3]][[8]], spp1 = "Black bear", spp2 = "Moose", risk_type = "Perc_Forest", season = "Summer")
  bear_moose_for_fall_out <- predprey_table(pred_prey_for_overlap[[3]][[9]], spp1 = "Black bear", spp2 = "Moose", risk_type = "Perc_Forest", season = "Fall")
  bear_moose_for_sprg_out <- predprey_table(pred_prey_for_overlap[[3]][[10]], spp1 = "Black bear", spp2 = "Moose", risk_type = "Perc_Forest", season = "Spring")
  #'  Black bear-white-tailed deer
  bear_wtd_for_smr_out <- predprey_table(pred_prey_for_overlap[[3]][[11]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Summer")
  bear_wtd_for_fall_out <- predprey_table(pred_prey_for_overlap[[3]][[12]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Fall")
  #'  Low risk <50 black bears; High risk >50 black bears
  bear_wtd_for_sprg_out1 <- predprey_table(pred_prey_for_overlap[[3]][[13]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Spring")
  bear_wtd_for_sprg_out4 <- predprey_table(pred_prey_for_overlap[[3]][[14]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Spring")
  bear_wtd_for_sprg_out <- rbind(bear_wtd_for_sprg_out1[1,], bear_wtd_for_sprg_out4[2,])
  #'  Black bear-prey results table
  bear_prey_for_out <- rbind(bear_md_for_smr_out, bear_md_for_fall_out, bear_md_for_sprg_out,
                             bear_elk_for_smr_out, bear_elk_for_fall_out, bear_elk_for_sprg_out,
                             bear_moose_for_smr_out, bear_moose_for_fall_out, bear_moose_for_sprg_out,
                             bear_wtd_for_smr_out, bear_wtd_for_fall_out, bear_wtd_for_sprg_out)
  
  ####  Bobcat - Prey % Forest Output  ####
  #'  Bobcat-mule deer
  bob_md_for_smr_out <- predprey_table(pred_prey_for_overlap[[4]][[1]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Summer")
  bob_md_for_fall_out <- predprey_table(pred_prey_for_overlap[[4]][[2]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Fall")
  bob_md_for_sprg_out <- predprey_table(pred_prey_for_overlap[[4]][[3]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Spring")
  #'  Bobcat-white-tailed deer
  bob_wtd_for_smr_out <- predprey_table(pred_prey_for_overlap[[4]][[4]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Summer")
  #'  Low risk <50 bobcat; High risk >50 bobcat
  bob_wtd_for_fall_out1 <- predprey_table(pred_prey_for_overlap[[4]][[5]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Fall")
  bob_wtd_for_fall_out4 <- predprey_table(pred_prey_for_overlap[[4]][[6]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Fall")
  bob_wtd_for_fall_out <- rbind(bob_wtd_for_fall_out1[1,], bob_wtd_for_fall_out4[2,])
  #'  Low risk <50 bobcat; High risk >50 bobcat
  bob_wtd_for_wtr_out1 <- predprey_table(pred_prey_for_overlap[[4]][[7]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Winter")
  bob_wtd_for_wtr_out4 <- predprey_table(pred_prey_for_overlap[[4]][[8]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Winter")
  bob_wtd_for_wtr_out <- rbind(bob_wtd_for_wtr_out1[1,], bob_wtd_for_wtr_out4[2,])
  bob_wtd_for_sprg_out <- predprey_table(pred_prey_for_overlap[[4]][[9]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Spring")
  #'  Bobcat-prey results table
  bob_prey_for_out <- rbind(bob_md_for_smr_out, bob_md_for_fall_out, bob_md_for_sprg_out,
                            bob_wtd_for_smr_out, bob_wtd_for_fall_out, bob_wtd_for_wtr_out, bob_wtd_for_sprg_out)
  
  ####  Coyote - Prey % Forest Output  ####
  #'  Coyote-mule deer
  coy_md_for_smr_out <- predprey_table(pred_prey_for_overlap[[5]][[1]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Summer")
  coy_md_for_fall_out <- predprey_table(pred_prey_for_overlap[[5]][[2]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Fall")
  coy_md_for_wtr_out <- predprey_table(pred_prey_for_overlap[[5]][[3]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Winter")
  coy_md_for_sprg_out <- predprey_table(pred_prey_for_overlap[[5]][[4]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "Perc_Forest", season = "Spring")
  #'  Coyote-white-tailed deer
  coy_wtd_for_smr_out <- predprey_table(pred_prey_for_overlap[[5]][[5]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Summer")
  coy_wtd_for_fall_out <- predprey_table(pred_prey_for_overlap[[5]][[6]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Fall")
  coy_wtd_for_wtr_out <- predprey_table(pred_prey_for_overlap[[5]][[7]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Winter")
  coy_wtd_for_sprg_out <- predprey_table(pred_prey_for_overlap[[5]][[8]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "Perc_Forest", season = "Spring")
  #'  Coyote-prey results table
  coy_prey_for_out <- rbind(coy_md_for_smr_out, coy_md_for_fall_out, coy_md_for_wtr_out, coy_md_for_sprg_out,
                            coy_wtd_for_smr_out, coy_wtd_for_fall_out, coy_wtd_for_wtr_out, coy_wtd_for_sprg_out)
  
  pred_prey_for_out_tbl <- rbind(coug_prey_for_out, wolf_prey_for_out, bear_prey_for_out, bob_prey_for_out, coy_prey_for_out) %>%
    mutate(ndet_predator = as.numeric(ndet_predator),
           ndet_prey = as.numeric(ndet_prey))
  # write.csv(pred_prey_for_out_tbl, file = paste0("./Outputs/Temporal Overlap/pred-prey_PercForest_overlap_tbl_", Sys.Date(), ".csv"))
  
  
  #'  Drop results that are based on small sample sizes (<20 observations for one spp; Niedballa et al. 2019)
  pred_prey_for_skinny <- pred_prey_for_out_tbl %>%
    filter(ndet_predator >= 20) %>%  
    filter(ndet_prey >= 20) 
  # write.csv(pred_prey_for_skinny, file = paste0("./Outputs/Temporal Overlap/pred-prey_PercForest_overlap_tbl_skinny_", Sys.Date(), ".csv"))
  #'  Drop results that no longer had something for comparison in actual csv file - too annoying to do it in R
  
  #'  ------------------------------------
  ####  Species-specific results tables  ####
  #'  ------------------------------------
  #'  Create tables from prey-prey overlap estimates in low vs high risk areas
  prey_table <- function(overlap_out, spp2, season, risk_type) {
    Dhat <- round(overlap_out$dhats_spp2.lowhigh, 2)
    l95 <- round(overlap_out$CI[2,1], 2)
    u95 <- round(overlap_out$CI[2,2], 2)
    Species <- spp2
    Season <- season
    Risk <- risk_type
    ndet_low <- overlap_out$ndet_spp2.lowrisk
    ndet_high <- overlap_out$ndet_spp2.highrisk
    df <- as.data.frame(cbind(Species, Season, Risk, Dhat, l95, u95, ndet_low, ndet_high))
    rownames(df) <- NULL
    # names(df)[names(df) == "Risk"] <- paste0(risk_type, "_type")
    df <- mutate(df, Dhat = as.numeric(Dhat),
                 l95 = as.numeric(l95),
                 u95 = as.numeric(u95))
    return(df)
  }
  #'  Coefficient of overlap for each species at low vs high risk cameras
  ####  Mule deer Ouptut  ####
  #'  Summer
  md_smr_tri_out <- prey_table(prey_overlap[[1]][[1]], spp2 = "Mule deer", season = "Summer", risk_type = "TRI")
  md_smr_for_out <- prey_table(prey_overlap[[1]][[2]], spp2 = "Mule deer", season = "Summer", risk_type = "Percent Forest")
  md_smr_coug_out <- prey_table(prey_overlap[[1]][[3]], spp2 = "Mule deer", season = "Summer", risk_type = "Cougar detected")
  md_smr_wolf_out <- prey_table(prey_overlap[[1]][[4]], spp2 = "Mule deer", season = "Summer", risk_type = "Wolf detected")
  md_smr_bear_out <- prey_table(prey_overlap[[1]][[5]], spp2 = "Mule deer", season = "Summer", risk_type = "Black bear detected")
  md_smr_bob_out <- prey_table(prey_overlap[[1]][[6]], spp2 = "Mule deer", season = "Summer", risk_type = "Bobcat detected")
  md_smr_coy_out <- prey_table(prey_overlap[[1]][[7]], spp2 = "Mule deer", season = "Summer", risk_type = "Coyote detected")
  #'  Fall
  md_fall_tri_out <- prey_table(prey_overlap[[1]][[8]], spp2 = "Mule deer", season = "Fall", risk_type = "TRI")
  md_fall_for_out <- prey_table(prey_overlap[[1]][[9]], spp2 = "Mule deer", season = "Fall", risk_type = "Percent Forest")
  md_fall_coug_out <- prey_table(prey_overlap[[1]][[10]], spp2 = "Mule deer", season = "Fall", risk_type = "Cougar detected")
  md_fall_wolf_out <- prey_table(prey_overlap[[1]][[11]], spp2 = "Mule deer", season = "Fall", risk_type = "Wolf detected")
  md_fall_bear_out <- prey_table(prey_overlap[[1]][[12]], spp2 = "Mule deer", season = "Fall", risk_type = "Black bear detected")
  md_fall_bob_out <- prey_table(prey_overlap[[1]][[13]], spp2 = "Mule deer", season = "Fall", risk_type = "Bobcat detected")
  md_fall_coy_out <- prey_table(prey_overlap[[1]][[14]], spp2 = "Mule deer", season = "Fall", risk_type = "Coyote detected")
  #'  Winter
  md_wtr_tri_out <- prey_table(prey_overlap[[1]][[15]], spp2 = "Mule deer", season = "Winter", risk_type = "TRI")
  md_wtr_for_out <- prey_table(prey_overlap[[1]][[16]], spp2 = "Mule deer", season = "Winter", risk_type = "Percent Forest")
  md_wtr_coug_out <- prey_table(prey_overlap[[1]][[17]], spp2 = "Mule deer", season = "Winter", risk_type = "Cougar detected")
  md_wtr_wolf_out <- prey_table(prey_overlap[[1]][[18]], spp2 = "Mule deer", season = "Winter", risk_type = "Wolf detected")
  md_wtr_bob_out <- prey_table(prey_overlap[[1]][[19]], spp2 = "Mule deer", season = "Winter", risk_type = "Bobcat detected")
  md_wtr_coy_out <- prey_table(prey_overlap[[1]][[20]], spp2 = "Mule deer", season = "Winter", risk_type = "Coyote detected")
  #'  Spring
  md_sprg_tri_out <- prey_table(prey_overlap[[1]][[21]], spp2 = "Mule deer", season = "Spring", risk_type = "TRI")
  md_sprg_for_out <- prey_table(prey_overlap[[1]][[22]], spp2 = "Mule deer", season = "Spring", risk_type = "Percent Forest")
  md_sprg_coug_out <- prey_table(prey_overlap[[1]][[23]], spp2 = "Mule deer", season = "Spring", risk_type = "Cougar detected")
  md_sprg_wolf_out <- prey_table(prey_overlap[[1]][[24]], spp2 = "Mule deer", season = "Spring", risk_type = "Wolf detected")
  md_sprg_bear_out <- prey_table(prey_overlap[[1]][[25]], spp2 = "Mule deer", season = "Spring", risk_type = "Black bear detected")
  md_sprg_bob_out <- prey_table(prey_overlap[[1]][[26]], spp2 = "Mule deer", season = "Spring", risk_type = "Bobcat detected")
  md_sprg_coy_out <- prey_table(prey_overlap[[1]][[27]], spp2 = "Mule deer", season = "Spring", risk_type = "Coyote detected")
  #'  Save mule deer overlap results
  md_overlap <- rbind(md_smr_tri_out, md_smr_for_out, md_smr_coug_out, md_smr_wolf_out, md_smr_bear_out, md_smr_bob_out, md_smr_coy_out,
                      md_fall_tri_out, md_fall_for_out, md_fall_coug_out, md_fall_wolf_out, md_fall_bear_out, md_fall_bob_out, md_fall_coy_out,
                      md_wtr_tri_out, md_wtr_for_out, md_wtr_coug_out, md_wtr_wolf_out, md_wtr_bob_out, md_wtr_coy_out,
                      md_sprg_tri_out, md_sprg_for_out, md_sprg_coug_out, md_sprg_wolf_out, md_sprg_bear_out, md_sprg_bob_out, md_sprg_coy_out)
  
  ####  Elk Ouptut  ####
  #'  Summer
  elk_smr_tri_out <- prey_table(prey_overlap[[2]][[1]], spp2 = "Elk", season = "Summer", risk_type = "TRI")
  elk_smr_for_out <- prey_table(prey_overlap[[2]][[2]], spp2 = "Elk", season = "Summer", risk_type = "Percent Forest")
  elk_smr_coug_out <- prey_table(prey_overlap[[2]][[3]], spp2 = "Elk", season = "Summer", risk_type = "Cougar detected")
  elk_smr_wolf_out <- prey_table(prey_overlap[[2]][[4]], spp2 = "Elk", season = "Summer", risk_type = "Wolf detected")
  elk_smr_bear_out <- prey_table(prey_overlap[[2]][[5]], spp2 = "Elk", season = "Summer", risk_type = "Black bear detected")
  #'  Fall
  elk_fall_tri_out <- prey_table(prey_overlap[[2]][[6]], spp2 = "Elk", season = "Fall", risk_type = "TRI")
  elk_fall_for_out <- prey_table(prey_overlap[[2]][[7]], spp2 = "Elk", season = "Fall", risk_type = "Percent Forest")
  elk_fall_coug_out <- prey_table(prey_overlap[[2]][[8]], spp2 = "Elk", season = "Fall", risk_type = "Cougar detected")
  elk_fall_bear_out <- prey_table(prey_overlap[[2]][[9]], spp2 = "Elk", season = "Fall", risk_type = "Black bear detected")
  #'  Winter
  elk_wtr_tri_out <- prey_table(prey_overlap[[2]][[10]], spp2 = "Elk", season = "Winter", risk_type = "TRI")
  elk_wtr_for_out <- prey_table(prey_overlap[[2]][[11]], spp2 = "Elk", season = "Winter", risk_type = "Percent Forest")
  elk_wtr_coug_out <- prey_table(prey_overlap[[2]][[12]], spp2 = "Elk", season = "Winter", risk_type = "Cougar detected")
  #'  Spring
  elk_sprg_tri_out <- prey_table(prey_overlap[[2]][[13]], spp2 = "Elk", season = "Spring", risk_type = "TRI")
  elk_sprg_for_out <- prey_table(prey_overlap[[2]][[14]], spp2 = "Elk", season = "Spring", risk_type = "Percent Forest")
  elk_sprg_coug_out <- prey_table(prey_overlap[[2]][[15]], spp2 = "Elk", season = "Spring", risk_type = "Cougar detected")
  elk_sprg_bear_out <- prey_table(prey_overlap[[2]][[16]], spp2 = "Elk", season = "Spring", risk_type = "Black bear detected")
  #'  Save elk overlap results
  elk_overlap <- rbind(elk_smr_tri_out, elk_smr_for_out, elk_smr_coug_out, elk_smr_wolf_out, elk_smr_bear_out,
                       elk_fall_tri_out, elk_fall_for_out, elk_fall_coug_out, elk_fall_bear_out,
                       elk_wtr_tri_out, elk_wtr_for_out, elk_wtr_coug_out, 
                       elk_sprg_tri_out, elk_sprg_for_out, elk_sprg_coug_out, elk_sprg_bear_out)
  
  ####  Moose Ouptut  ####
  #'  Summer
  moose_smr_tri_out <- prey_table(prey_overlap[[3]][[1]], spp2 = "Moose", season = "Summer", risk_type = "TRI")
  moose_smr_for_out <- prey_table(prey_overlap[[3]][[2]], spp2 = "Moose", season = "Summer", risk_type = "Percent Forest")
  moose_smr_coug_out <- prey_table(prey_overlap[[3]][[3]], spp2 = "Moose", season = "Summer", risk_type = "Cougar detected")
  moose_smr_wolf_out <- prey_table(prey_overlap[[3]][[4]], spp2 = "Moose", season = "Summer", risk_type = "Wolf detected")
  moose_smr_bear_out <- prey_table(prey_overlap[[3]][[5]], spp2 = "Moose", season = "Summer", risk_type = "Black bear detected")
  #'  Fall
  moose_fall_tri_out <- prey_table(prey_overlap[[3]][[6]], spp2 = "Moose", season = "Fall", risk_type = "TRI")
  moose_fall_for_out <- prey_table(prey_overlap[[3]][[7]], spp2 = "Moose", season = "Fall", risk_type = "Percent Forest")
  moose_fall_coug_out <- prey_table(prey_overlap[[3]][[8]], spp2 = "Moose", season = "Fall", risk_type = "Cougar detected")
  moose_fall_wolf_out <- prey_table(prey_overlap[[3]][[9]], spp2 = "Moose", season = "Fall", risk_type = "Wolf detected")
  moose_fall_bear_out <- prey_table(prey_overlap[[3]][[10]], spp2 = "Moose", season = "Fall", risk_type = "Black bear detected")
  #'  Winter
  moose_wtr_tri_out <- prey_table(prey_overlap[[3]][[11]], spp2 = "Moose", season = "Winter", risk_type = "TRI")
  moose_wtr_for_out <- prey_table(prey_overlap[[3]][[12]], spp2 = "Moose", season = "Winter", risk_type = "Percent Forest")
  moose_wtr_coug_out <- prey_table(prey_overlap[[3]][[13]], spp2 = "Moose", season = "Winter", risk_type = "Cougar detected")
  moose_wtr_wolf_out <- prey_table(prey_overlap[[3]][[14]], spp2 = "Moose", season = "Winter", risk_type = "Wolf detected")
  #'  Spring
  moose_sprg_tri_out <- prey_table(prey_overlap[[3]][[15]], spp2 = "Moose", season = "Spring", risk_type = "TRI")
  moose_sprg_for_out <- prey_table(prey_overlap[[3]][[16]], spp2 = "Moose", season = "Spring", risk_type = "Percent Forest")
  moose_sprg_coug_out <- prey_table(prey_overlap[[3]][[17]], spp2 = "Moose", season = "Spring", risk_type = "Cougar detected")
  moose_sprg_wolf_out <- prey_table(prey_overlap[[3]][[18]], spp2 = "Moose", season = "Spring", risk_type = "Wolf detected")
  moose_sprg_bear_out <- prey_table(prey_overlap[[3]][[19]], spp2 = "Moose", season = "Spring", risk_type = "Black bear detected")
  #'  Save elk overlap results
  moose_overlap <- rbind(moose_smr_tri_out, moose_smr_for_out, moose_smr_coug_out, moose_smr_wolf_out, moose_smr_bear_out,
                         moose_fall_tri_out, moose_fall_for_out, moose_fall_coug_out, moose_fall_wolf_out, moose_fall_bear_out,
                         moose_wtr_tri_out, moose_wtr_for_out, moose_wtr_coug_out,  moose_wtr_wolf_out, 
                         moose_sprg_tri_out, moose_sprg_for_out, moose_sprg_coug_out, moose_sprg_wolf_out, moose_sprg_bear_out)
  
  ####  White-tailed deer Ouptut  ####
  #'  Summer
  wtd_smr_tri_out <- prey_table(prey_overlap[[4]][[1]], spp2 = "White-tailed deer", season = "Summer", risk_type = "TRI")
  wtd_smr_for_out <- prey_table(prey_overlap[[4]][[2]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Percent Forest")
  wtd_smr_coug_out <- prey_table(prey_overlap[[4]][[3]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Cougar detected")
  wtd_smr_wolf_out <- prey_table(prey_overlap[[4]][[4]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Wolf detected")
  wtd_smr_bear_out <- prey_table(prey_overlap[[4]][[5]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Black bear detected")
  wtd_smr_bob_out <- prey_table(prey_overlap[[4]][[6]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Bobcat detected")
  wtd_smr_coy_out <- prey_table(prey_overlap[[4]][[7]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Coyote detected")
  #'  Fall
  wtd_fall_tri_out <- prey_table(prey_overlap[[4]][[8]], spp2 = "White-tailed deer", season = "Fall", risk_type = "TRI")
  wtd_fall_for_out <- prey_table(prey_overlap[[4]][[9]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Percent Forest")
  wtd_fall_coug_out <- prey_table(prey_overlap[[4]][[10]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Cougar detected")
  wtd_fall_wolf_out <- prey_table(prey_overlap[[4]][[11]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Wolf detected")
  wtd_fall_bear_out <- prey_table(prey_overlap[[4]][[12]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Black bear detected")
  wtd_fall_bob_out <- prey_table(prey_overlap[[4]][[13]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Bobcat detected")
  wtd_fall_coy_out <- prey_table(prey_overlap[[4]][[14]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Coyote detected")
  #'  Winter
  wtd_wtr_tri_out <- prey_table(prey_overlap[[4]][[15]], spp2 = "White-tailed deer", season = "Winter", risk_type = "TRI")
  wtd_wtr_for_out <- prey_table(prey_overlap[[4]][[16]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Percent Forest")
  wtd_wtr_coug_out <- prey_table(prey_overlap[[4]][[17]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Cougar detected")
  wtd_wtr_wolf_out <- prey_table(prey_overlap[[4]][[18]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Wolf detected")
  wtd_wtr_bob_out <- prey_table(prey_overlap[[4]][[19]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Bobcat detected")
  wtd_wtr_coy_out <- prey_table(prey_overlap[[4]][[20]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Coyote detected")
  #'  Spring
  wtd_sprg_tri_out <- prey_table(prey_overlap[[4]][[21]], spp2 = "White-tailed deer", season = "Spring", risk_type = "TRI")
  wtd_sprg_for_out <- prey_table(prey_overlap[[4]][[22]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Percent Forest")
  wtd_sprg_coug_out <- prey_table(prey_overlap[[4]][[23]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Cougar detected")
  wtd_sprg_wolf_out <- prey_table(prey_overlap[[4]][[24]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Wolf detected")
  wtd_sprg_bear_out <- prey_table(prey_overlap[[4]][[25]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Black bear detected")
  wtd_sprg_bob_out <- prey_table(prey_overlap[[4]][[26]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Bobcat detected")
  wtd_sprg_coy_out <- prey_table(prey_overlap[[4]][[27]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Coyote detected")
  #'  Save mule deer overlap results
  wtd_overlap <- rbind(wtd_smr_tri_out, wtd_smr_for_out, wtd_smr_coug_out, wtd_smr_wolf_out, wtd_smr_bear_out, wtd_smr_bob_out, wtd_smr_coy_out,
                       wtd_fall_tri_out, wtd_fall_for_out, wtd_fall_coug_out, wtd_fall_wolf_out, wtd_fall_bear_out, wtd_fall_bob_out, wtd_fall_coy_out,
                       wtd_wtr_tri_out, wtd_wtr_for_out, wtd_wtr_coug_out, wtd_wtr_wolf_out, wtd_wtr_bob_out, wtd_wtr_coy_out,
                       wtd_sprg_tri_out, wtd_sprg_for_out, wtd_sprg_coug_out, wtd_sprg_wolf_out, wtd_sprg_bear_out, wtd_sprg_bob_out, wtd_sprg_coy_out)
  
  prey_out_tbl <- rbind(md_overlap, elk_overlap, moose_overlap, wtd_overlap) %>%
    mutate(ndet_low = as.numeric(ndet_low),
           ndet_high = as.numeric(ndet_high))
  # write.csv(prey_out_tbl, file = paste0("./Outputs/Temporal Overlap/prey_overlap_tbl_", Sys.Date(), ".csv"))
  
  #'  Drop results that are based on small sample sizes (<20 observations for one curve; Niedballa et al. 2019)
  prey_tbl_skinny <- prey_out_tbl %>%
    filter(ndet_low >= 20) %>%  
    filter(ndet_high >= 20)
  # write.csv(prey_tbl_skinny, file = paste0("./Outputs/Temporal Overlap/prey_overlap_tbl_skinny_", Sys.Date(), ".csv"))
  
  
  ####  Black bear Output  ####
  #'  Summer
  bear_smr_tri_out <- prey_table(pred_overlap[[1]][[1]], spp2 = "Black bear", season = "Summer", risk_type = "TRI")
  bear_smr_for_out <- prey_table(pred_overlap[[1]][[2]], spp2 = "Black bear", season = "Summer", risk_type = "Percent Forest")
  #'  Fall
  bear_fall_tri_out <- prey_table(pred_overlap[[1]][[3]], spp2 = "Black bear", season = "Fall", risk_type = "TRI")
  bear_fall_for_out <- prey_table(pred_overlap[[1]][[4]], spp2 = "Black bear", season = "Fall", risk_type = "Percent Forest")
  #'  Spring
  bear_sprg_tri_out <- prey_table(pred_overlap[[1]][[5]], spp2 = "Black bear", season = "Spring", risk_type = "TRI")
  bear_sprg_for_out <- prey_table(pred_overlap[[1]][[6]], spp2 = "Black bear", season = "Spring", risk_type = "Percent Forest")
  
  bear_out_tbl <- rbind(bear_smr_tri_out, bear_smr_for_out, bear_fall_tri_out, bear_fall_for_out,
                        bear_sprg_tri_out, bear_sprg_for_out)
  
  ####  Bobcat Output  ####
  #'  Summer
  bob_smr_tri_out <- prey_table(pred_overlap[[2]][[1]], spp2 = "Bobcat", season = "Summer", risk_type = "TRI")
  bob_smr_for_out <- prey_table(pred_overlap[[2]][[2]], spp2 = "Bobcat", season = "Summer", risk_type = "Percent Forest")
  #'  Fall
  bob_fall_tri_out <- prey_table(pred_overlap[[2]][[3]], spp2 = "Bobcat", season = "Fall", risk_type = "TRI")
  bob_fall_for_out <- prey_table(pred_overlap[[2]][[4]], spp2 = "Bobcat", season = "Fall", risk_type = "Percent Forest")
  #'  Winter
  bob_wtr_tri_out <- prey_table(pred_overlap[[2]][[5]], spp2 = "Bobcat", season = "Winter", risk_type = "TRI")
  bob_wtr_for_out <- prey_table(pred_overlap[[2]][[6]], spp2 = "Bobcat", season = "Winter", risk_type = "Percent Forest")
  #'  Spring
  bob_sprg_tri_out <- prey_table(pred_overlap[[2]][[7]], spp2 = "Bobcat", season = "Spring", risk_type = "TRI")
  bob_sprg_for_out <- prey_table(pred_overlap[[2]][[8]], spp2 = "Bobcat", season = "Spring", risk_type = "Percent Forest")
  
  bob_out_tbl <- rbind(bob_smr_tri_out, bob_smr_for_out, bob_fall_tri_out, bob_fall_for_out,
                       bob_wtr_tri_out, bob_wtr_for_out, bob_sprg_tri_out, bob_sprg_for_out)
  
  ####  Cougar Output  ####
  #'  Summer
  coug_smr_tri_out <- prey_table(pred_overlap[[3]][[1]], spp2 = "Cougar", season = "Summer", risk_type = "TRI")
  coug_smr_for_out <- prey_table(pred_overlap[[3]][[2]], spp2 = "Cougar", season = "Summer", risk_type = "Percent Forest")
  #'  Fall
  coug_fall_tri_out <- prey_table(pred_overlap[[3]][[3]], spp2 = "Cougar", season = "Fall", risk_type = "TRI")
  coug_fall_for_out <- prey_table(pred_overlap[[3]][[4]], spp2 = "Cougar", season = "Fall", risk_type = "Percent Forest")
  #'  Winter
  coug_wtr_tri_out <- prey_table(pred_overlap[[3]][[5]], spp2 = "Cougar", season = "Winter", risk_type = "TRI")
  coug_wtr_for_out <- prey_table(pred_overlap[[3]][[6]], spp2 = "Cougar", season = "Winter", risk_type = "Percent Forest")
  #'  Spring
  coug_sprg_tri_out <- prey_table(pred_overlap[[3]][[7]], spp2 = "Cougar", season = "Spring", risk_type = "TRI")
  coug_sprg_for_out <- prey_table(pred_overlap[[3]][[8]], spp2 = "Cougar", season = "Spring", risk_type = "Percent Forest")
  
  coug_out_tbl <- rbind(coug_smr_tri_out, coug_smr_for_out, coug_fall_tri_out, coug_fall_for_out,
                        coug_wtr_tri_out, coug_wtr_for_out, coug_sprg_tri_out, coug_sprg_for_out)
  
  ####  Coyote Output  ####
  #'  Summer
  coy_smr_tri_out <- prey_table(pred_overlap[[4]][[1]], spp2 = "Coyote", season = "Summer", risk_type = "TRI")
  coy_smr_for_out <- prey_table(pred_overlap[[4]][[2]], spp2 = "Coyote", season = "Summer", risk_type = "Percent Forest")
  #'  Fall
  coy_fall_tri_out <- prey_table(pred_overlap[[4]][[3]], spp2 = "Coyote", season = "Fall", risk_type = "TRI")
  coy_fall_for_out <- prey_table(pred_overlap[[4]][[4]], spp2 = "Coyote", season = "Fall", risk_type = "Percent Forest")
  #'  Winter
  coy_wtr_tri_out <- prey_table(pred_overlap[[4]][[5]], spp2 = "Coyote", season = "Winter", risk_type = "TRI")
  coy_wtr_for_out <- prey_table(pred_overlap[[4]][[6]], spp2 = "Coyote", season = "Winter", risk_type = "Percent Forest")
  #'  Spring
  coy_sprg_tri_out <- prey_table(pred_overlap[[4]][[7]], spp2 = "Coyote", season = "Spring", risk_type = "TRI")
  coy_sprg_for_out <- prey_table(pred_overlap[[4]][[8]], spp2 = "Coyote", season = "Spring", risk_type = "Percent Forest")
  
  coy_out_tbl <- rbind(coy_smr_tri_out, coy_smr_for_out, coy_fall_tri_out, coy_fall_for_out,
                       coy_wtr_tri_out, coy_wtr_for_out, coy_sprg_tri_out, coy_sprg_for_out)
  
  ####  Wolf Output  ####
  #'  Summer
  wolf_smr_tri_out <- prey_table(pred_overlap[[5]][[1]], spp2 = "Wolf", season = "Summer", risk_type = "TRI")
  wolf_smr_for_out <- prey_table(pred_overlap[[5]][[2]], spp2 = "Wolf", season = "Summer", risk_type = "Percent Forest")
  #'  Fall
  wolf_fall_tri_out <- prey_table(pred_overlap[[5]][[3]], spp2 = "Wolf", season = "Fall", risk_type = "TRI")
  wolf_fall_for_out <- prey_table(pred_overlap[[5]][[4]], spp2 = "Wolf", season = "Fall", risk_type = "Percent Forest")
  #'  Winter
  wolf_wtr_tri_out <- prey_table(pred_overlap[[5]][[5]], spp2 = "Wolf", season = "Winter", risk_type = "TRI")
  wolf_wtr_for_out <- prey_table(pred_overlap[[5]][[6]], spp2 = "Wolf", season = "Winter", risk_type = "Percent Forest")
  #'  Spring
  wolf_sprg_tri_out <- prey_table(pred_overlap[[5]][[7]], spp2 = "Wolf", season = "Spring", risk_type = "TRI")
  wolf_sprg_for_out <- prey_table(pred_overlap[[5]][[8]], spp2 = "Wolf", season = "Spring", risk_type = "Percent Forest")
  
  wolf_out_tbl <- rbind(wolf_smr_tri_out, wolf_smr_for_out, wolf_fall_tri_out, wolf_fall_for_out,
                        wolf_wtr_tri_out, wolf_wtr_for_out, wolf_sprg_tri_out, wolf_sprg_for_out)
  
  pred_out_tbl <- rbind(bear_out_tbl, bob_out_tbl, coug_out_tbl, coy_out_tbl, wolf_out_tbl) %>%
    mutate(ndet_low = as.numeric(ndet_low),
           ndet_high = as.numeric(ndet_high))
  # write.csv(pred_out_tbl, file = paste0("./Outputs/Temporal Overlap/pred_overlap_tbl_", Sys.Date(), ".csv"))
  
  #'  Drop results that are based on small sample sizes (<20 observations for one curve; Niedballa et al. 2019)
  pred_tbl_skinny <- pred_out_tbl %>%
    filter(ndet_low >= 20) %>%  
    filter(ndet_high >= 20)
  # write.csv(pred_tbl_skinny, file = paste0("./Outputs/Temporal Overlap/pred_overlap_tbl_skinny_", Sys.Date(), ".csv"))
  
    
  #'  ---------------------------------
  ####  Coefficient Comparison Plots  ####
  #'  ---------------------------------
  #'  Compare coefficient of overlapping estimates for each species-pairing or
  #'  individual species at cameras with low versus high background risk
  
  #'  Load data in table format
  pred_prey_tri_overlap_tbl <- read.csv("./Outputs/Temporal Overlap/pred-prey_TRI_overlap_tbl_skinny_2022-10-15.csv") #pred-prey_TRI_overlap_tbl_2022-09-19
  pred_prey_for_overlap_tbl <- read.csv("./Outputs/Temporal Overlap/pred-prey_PercForest_overlap_tbl_skinny_2022-10-15.csv") #pred-prey_PercForest_overlap_tbl_2022-09-19
  prey_overlap_tbl <- read.csv("./Outputs/Temporal Overlap/prey_overlap_tbl_skinny_2022-10-15.csv") %>% #prey_overlap_tbl_2022-09-26
    mutate(Risk = ifelse(Risk == "Black bear detected", "Bear", Risk),
           Risk = ifelse(Risk == "Bobcat detected", "Bobcat", Risk),
           Risk = ifelse(Risk == "Coyote detected", "Coyote", Risk),
           Risk = ifelse(Risk == "Cougar detected", "Cougar", Risk),
           Risk = ifelse(Risk == "Wolf detected", "Wolf", Risk))
  pred_overlap_tbl <- read.csv("./Outputs/Temporal Overlap/pred_overlap_tbl_skinny_2022-10-15.csv") %>%
    mutate(Season = factor(Season, levels = c("Summer", "Fall", "Winter", "Spring")))
  
  #'  Set seasonal factor levels
  pred_prey_tri_overlap_tbl$season <- factor(pred_prey_tri_overlap_tbl$season, levels = c("Summer", "Fall", "Winter", "Spring"))
  pred_prey_tri_overlap_tbl$TRI_level <- factor(pred_prey_tri_overlap_tbl$TRI_level, levels = c("Low", "High"))
  pred_prey_for_overlap_tbl$season <- factor(pred_prey_for_overlap_tbl$season, levels = c("Summer", "Fall", "Winter", "Spring"))
  pred_prey_for_overlap_tbl$Perc_Forest_level <- factor(pred_prey_for_overlap_tbl$Perc_Forest_level, levels = c("Low", "High"))
  prey_overlap_tbl$Risk <- factor(prey_overlap_tbl$Risk, levels = c("TRI", "Percent Forest", "Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))
  pred_overlap_tbl$Risk <- factor(pred_overlap_tbl$Risk, levels = c("TRI", "Percent Forest"))
    
  ####  Effect of TRI on predator-prey temporal overlap  ####
  #'  By species
  pred_md_tri_coeff_plot <- ggplot(pred_prey_tri_overlap_tbl[pred_prey_tri_overlap_tbl$prey == "Mule Deer",], aes(x = predator, y = Dhat, group = TRI_level)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = TRI_level), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Terrain ruggedness")) + 
    ggtitle("Mule deer - predator temporal overlap") +
    xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~season, scales = "free", space = "free")
  pred_md_tri_coeff_plot
  
  pred_elk_tri_coeff_plot <- ggplot(pred_prey_tri_overlap_tbl[pred_prey_tri_overlap_tbl$prey == "Elk",], aes(x = predator, y = Dhat, group = TRI_level)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = TRI_level), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Terrain ruggedness")) + 
    ggtitle("Elk - predator temporal overlap") +
    xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~season, scales = "free", space = "free")
  pred_elk_tri_coeff_plot
  
  pred_moose_tri_coeff_plot <- ggplot(pred_prey_tri_overlap_tbl[pred_prey_tri_overlap_tbl$prey == "Moose",], aes(x = predator, y = Dhat, group = TRI_level)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = TRI_level), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Terrain ruggedness")) + 
    ggtitle("Moose - predator temporal overlap") +
    xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~season, scales = "free", space = "free")
  pred_moose_tri_coeff_plot
  
  pred_wtd_tri_coeff_plot <- ggplot(pred_prey_tri_overlap_tbl[pred_prey_tri_overlap_tbl$prey == "White-tailed Deer",], aes(x = predator, y = Dhat, group = TRI_level)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = TRI_level), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Terrain ruggedness")) + 
    ggtitle("White-tailed deer - predator temporal overlap") +
    xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~season, scales = "free", space = "free")
  pred_wtd_tri_coeff_plot
  
  #' Save 'em - predator-prey overlap coefficients
  ggsave(pred_md_tri_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_md_TRI_bckgrd_risk_skinny.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_elk_tri_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_elk_TRI_bckgrd_risk_skinny.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_moose_tri_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_moose_TRI_bckgrd_risk_skinny.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_wtd_tri_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_wtd_TRI_bckgrd_risk_skinny.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')

  
  #' #'  By season
  #' pred_prey_tri_smr_coeff_plot <- ggplot(pred_prey_tri_overlap_tbl[pred_prey_tri_overlap_tbl$season == "Summer",], aes(x = predator, y = Dhat, group = TRI_level)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = predator, shape = TRI_level), size = 2.75, position = position_dodge(width = 0.4)) +
  #'   scale_color_bright() + 
  #'   ylim(0,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none", shape = guide_legend(title = "Terrain ruggedness")) + 
  #'   ggtitle("Summer predator-prey temporal overlap") +
  #'   xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~prey, scales = "free", space = "free")
  #' pred_prey_tri_smr_coeff_plot
  #' 
  #' pred_prey_tri_fall_coeff_plot <- ggplot(pred_prey_tri_overlap_tbl[pred_prey_tri_overlap_tbl$season == "Fall",], aes(x = predator, y = Dhat, group = TRI_level)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = predator, shape = TRI_level), size = 2.75, position = position_dodge(width = 0.4)) +
  #'   scale_color_bright() + 
  #'   ylim(0,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none", shape = guide_legend(title = "Terrain ruggedness")) + 
  #'   ggtitle("Fall predator-prey temporal overlap") +
  #'   xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~prey, scales = "free", space = "free")
  #' pred_prey_tri_fall_coeff_plot
  #' 
  #' pred_prey_tri_wtr_coeff_plot <- ggplot(pred_prey_tri_overlap_tbl[pred_prey_tri_overlap_tbl$season == "Winter",], aes(x = predator, y = Dhat, group = TRI_level)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = predator, shape = TRI_level), size = 2.75, position = position_dodge(width = 0.4)) +
  #'   scale_color_bright() + 
  #'   ylim(0,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none", shape = guide_legend(title = "Terrain ruggedness")) + 
  #'   ggtitle("Winter predator-prey temporal overlap") +
  #'   xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~prey, scales = "free", space = "free")
  #' pred_prey_tri_wtr_coeff_plot
  #' 
  #' pred_prey_tri_sprg_coeff_plot <- ggplot(pred_prey_tri_overlap_tbl[pred_prey_tri_overlap_tbl$season == "Spring",], aes(x = predator, y = Dhat, group = TRI_level)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = predator, shape = TRI_level), size = 2.75, position = position_dodge(width = 0.4)) +
  #'   scale_color_bright() + 
  #'   ylim(0,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none", shape = guide_legend(title = "Terrain ruggedness")) + 
  #'   ggtitle("Spring predator-prey temporal overlap") +
  #'   xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~prey, scales = "free", space = "free")
  #' pred_prey_tri_sprg_coeff_plot
  #' 
  #' #' Save 'em - predator-prey overlap coefficients
  #' ggsave(pred_prey_tri_smr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_TRI_smr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(pred_prey_tri_fall_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_TRI_fall_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(pred_prey_tri_wtr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_TRI_wtr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(pred_prey_tri_sprg_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_TRI_sprg_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ####  Effect of % forest on predator-prey temporal overlap  ####
  #'  By species
  pred_md_for_coeff_plot <- ggplot(pred_prey_for_overlap_tbl[pred_prey_for_overlap_tbl$prey == "Mule Deer",], aes(x = predator, y = Dhat, group = Perc_Forest_level)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = Perc_Forest_level), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Percent forested habitat")) + 
    ggtitle("Mule deer - predator temporal overlap") +
    xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~season, scales = "free", space = "free")
  pred_md_for_coeff_plot
  
  pred_elk_for_coeff_plot <- ggplot(pred_prey_for_overlap_tbl[pred_prey_for_overlap_tbl$prey == "Elk",], aes(x = predator, y = Dhat, group = Perc_Forest_level)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = Perc_Forest_level), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Percent forested habitat")) + 
    ggtitle("Elk - predator temporal overlap") +
    xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~season, scales = "free", space = "free")
  pred_elk_for_coeff_plot
  
  pred_moose_for_coeff_plot <- ggplot(pred_prey_for_overlap_tbl[pred_prey_for_overlap_tbl$prey == "Moose",], aes(x = predator, y = Dhat, group = Perc_Forest_level)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = Perc_Forest_level), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Percent forested habitat")) + 
    ggtitle("Moose - predator temporal overlap") +
    xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~season, scales = "free", space = "free")
  pred_moose_for_coeff_plot
  
  pred_wtd_for_coeff_plot <- ggplot(pred_prey_for_overlap_tbl[pred_prey_for_overlap_tbl$prey == "White-tailed Deer",], aes(x = predator, y = Dhat, group = Perc_Forest_level)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = Perc_Forest_level), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Percent forested habitat")) + 
    ggtitle("White-tailed deer - predator temporal overlap") +
    xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~season, scales = "free", space = "free")
  pred_wtd_for_coeff_plot
  
  #' Save 'em - predator-prey overlap coefficients
  ggsave(pred_md_for_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_md_PercForest_bckgrd_risk_skinny.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_elk_for_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_elk_PercForest_bckgrd_risk_skinny.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_moose_for_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_moose_PercForest_bckgrd_risk_skinny.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_wtd_for_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_wtd_PercForest_bckgrd_risk_skinny.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  
  #' #'  By season
  #' pred_prey_for_smr_coeff_plot <- ggplot(pred_prey_for_overlap_tbl[pred_prey_for_overlap_tbl$season == "Summer",], aes(x = predator, y = Dhat, group = Perc_Forest_level)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = predator, shape = Perc_Forest_level), size = 2.75, position = position_dodge(width = 0.4)) +
  #'   scale_color_bright() + 
  #'   ylim(0,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none", shape = guide_legend(title = "Percent forested habitat")) + 
  #'   ggtitle("Effect of percent forest on summer predator-prey diel activity patterns") +
  #'   xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~prey, scales = "free", space = "free")
  #' pred_prey_for_smr_coeff_plot
  #' 
  #' pred_prey_for_fall_coeff_plot <- ggplot(pred_prey_for_overlap_tbl[pred_prey_for_overlap_tbl$season == "Fall",], aes(x = predator, y = Dhat, group = Perc_Forest_level)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = predator, shape = Perc_Forest_level), size = 2.75, position = position_dodge(width = 0.4)) +
  #'   scale_color_bright() + 
  #'   ylim(0,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none", shape = guide_legend(title = "Percent forested habitat")) + 
  #'   ggtitle("Effect of percent forest on fall predator-prey diel activity patterns") +
  #'   xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~prey, scales = "free", space = "free")
  #' pred_prey_for_fall_coeff_plot
  #' 
  #' pred_prey_for_wtr_coeff_plot <- ggplot(pred_prey_for_overlap_tbl[pred_prey_for_overlap_tbl$season == "Winter",], aes(x = predator, y = Dhat, group = Perc_Forest_level)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = predator, shape = Perc_Forest_level), size = 2.75, position = position_dodge(width = 0.4)) +
  #'   scale_color_bright() + 
  #'   ylim(0,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none", shape = guide_legend(title = "Percent forested habitat")) + 
  #'   ggtitle("Effect of percent forest on winter predator-prey diel activity patterns") +
  #'   xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~prey, scales = "free", space = "free")
  #' pred_prey_for_wtr_coeff_plot
  #' 
  #' pred_prey_for_sprg_coeff_plot <- ggplot(pred_prey_for_overlap_tbl[pred_prey_for_overlap_tbl$season == "Spring",], aes(x = predator, y = Dhat, group = Perc_Forest_level)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = predator, shape = Perc_Forest_level), size = 2.75, position = position_dodge(width = 0.4)) +
  #'   scale_color_bright() + 
  #'   ylim(0,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none", shape = guide_legend(title = "Percent forested habitat")) + 
  #'   ggtitle("Effect of percent forest on spring predator-prey diel activity patterns") +
  #'   xlab("Species pairing") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~prey, scales = "free", space = "free")
  #' pred_prey_for_sprg_coeff_plot
  #' 
  #' #' Save 'em - predator-prey overlap coefficients
  #' ggsave(pred_prey_for_smr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_PercForest_smr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(pred_prey_for_fall_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_PercForest_fall_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(pred_prey_for_wtr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_PercForest_wtr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(pred_prey_for_sprg_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_PercForest_sprg_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  
    
  ####  Effect of habitat complexity and predator ID on prey temporal overlap  ####
  #'  Pull out percent forest results
  prey_overlap_for_tbl <- filter(prey_overlap_tbl, Risk == "Percent Forest")
  prey_overlap_tbl <- filter(prey_overlap_tbl, Risk != "Percent Forest")
  
  #'  By Species and type of background risk
  prey_md_byrisk_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Species == "Mule deer",], aes(x = Season, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Season), width = 0.3) +
    geom_point(stat = "identity", aes(col = Season), size = 2.75) +
    scale_color_bright() + 
    ylim(0.5,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16), legend.text = element_text(size = 16)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Mule deer temporal overlap") +
    xlab("Season") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~Risk, scales = "free", space = "free")
  prey_md_byrisk_coeff_plot
  
  prey_elk_byrisk_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Species == "Elk",], aes(x = Season, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Season), width = 0.3) +
    geom_point(stat = "identity", aes(col = Season), size = 2.75) +
    scale_color_bright() + 
    ylim(0.5,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16), legend.text = element_text(size = 16)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Elk temporal overlap") +
    xlab("Season") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~Risk, scales = "free", space = "free")
  prey_elk_byrisk_coeff_plot
  
  prey_moose_byrisk_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Species == "Moose",], aes(x = Season, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Season), width = 0.3) +
    geom_point(stat = "identity", aes(col = Season), size = 2.75) +
    scale_color_bright() + 
    ylim(0.5,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16), legend.text = element_text(size = 16)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Moose temporal overlap") +
    xlab("Season") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~Risk, scales = "free", space = "free")
  prey_moose_byrisk_coeff_plot
  
  prey_wtd_byrisk_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Species == "White-tailed deer",], aes(x = Season, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Season), width = 0.3) +
    geom_point(stat = "identity", aes(col = Season), size = 2.75) +
    scale_color_bright() + 
    ylim(0.5,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16), legend.text = element_text(size = 16)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("White-tailed deer temporal overlap") +
    xlab("Season") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~Risk, scales = "free", space = "free")
  prey_wtd_byrisk_coeff_plot
  
  #' Save 'em - predator-prey overlap coefficients
  ggsave(prey_md_byrisk_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_md_byrisk_bckgrd_risk_skinny.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(prey_elk_byrisk_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_elk_byrisk_bckgrd_risk_skinny.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(prey_moose_byrisk_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_moose_byrisk_bckgrd_risk_skinny.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(prey_wtd_byrisk_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_wtd_byrisk_bckgrd_risk_skinny.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')

  
  #'  Just the percent forest effect
  prey_forest_coeff_plot <- ggplot(prey_overlap_for_tbl, aes(x = Season, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Season), width = 0.3) +
    geom_point(stat = "identity", aes(col = Season), size = 2.75) +
    scale_color_bright() + 
    ylim(0.5,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16), legend.text = element_text(size = 16)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Percent forest effect on temporal overlap") +
    xlab("Season") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~Species, scales = "free", space = "free")
  prey_forest_coeff_plot
  ggsave(prey_forest_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_prey_forest_bckgrd_risk_skinny.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  
  #' #'  By species
  #' prey_md_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Species == "Mule deer",], aes(x = Risk, y = Dhat)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
  #'   geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
  #'   scale_color_bright() + 
  #'   ylim(0.5,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none") + 
  #'   ggtitle("Mule deer temporal overlap") +
  #'   xlab("Type of background risk") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~Season, scales = "free", space = "free")
  #' prey_md_coeff_plot
  #' 
  #' prey_elk_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Species == "Elk",], aes(x = Risk, y = Dhat)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
  #'   geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
  #'   scale_color_bright() + 
  #'   ylim(0.5,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none") + 
  #'   ggtitle("Elk temporal overlap") +
  #'   xlab("Type of background risk") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~Season, scales = "free", space = "free")
  #' prey_elk_coeff_plot
  #' 
  #' prey_moose_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Species == "Moose",], aes(x = Risk, y = Dhat)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
  #'   geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
  #'   scale_color_bright() + 
  #'   ylim(0.5,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none") + 
  #'   ggtitle("Moose temporal overlap") +
  #'   xlab("Type of background risk") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~Season, scales = "free", space = "free")
  #' prey_moose_coeff_plot
  #' 
  #' prey_wtd_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Species == "White-tailed deer",], aes(x = Risk, y = Dhat)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
  #'   geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
  #'   scale_color_bright() + 
  #'   ylim(0.5,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none") + 
  #'   ggtitle("White-tailed deer temporal overlap") +
  #'   xlab("Type of background risk") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~Season, scales = "free", space = "free")
  #' prey_wtd_coeff_plot
  #' 
  #' #' Save 'em - prey overlap coefficients
  #' ggsave(prey_md_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_md_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(prey_elk_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_elk_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(prey_moose_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_moose_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(prey_wtd_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_wtd_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' 
  #' 
  #' #'  By season
  #' prey_smr_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Season == "Summer",], aes(x = Risk, y = Dhat)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
  #'   geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
  #'   scale_color_bright() + 
  #'   ylim(0.5,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none") + 
  #'   ggtitle("Prey activity in Summer") +
  #'   xlab("Type of background risk") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~Species, scales = "free", space = "free")
  #' prey_smr_coeff_plot
  #' 
  #' prey_fall_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Season == "Fall",], aes(x = Risk, y = Dhat)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
  #'   geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
  #'   scale_color_bright() + 
  #'   ylim(0.5,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none") + 
  #'   ggtitle("Prey activity in Fall") +
  #'   xlab("Type of background risk") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~Species, scales = "free", space = "free")
  #' prey_fall_coeff_plot
  #' 
  #' prey_wtr_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Season == "Winter",], aes(x = Risk, y = Dhat)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
  #'   geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
  #'   scale_color_bright() + 
  #'   ylim(0.5,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none") + 
  #'   ggtitle("Prey activity in Winter") +
  #'   xlab("Type of background risk") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~Species, scales = "free", space = "free")
  #' prey_wtr_coeff_plot
  #' 
  #' prey_sprg_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Season == "Spring",], aes(x = Risk, y = Dhat)) +
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
  #'   geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
  #'   scale_color_bright() + 
  #'   ylim(0.5,1) + theme_bw() +
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14), legend.text = element_text(size = 14)) + 
  #'   theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
  #'   guides(color = "none") + 
  #'   ggtitle("Prey activity in Spring") +
  #'   xlab("Type of background risk") + ylab("Coefficient of overlap (\u0394)") +
  #'   facet_grid(~Species, scales = "free", space = "free")
  #' prey_sprg_coeff_plot
  #' 
  #' #' Save 'em - prey overlap coefficients
  #' ggsave(prey_smr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_prey_smr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(prey_fall_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_prey_fall_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(prey_wtr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_prey_wtr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' ggsave(prey_sprg_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_prey_sprg_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ####  Effect of habitat complexity on predator temporal overlap  ####
  #'  Pull out percent forest results
  pred_overlap_for_tbl <- filter(pred_overlap_tbl, Risk == "Percent Forest")
  pred_overlap_tri_tbl <- filter(pred_overlap_tbl, Risk != "Percent Forest")
  
  #'  By type of background risk
  pred_for_coeff_plot <- ggplot(pred_overlap_for_tbl, aes(x = Season, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Season), width = 0.3) +
    geom_point(stat = "identity", aes(col = Season), size = 2.75) +
    scale_color_bright() + 
    ylim(0.5,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16), legend.text = element_text(size = 16)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Percent forest effect on temporal overlap") +
    xlab("Season") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~Species, scales = "free", space = "free")
  pred_for_coeff_plot
  ggsave(pred_for_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_forest_skinny.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  pred_tri_coeff_plot <- ggplot(pred_overlap_tri_tbl, aes(x = Season, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Season), width = 0.3) +
    geom_point(stat = "identity", aes(col = Season), size = 2.75) +
    scale_color_bright() + 
    ylim(0.5,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16), legend.text = element_text(size = 16)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Percent forest effect on temporal overlap") +
    xlab("Season") + ylab("Coefficient of overlap (\u0394)") +
    facet_grid(~Species, scales = "free", space = "free")
  pred_tri_coeff_plot
  ggsave(pred_tri_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_tri_skinny.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  
  #'  ------------------------------------------
  ####  Predator-Prey Activity Overlap Curves  ####
  #'  ------------------------------------------
  #'  Plot activity curves for each species-pairing at cameras with low versus 
  #'  high levels of background risk
  overlap_pred_prey_plots <- function(dat, name1, name2, name3, dhat, y_up, season) {
    #'  Sample sizes for predators[1] and prey[2] where background risk is low or high
    n1low <- dat[[7]]; n1high <- dat[[9]]
    n2low <- dat[[8]]; n2high <- dat[[10]]
    spp1low <- paste0(name1, " (n = ", n1low, ")"); spp1high <- paste0(name1, " (n = ", n1high, ")")
    spp2low <- paste0(name2, " (n = ", n2low, ")"); spp2high <- paste0(name2, " (n = ", n2high, ")")
    #'  Temporal overlap between predators and prey where background risk is low or high
    dhatlow <- dhat[1,6]; dhatlowl <- dhat[1,7]; dhatlowu<- dhat[1,8]
    dhathigh <- dhat[2,6]; dhathighl <- dhat[2,7]; dhathighu<- dhat[2,8]
    #'  Density data for overlap plots
    overdensity <- dat[[11]]
    #'  Separate data sets based on whether background risk is low or high
    lowrisk <- overdensity[overdensity$BackgroundRisk == "Low",]
    highrisk <- overdensity[overdensity$BackgroundRisk == "High",]
    
    overlap_low <- ggplot(lowrisk, aes(x, densityA, colour = Species.x)) +
      geom_line(lwd = 0.75) + 
      geom_line(aes(x, densityB, colour = Species.y), lwd = 0.75) +  
      geom_area(aes(y = pmin(densityA, densityB)),
                alpha = 0.3, color = NA) +
      # geom_line(aes(x, y, colour =  Species.z), linetype = "dashed", lwd = 0.75) +  
      scale_x_continuous(breaks = c(0, 1.57, 3.0, 4.71, 6.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = pi/2, linetype="dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype="dotted") +
      theme_bw() +
      theme(legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      ylim(0, y_up) +
      labs(x = "Time of day", y = "Density", color = paste0("\u0394 = ", dhatlow, " (", dhatlowl, " - ", dhatlowu, ")"), title = paste0("Low ", name3, ", ", season)) + 
      scale_color_manual(labels = c(spp1low, spp2low), values = c("red", "blue")) +
      theme(text = element_text(size = 14), legend.text = element_text(size = 14)) 
    plot(overlap_low)
    
    overlap_high <- ggplot(highrisk, aes(x, densityA, colour = Species.x)) +
      geom_line(lwd = 0.75) + 
      geom_line(aes(x, densityB, colour = Species.y), lwd = 0.75) +  
      geom_area(aes(y = pmin(densityA, densityB)),
                alpha = 0.3, color = NA) +
      scale_x_continuous(breaks = c(0, 1.57, 3.0, 4.71, 6.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = pi/2, linetype="dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype="dotted") +
      theme_bw() +
      theme(legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      ylim(0, y_up) +
      labs(x = "Time of day", y = "Density", color = paste0("\u0394 = ", dhathigh, " (", dhathighl, " - ", dhathighu, ")"), title = paste0("High ", name3, ", ", season)) + 
      scale_color_manual(labels = c(spp1high, spp2high), values = c("red", "blue")) +
      theme(text = element_text(size = 14), legend.text = element_text(size = 14)) 
    plot(overlap_high)
    
    plots <- list(overlap_low, overlap_high)
    return(plots)
  }
  ####  Predator-Prey TRI Overlap Plots  ####
  #'  Keep track of list positions when dhat1 and dhat4 are being combined
  #'  Dhat1 for sample sizes <50, Dhat4 for sample sizes >50, fig [[1]] = low, fig [[2]] = high
  ####  Cougar - mule deer TRI  ####
  coug_md_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[1]], name1 = "Cougar", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = coug_md_tri_smr_out, y_up = 0.6, season = "Summer")
  (coug_md_tri_smr_overlap_plot <- coug_md_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_md_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))  
  ggsave(coug_md_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_md_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_md_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[2]], name1 = "Cougar", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = coug_md_tri_fall_out, y_up = 0.6, season = "Fall")
  # (coug_md_tri_fall_overlap_plot <- coug_md_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_md_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_md_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_md_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_md_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[3]], name1 = "Cougar", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = coug_md_tri_wtr_out, y_up = 0.6, season = "Winter")
  # (coug_md_tri_wtr_overlap_plot <- coug_md_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_md_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_md_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_md_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_md_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[4]], name1 = "Cougar", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = coug_md_tri_sprg_out, y_up = 0.6, season  = "Spring")
  (coug_md_tri_sprg_overlap_plot <- coug_md_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_md_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_md_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_md_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Cougar - elk TRI  ####
  #' Low risk >50 cougars & elk; High risk <50 cougars & elk
  coug_elk_tri_smr_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[5]], name1 = "Cougar", name2 = "Elk", name3 = "terrain ruggedness", dhat = coug_elk_tri_smr_out, y_up = 0.6, season = "Summer")
  coug_elk_tri_smr_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[6]], name1 = "Cougar", name2 = "Elk", name3 = "terrain ruggedness", dhat = coug_elk_tri_smr_out, y_up = 0.6, season = "Summer")
  (coug_elk_tri_smr_overlap_plot <- coug_elk_tri_smr_overPlot_g4[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_elk_tri_smr_overPlot_g1[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_elk_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_elk_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_elk_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[7]], name1 = "Cougar", name2 = "Elk", name3 = "terrain ruggedness", dhat = coug_elk_tri_fall_out, y_up = 0.6, season = "Fall")
  # (coug_elk_tri_fall_overlap_plot <- coug_elk_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_elk_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_elk_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_elk_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_elk_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[8]], name1 = "Cougar", name2 = "Elk", name3 = "terrain ruggedness", dhat = coug_elk_tri_wtr_out, y_up = 0.6, season = "Winter")
  # (coug_elk_tri_wtr_overlap_plot <- coug_elk_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_elk_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_elk_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_elk_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_elk_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[9]], name1 = "Cougar", name2 = "Elk", name3 = "terrain ruggedness", dhat = coug_elk_tri_sprg_out, y_up = 0.6, season  = "Spring")
  # (coug_elk_tri_sprg_overlap_plot <- coug_elk_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_elk_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_elk_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_elk_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Cougar - moose TRI  ####
  coug_moose_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[10]], name1 = "Cougar", name2 = "Moose", name3 = "terrain ruggedness", dhat = coug_moose_tri_smr_out, y_up = 0.6, season = "Summer")
  (coug_moose_tri_smr_overlap_plot <- coug_moose_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_moose_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_moose_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_moose_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_moose_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[11]], name1 = "Cougar", name2 = "Moose", name3 = "terrain ruggedness", dhat = coug_moose_tri_fall_out, y_up = 0.6, season = "Fall")
  # (coug_moose_tri_fall_overlap_plot <- coug_moose_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_moose_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_moose_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_moose_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_moose_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[12]], name1 = "Cougar", name2 = "Moose", name3 = "terrain ruggedness", dhat = coug_moose_tri_wtr_out, y_up = 0.6, season = "Winter")
  # (coug_moose_tri_wtr_overlap_plot <- coug_moose_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_moose_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_moose_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_moose_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_moose_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[13]], name1 = "Cougar", name2 = "Moose", name3 = "terrain ruggedness", dhat = coug_moose_tri_sprg_out, y_up = 0.6, season  = "Spring")
  (coug_moose_tri_sprg_overlap_plot <- coug_moose_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + coug_moose_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_moose_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_moose_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Cougar - white-tailed deer TRI  ####
  coug_wtd_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[14]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coug_wtd_tri_smr_out, y_up = 0.6, season = "Summer")
  (coug_wtd_tri_smr_overlap_plot <- coug_wtd_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coug_wtd_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_wtd_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtd_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_wtd_tri_fall_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[15]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coug_wtd_tri_fall_out, y_up = 0.6, season = "Fall")
  coug_wtd_tri_fall_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[16]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coug_wtd_tri_fall_out, y_up = 0.6, season = "Fall")
  (coug_wtd_tri_fall_overlap_plot <- coug_wtd_tri_fall_overPlot_g4[[1]] + theme(legend.position = c(0.40, 0.895)) + coug_wtd_tri_fall_overPlot_g1[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_wtd_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtd_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_wtd_tri_wtr_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[17]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coug_wtd_tri_wtr_out, y_up = 0.6, season = "Winter")
  coug_wtd_tri_wtr_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[18]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coug_wtd_tri_wtr_out, y_up = 0.6, season = "Winter")
  (coug_wtd_tri_wtr_overlap_plot <- coug_wtd_tri_wtr_overPlot_g4[[1]] + theme(legend.position = c(0.40, 0.895)) + coug_wtd_tri_wtr_overPlot_g1[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_wtd_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtd_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #' #'  Low risk >50 cougars; High risk < 50 cougars
  #' coug_wtd_tri_sprg_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[19]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coug_wtd_tri_sprg_out, y_up = 0.6, season  = "Spring")
  #' coug_wtd_tri_sprg_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[1]][[20]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coug_wtd_tri_sprg_out, y_up = 0.6, season  = "Spring")
  #' (coug_wtd_tri_sprg_overlap_plot <- coug_wtd_tri_sprg_overPlot_g4[[1]] + theme(legend.position = c(0.40, 0.895)) + coug_wtd_tri_sprg_overPlot_g1[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  #' ggsave(coug_wtd_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtd_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Wolf - mule deer TRI  ####
  wolf_md_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[1]], name1 = "Wolf", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = wolf_md_tri_smr_out, y_up = 0.6, season = "Summer")
  (wolf_md_tri_smr_overlap_plot <- wolf_md_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_md_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(wolf_md_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_md_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_md_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[2]], name1 = "Wolf", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = wolf_md_tri_fall_out, y_up = 0.6, season = "Fall")
  # (wolf_md_tri_fall_overlap_plot <- wolf_md_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_md_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_md_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_md_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_md_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[3]], name1 = "Wolf", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = wolf_md_tri_wtr_out, y_up = 0.6, season = "Winter")
  # (wolf_md_tri_wtr_overlap_plot <- wolf_md_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_md_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_md_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_md_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_md_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[4]], name1 = "Wolf", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = wolf_md_tri_sprg_out, y_up = 0.65, season  = "Spring")
  # (wolf_md_tri_sprg_overlap_plot <- wolf_md_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_md_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_md_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_md_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Wolf - elk TRI  ####
  # wolf_elk_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[5]], name1 = "Wolf", name2 = "Elk", name3 = "terrain ruggedness", dhat = wolf_elk_tri_smr_out, y_up = 0.6, season = "Summer")
  # (wolf_elk_tri_smr_overlap_plot <- wolf_elk_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_elk_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_elk_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_elk_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Wolf - moose TRI  ####
  wolf_moose_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[6]], name1 = "Wolf", name2 = "Moose", name3 = "terrain ruggedness", dhat = wolf_moose_tri_smr_out, y_up = 0.6, season = "Summer")
  (wolf_moose_tri_smr_overlap_plot <- wolf_moose_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_moose_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(wolf_moose_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_moose_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_moose_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[7]], name1 = "Wolf", name2 = "Moose", name3 = "terrain ruggedness", dhat = wolf_moose_tri_fall_out, y_up = 0.6, season = "Fall")
  # (wolf_moose_tri_fall_overlap_plot <- wolf_moose_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_moose_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_moose_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_moose_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_moose_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[8]], name1 = "Wolf", name2 = "Moose", name3 = "terrain ruggedness", dhat = wolf_moose_tri_wtr_out, y_up = 0.6, season = "Winter")
  # (wolf_moose_tri_wtr_overlap_plot <- wolf_moose_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_moose_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_moose_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_moose_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_moose_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[9]], name1 = "Wolf", name2 = "Moose", name3 = "terrain ruggedness", dhat = wolf_moose_tri_sprg_out, y_up = 0.6, season  = "Spring")
  # (wolf_moose_tri_sprg_overlap_plot <- wolf_moose_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.35, 0.895)) + wolf_moose_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.35, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_moose_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_moose_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Wolf - white-tailed deer TRI  ####
  wolf_wtd_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[10]], name1 = "Wolf", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = wolf_wtd_tri_smr_out, y_up = 0.6, season = "Summer")
  (wolf_wtd_tri_smr_overlap_plot <- wolf_wtd_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + wolf_wtd_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(wolf_wtd_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtd_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_wtd_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[11]], name1 = "Wolf", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = wolf_wtd_tri_fall_out, y_up = 0.6, season = "Fall")
  # (wolf_wtd_tri_fall_overlap_plot <- wolf_wtd_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + wolf_wtd_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_wtd_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtd_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_wtd_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[12]], name1 = "Wolf", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = wolf_wtd_tri_wtr_out, y_up = 0.6, season = "Winter")
  # (wolf_wtd_tri_wtr_overlap_plot <- wolf_wtd_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + wolf_wtd_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_wtd_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtd_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_wtd_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[2]][[13]], name1 = "Wolf", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = wolf_wtd_tri_sprg_out, y_up = 0.65, season  = "Spring")
  # (wolf_wtd_tri_sprg_overlap_plot <- wolf_wtd_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + wolf_wtd_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_wtd_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtd_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ####  Black bear - mule deer TRI  ####
  bear_md_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[1]], name1 = "Black bear", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = bear_md_tri_smr_out, y_up = 0.6, season = "Summer")
  (bear_md_tri_smr_overlap_plot <- bear_md_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_md_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_md_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_md_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk 50 black bears; High risk <50 black bears
  bear_md_tri_fall_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[2]], name1 = "Black bear", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = bear_md_tri_fall_out, y_up = 0.6, season = "Fall")
  bear_md_tri_fall_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[3]], name1 = "Black bear", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = bear_md_tri_fall_out, y_up = 0.6, season = "Fall")
  (bear_md_tri_fall_overlap_plot <- bear_md_tri_fall_overPlot_g4[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_md_tri_fall_overPlot_g1[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_md_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_md_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bear_md_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[4]], name1 = "Black bear", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = bear_md_tri_sprg_out, y_up = 0.65, season  = "Spring")
  (bear_md_tri_sprg_overlap_plot <- bear_md_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_md_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_md_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_md_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Black bear - elk TRI  ####
  bear_elk_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[5]], name1 = "Black bear", name2 = "Elk", name3 = "terrain ruggedness", dhat = bear_elk_tri_smr_out, y_up = 0.6, season = "Summer")
  (bear_elk_tri_smr_overlap_plot <- bear_elk_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_elk_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_elk_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_elk_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # bear_elk_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[6]], name1 = "Black bear", name2 = "Elk", name3 = "terrain ruggedness", dhat = bear_elk_tri_fall_out, y_up = 0.6, season = "Fall")
  # (bear_elk_tri_fall_overlap_plot <- bear_elk_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_elk_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(bear_elk_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_elk_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # bear_elk_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[7]], name1 = "Black bear", name2 = "Elk", name3 = "terrain ruggedness", dhat = bear_elk_tri_sprg_out, y_up = 0.65, season  = "Spring")
  # (bear_elk_tri_sprg_overlap_plot <- bear_elk_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_elk_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(bear_elk_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_elk_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Black bear - moose TRI  ####
  bear_moose_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[8]], name1 = "Black bear", name2 = "Moose", name3 = "terrain ruggedness", dhat = bear_moose_tri_smr_out, y_up = 0.6, season = "Summer")
  (bear_moose_tri_smr_overlap_plot <- bear_moose_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_moose_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_moose_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_moose_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # bear_moose_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[9]], name1 = "Black bear", name2 = "Moose", name3 = "terrain ruggedness", dhat = bear_moose_tri_fall_out, y_up = 0.6, season = "Fall")
  # (bear_moose_tri_fall_overlap_plot <- bear_moose_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_moose_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(bear_moose_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_moose_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bear_moose_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[10]], name1 = "Black bear", name2 = "Moose", name3 = "terrain ruggedness", dhat = bear_moose_tri_sprg_out, y_up = 0.65, season  = "Spring")
  (bear_moose_tri_sprg_overlap_plot <- bear_moose_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_moose_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_moose_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_moose_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Black bear - white-tailed deer TRI  ####
  bear_wtd_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[11]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bear_wtd_tri_smr_out, y_up = 0.6, season = "Summer")
  (bear_wtd_tri_smr_overlap_plot <- bear_wtd_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + bear_wtd_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_wtd_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_wtd_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #' bear_wtd_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[12]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bear_wtd_tri_fall_out, y_up = 0.6, season = "Fall")
  #' (bear_wtd_tri_fall_overlap_plot <- bear_wtd_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + bear_wtd_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  #' ggsave(bear_wtd_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_wtd_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #' #'  Low risk >50 black bears; High risk <50 black bears
  #' bear_wtd_tri_sprg_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[13]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bear_wtd_tri_sprg_out, y_up = 0.65, season  = "Spring")
  #' bear_wtd_tri_sprg_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[3]][[14]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bear_wtd_tri_sprg_out, y_up = 0.65, season  = "Spring")
  #' (bear_wtd_tri_sprg_overlap_plot <- bear_wtd_tri_sprg_overPlot_g4[[1]] + theme(legend.position = c(0.40, 0.895)) + bear_wtd_tri_sprg_overPlot_g1[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  #' ggsave(bear_wtd_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_wtd_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Bobcat - mule deer TRI  ####
  bob_md_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[1]], name1 = "Bobcat", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = bob_md_tri_smr_out, y_up = 0.6, season = "Summer")
  (bob_md_tri_smr_overlap_plot <- bob_md_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.28, 0.895)) + bob_md_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.28, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_md_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_md_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bob_md_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[2]], name1 = "Bobcat", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = bob_md_tri_fall_out, y_up = 0.6, season = "Fall")
  (bob_md_tri_fall_overlap_plot <- bob_md_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.28, 0.895)) + bob_md_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.28, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_md_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_md_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # bob_md_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[3]], name1 = "Bobcat", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = bob_md_tri_wtr_out, y_up = 0.6, season = "Winter")
  # (bob_md_tri_wtr_overlap_plot <- bob_md_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.28, 0.895)) + bob_md_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.28, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(bob_md_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_md_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bob_md_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[4]], name1 = "Bobcat", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = bob_md_tri_sprg_out, y_up = 0.65, season  = "Spring")
  (bob_md_tri_sprg_overlap_plot <- bob_md_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.28, 0.895)) + bob_md_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.28, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_md_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_md_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Bobcat - white-tailed deer TRI  ####
  bob_wtd_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[5]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bob_wtd_tri_smr_out, y_up = 0.6, season = "Summer")
  (bob_wtd_tri_smr_overlap_plot <- bob_wtd_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + bob_wtd_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_wtd_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtd_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk >50 bobcat; High risk <50 bobcat
  bob_wtd_tri_fall_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[6]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bob_wtd_tri_fall_out, y_up = 0.6, season = "Fall")
  bob_wtd_tri_fall_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[7]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bob_wtd_tri_fall_out, y_up = 0.6, season = "Fall")
  (bob_wtd_tri_fall_overlap_plot <- bob_wtd_tri_fall_overPlot_g4[[1]] + theme(legend.position = c(0.40, 0.895)) + bob_wtd_tri_fall_overPlot_g1[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_wtd_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtd_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk >50 bobcat; High risk <50 bobcat
  bob_wtd_tri_wtr_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[8]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bob_wtd_tri_wtr_out, y_up = 0.6, season = "Winter")
  bob_wtd_tri_wtr_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[9]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bob_wtd_tri_wtr_out, y_up = 0.6, season = "Winter")
  (bob_wtd_tri_wtr_overlap_plot <- bob_wtd_tri_wtr_overPlot_g4[[1]] + theme(legend.position = c(0.40, 0.895)) + bob_wtd_tri_wtr_overPlot_g1[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_wtd_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtd_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bob_wtd_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[4]][[10]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = bob_wtd_tri_sprg_out, y_up = 0.65, season  = "Spring")
  (bob_wtd_tri_sprg_overlap_plot <- bob_wtd_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + bob_wtd_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_wtd_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtd_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Coyote - mule deer TRI  ####
  coy_md_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[5]][[1]], name1 = "Coyote", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = coy_md_tri_smr_out, y_up = 0.6, season = "Summer")
  (coy_md_tri_smr_overlap_plot <- coy_md_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coy_md_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_md_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_md_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_md_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[5]][[2]], name1 = "Coyote", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = coy_md_tri_fall_out, y_up = 0.6, season = "Fall")
  (coy_md_tri_fall_overlap_plot <- coy_md_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coy_md_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_md_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_md_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_md_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[5]][[3]], name1 = "Coyote", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = coy_md_tri_wtr_out, y_up = 0.6, season = "Winter")
  (coy_md_tri_wtr_overlap_plot <- coy_md_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coy_md_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_md_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_md_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_md_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[5]][[4]], name1 = "Coyote", name2 = "Mule deer", name3 = "terrain ruggedness", dhat = coy_md_tri_sprg_out, y_up = 0.65, season  = "Spring")
  (coy_md_tri_sprg_overlap_plot <- coy_md_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coy_md_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_md_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_md_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Coyote - white-tailed deer TRI  ####
  coy_wtd_tri_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[5]][[5]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coy_wtd_tri_smr_out, y_up = 0.6, season = "Summer")
  (coy_wtd_tri_smr_overlap_plot <- coy_wtd_tri_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coy_wtd_tri_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_wtd_tri_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtd_tri_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtd_tri_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[5]][[6]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coy_wtd_tri_fall_out, y_up = 0.6, season = "Fall")
  (coy_wtd_tri_fall_overlap_plot <- coy_wtd_tri_fall_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coy_wtd_tri_fall_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_wtd_tri_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtd_tri_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtd_tri_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[5]][[7]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coy_wtd_tri_wtr_out, y_up = 0.6, season = "Winter")
  (coy_wtd_tri_wtr_overlap_plot <- coy_wtd_tri_wtr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coy_wtd_tri_wtr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_wtd_tri_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtd_tri_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtd_tri_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_tri_overlap[[5]][[8]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "terrain ruggedness", dhat = coy_wtd_tri_sprg_out, y_up = 0.65, season  = "Spring")
  (coy_wtd_tri_sprg_overlap_plot <- coy_wtd_tri_sprg_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coy_wtd_tri_sprg_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_wtd_tri_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtd_tri_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ####  Predator-Prey % Forest Overlap Plots  ####
  #'  Keep track of list positions when dhat1 and dhat4 are being combined
  #'  Dhat1 for sample sizes <50, Dhat4 for sample sizes >50, fig [[1]] = low, fig [[2]] = high
  ####  Cougar - mule deer % forest  ####
  coug_md_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[1]], name1 = "Cougar", name2 = "Mule deer", name3 = "percent forest", dhat = coug_md_for_smr_out, y_up = 0.6, season = "Summer")
  (coug_md_for_smr_overlap_plot <- coug_md_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_md_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_md_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_md_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_md_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[2]], name1 = "Cougar", name2 = "Mule deer", name3 = "percent forest", dhat = coug_md_for_fall_out, y_up = 0.6, season = "Fall")
  # (coug_md_for_fall_overlap_plot <- coug_md_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_md_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_md_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_md_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_md_for_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[3]], name1 = "Cougar", name2 = "Mule deer", name3 = "percent forest", dhat = coug_md_for_wtr_out, y_up = 0.6, season = "Winter")
  # (coug_md_for_wtr_overlap_plot <- coug_md_for_wtr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_md_for_wtr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_md_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_md_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_md_for_sprg_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[4]], name1 = "Cougar", name2 = "Mule deer", name3 = "percent forest", dhat = coug_md_for_sprg_out, y_up = 0.6, season  = "Spring")
  coug_md_for_sprg_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[5]], name1 = "Cougar", name2 = "Mule deer", name3 = "percent forest", dhat = coug_md_for_sprg_out, y_up = 0.6, season  = "Spring")
  (coug_md_for_sprg_overlap_plot <- coug_md_for_sprg_overPlot_g4[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_md_for_sprg_overPlot_g1[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_md_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_md_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Cougar - elk % forest  ####
  #' Low risk <50 cougars; High risk >50 cougars
  coug_elk_for_smr_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[6]], name1 = "Cougar", name2 = "Elk", name3 = "percent forest", dhat = coug_elk_for_smr_out, y_up = 0.6, season = "Summer")
  coug_elk_for_smr_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[7]], name1 = "Cougar", name2 = "Elk", name3 = "percent forest", dhat = coug_elk_for_smr_out, y_up = 0.6, season = "Summer")
  (coug_elk_for_smr_overlap_plot <- coug_elk_for_smr_overPlot_g1[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_elk_for_smr_overPlot_g4[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_elk_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_elk_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_elk_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[8]], name1 = "Cougar", name2 = "Elk", name3 = "percent forest", dhat = coug_elk_for_fall_out, y_up = 0.6, season = "Fall")
  # (coug_elk_for_fall_overlap_plot <- coug_elk_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_elk_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_elk_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_elk_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_elk_for_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[9]], name1 = "Cougar", name2 = "Elk", name3 = "percent forest", dhat = coug_elk_for_wtr_out, y_up = 0.6, season = "Winter")
  # (coug_elk_for_wtr_overlap_plot <- coug_elk_for_wtr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_elk_for_wtr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_elk_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_elk_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_elk_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[10]], name1 = "Cougar", name2 = "Elk", name3 = "percent forest", dhat = coug_elk_for_sprg_out, y_up = 0.6, season = "Spring")
  (coug_elk_for_sprg_overlap_plot <- coug_elk_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_elk_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_elk_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_elk_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Cougar - moose % forest  ####
  coug_moose_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[11]], name1 = "Cougar", name2 = "Moose", name3 = "percent forest", dhat = coug_moose_for_smr_out, y_up = 0.6, season = "Summer")
  (coug_moose_for_smr_overlap_plot <- coug_moose_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_moose_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_moose_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_moose_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_moose_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[12]], name1 = "Cougar", name2 = "Moose", name3 = "percent forest", dhat = coug_moose_for_fall_out, y_up = 0.6, season = "Fall")
  (coug_moose_for_fall_overlap_plot <- coug_moose_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_moose_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_moose_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_moose_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coug_moose_for_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[13]], name1 = "Cougar", name2 = "Moose", name3 = "percent forest", dhat = coug_moose_for_wtr_out, y_up = 0.6, season = "Winter")
  # (coug_moose_for_wtr_overlap_plot <- coug_moose_for_wtr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_moose_for_wtr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coug_moose_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_moose_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_moose_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[14]], name1 = "Cougar", name2 = "Moose", name3 = "percent forest", dhat = coug_moose_for_sprg_out, y_up = 0.6, season = "Spring")
  (coug_moose_for_sprg_overlap_plot <- coug_moose_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coug_moose_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_moose_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_moose_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Cougar - white-tailed deer % forest  ####
  coug_wtd_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[15]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "percent forest", dhat = coug_wtd_for_smr_out, y_up = 0.6, season = "Summer")
  (coug_wtd_for_smr_overlap_plot <- coug_wtd_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coug_wtd_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_wtd_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtd_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk <50 cougars; High risk >50 cougars  
  coug_wtd_for_fall_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[16]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "percent forest", dhat = coug_wtd_for_fall_out, y_up = 0.6, season = "Fall")
  coug_wtd_for_fall_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[17]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "percent forest", dhat = coug_wtd_for_fall_out, y_up = 0.6, season = "Fall")
  (coug_wtd_for_fall_overlap_plot <- coug_wtd_for_fall_overPlot_g1[[1]] + theme(legend.position = c(0.40, 0.895)) + coug_wtd_for_fall_overPlot_g4[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_wtd_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtd_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_wtd_for_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[18]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "percent forest", dhat = coug_wtd_for_wtr_out, y_up = 0.6, season = "Winter")
  (coug_wtd_for_wtr_overlap_plot <- coug_wtd_for_wtr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coug_wtd_for_wtr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_wtd_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtd_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coug_wtd_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[1]][[19]], name1 = "Cougar", name2 = "White-tailed deer", name3 = "percent forest", dhat = coug_wtd_for_sprg_out, y_up = 0.6, season = "Spring")
  (coug_wtd_for_sprg_overlap_plot <- coug_wtd_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coug_wtd_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coug_wtd_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtd_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Wolf - mule deer % forest  ####
  # wolf_md_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[1]], name1 = "Wolf", name2 = "Mule deer", name3 = "percent forest", dhat = wolf_md_for_smr_out, y_up = 0.6, season = "Summer")
  # (wolf_md_for_smr_overlap_plot <- wolf_md_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + wolf_md_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_md_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_md_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_md_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[2]], name1 = "Wolf", name2 = "Mule deer", name3 = "percent forest", dhat = wolf_md_for_fall_out, y_up = 0.6, season = "Fall")
  # (wolf_md_for_fall_overlap_plot <- wolf_md_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + wolf_md_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_md_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_md_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_md_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[3]], name1 = "Wolf", name2 = "Mule deer", name3 = "percent forest", dhat = wolf_md_for_sprg_out, y_up = 0.6, season = "Spring")
  # (wolf_md_for_sprg_overlap_plot <- wolf_md_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + wolf_md_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_md_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_md_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Wolf - elk % forest  ####
  # wolf_elk_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[4]], name1 = "Wolf", name2 = "Elk", name3 = "percent forest", dhat = wolf_elk_for_smr_out, y_up = 0.6, season = "Summer")
  # (wolf_elk_for_smr_overlap_plot <- wolf_elk_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + wolf_elk_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_elk_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_elk_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Wolf - moose % forest  ####
  wolf_moose_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[5]], name1 = "Wolf", name2 = "Moose", name3 = "percent forest", dhat = wolf_moose_for_smr_out, y_up = 0.6, season = "Summer")
  (wolf_moose_for_smr_overlap_plot <- wolf_moose_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + wolf_moose_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(wolf_moose_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_moose_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_moose_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[6]], name1 = "Wolf", name2 = "Moose", name3 = "percent forest", dhat = wolf_moose_for_fall_out, y_up = 0.6, season = "Fall")
  # (wolf_moose_for_fall_overlap_plot <- wolf_moose_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + wolf_moose_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_moose_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_moose_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_moose_for_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[7]], name1 = "Wolf", name2 = "Moose", name3 = "percent forest", dhat = wolf_moose_for_wtr_out, y_up = 0.6, season = "Winter")
  # (wolf_moose_for_wtr_overlap_plot <- wolf_moose_for_wtr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + wolf_moose_for_wtr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_moose_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_moose_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Wolf - white-tailed deer % forest  ####
  # wolf_wtd_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[8]], name1 = "Wolf", name2 = "White-tailed deer", name3 = "percent forest", dhat = wolf_wtd_for_smr_out, y_up = 0.6, season = "Summer")
  # (wolf_wtd_for_smr_overlap_plot <- wolf_wtd_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + wolf_wtd_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_wtd_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtd_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_wtd_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[9]], name1 = "Wolf", name2 = "White-tailed deer", name3 = "percent forest", dhat = wolf_wtd_for_fall_out, y_up = 0.6, season = "Fall")
  # (wolf_wtd_for_fall_overlap_plot <- wolf_wtd_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + wolf_wtd_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_wtd_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtd_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # wolf_wtd_for_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[2]][[10]], name1 = "Wolf", name2 = "White-tailed deer", name3 = "percent forest", dhat = wolf_wtd_for_wtr_out, y_up = 0.6, season = "Winter")
  # (wolf_wtd_for_wtr_overlap_plot <- wolf_wtd_for_wtr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + wolf_wtd_for_wtr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(wolf_wtd_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtd_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Black  bear - mule deer % forest  ####
  bear_md_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[1]], name1 = "Black bear", name2 = "Mule deer", name3 = "percent forest", dhat = bear_md_for_smr_out, y_up = 0.6, season = "Summer")
  (bear_md_for_smr_overlap_plot <- bear_md_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_md_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_md_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_md_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk >50 black bears; High risk <50 black bears
  bear_md_for_fall_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[2]], name1 = "Black bear", name2 = "Mule deer", name3 = "percent forest", dhat = bear_md_for_fall_out, y_up = 0.6, season = "Fall")
  bear_md_for_fall_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[3]], name1 = "Black bear", name2 = "Mule deer", name3 = "percent forest", dhat = bear_md_for_fall_out, y_up = 0.6, season = "Fall")
  (bear_md_for_fall_overlap_plot <- bear_md_for_fall_overPlot_g4[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_md_for_fall_overPlot_g1[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_md_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_md_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bear_md_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[4]], name1 = "Black bear", name2 = "Mule deer", name3 = "percent forest", dhat = bear_md_for_sprg_out, y_up = 0.6, season = "Spring")
  (bear_md_for_sprg_overlap_plot <- bear_md_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_md_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_md_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_md_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Black  bear - elk % forest  ####
  bear_elk_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[5]], name1 = "Black bear", name2 = "Elk", name3 = "percent forest", dhat = bear_elk_for_smr_out, y_up = 0.6, season = "Summer")
  (bear_elk_for_smr_overlap_plot <- bear_elk_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_elk_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_elk_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_elk_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # bear_elk_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[6]], name1 = "Black bear", name2 = "Elk", name3 = "percent forest", dhat = bear_elk_for_fall_out, y_up = 0.6, season = "Fall")
  # (bear_elk_for_fall_overlap_plot <- bear_elk_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_elk_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(bear_elk_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_elk_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # bear_elk_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[7]], name1 = "Black bear", name2 = "Elk", name3 = "percent forest", dhat = bear_elk_for_sprg_out, y_up = 0.6, season = "Spring")
  # (bear_elk_for_sprg_overlap_plot <- bear_elk_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_elk_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(bear_elk_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_elk_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Black  bear - moose % forest  ####
  bear_moose_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[8]], name1 = "Black bear", name2 = "Moose", name3 = "percent forest", dhat = bear_moose_for_smr_out, y_up = 0.6, season = "Summer")
  (bear_moose_for_smr_overlap_plot <- bear_moose_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_moose_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_moose_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_moose_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # bear_moose_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[9]], name1 = "Black bear", name2 = "Moose", name3 = "percent forest", dhat = bear_moose_for_fall_out, y_up = 0.65, season = "Fall")
  # (bear_moose_for_fall_overlap_plot <- bear_moose_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_moose_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(bear_moose_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_moose_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bear_moose_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[10]], name1 = "Black bear", name2 = "Moose", name3 = "percent forest", dhat = bear_moose_for_sprg_out, y_up = 0.6, season = "Spring")
  (bear_moose_for_sprg_overlap_plot <- bear_moose_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bear_moose_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_moose_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_moose_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Black  bear - white-tailed deer % forest  ####
  bear_wtd_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[11]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "percent forest", dhat = bear_wtd_for_smr_out, y_up = 0.6, season = "Summer")
  (bear_wtd_for_smr_overlap_plot <- bear_wtd_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + bear_wtd_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_wtd_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_wtd_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bear_wtd_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[12]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "percent forest", dhat = bear_wtd_for_fall_out, y_up = 0.65, season = "Fall")
  (bear_wtd_for_fall_overlap_plot <- bear_wtd_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + bear_wtd_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_wtd_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_wtd_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk <50 black bears; High risk >50 black bears
  bear_wtd_for_sprg_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[13]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "percent forest", dhat = bear_wtd_for_sprg_out, y_up = 0.6, season = "Spring")
  bear_wtd_for_sprg_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_for_overlap[[3]][[14]], name1 = "Black bear", name2 = "White-tailed deer", name3 = "percent forest", dhat = bear_wtd_for_sprg_out, y_up = 0.6, season = "Spring")
  (bear_wtd_for_sprg_overlap_plot <- bear_wtd_for_sprg_overPlot_g1[[1]] + theme(legend.position = c(0.40, 0.895)) + bear_wtd_for_sprg_overPlot_g4[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bear_wtd_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_wtd_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Bobcat - mule deer % forest  ####
  bob_md_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[1]], name1 = "Bobcat", name2 = "Mule deer", name3 = "percent forest", dhat = bob_md_for_smr_out, y_up = 0.6, season = "Summer")
  (bob_md_for_smr_overlap_plot <- bob_md_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bob_md_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_md_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_md_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bob_md_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[2]], name1 = "Bobcat", name2 = "Mule deer", name3 = "percent forest", dhat = bob_md_for_fall_out, y_up = 0.65, season = "Fall")
  (bob_md_for_fall_overlap_plot <- bob_md_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bob_md_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_md_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_md_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bob_md_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[3]], name1 = "Bobcat", name2 = "Mule deer", name3 = "percent forest", dhat = bob_md_for_sprg_out, y_up = 0.6, season = "Spring")
  (bob_md_for_sprg_overlap_plot <- bob_md_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + bob_md_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_md_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_md_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Bobcat - white-tailed deer % forest  ####
  bob_wtd_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[4]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "percent forest", dhat = bob_wtd_for_smr_out, y_up = 0.6, season = "Summer")
  (bob_wtd_for_smr_overlap_plot <- bob_wtd_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + bob_wtd_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_wtd_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtd_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk <50 bobcat; High risk >50 bobcat
  bob_wtd_for_fall_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[5]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "percent forest", dhat = bob_wtd_for_fall_out, y_up = 0.65, season = "Fall")
  bob_wtd_for_fall_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[6]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "percent forest", dhat = bob_wtd_for_fall_out, y_up = 0.65, season = "Fall")
  (bob_wtd_for_fall_overlap_plot <- bob_wtd_for_fall_overPlot_g1[[1]] + theme(legend.position = c(0.40, 0.895)) + bob_wtd_for_fall_overPlot_g4[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_wtd_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtd_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  #'  Low risk <50 bobcats; High risk >50 bobcats
  bob_wtd_for_wtr_overPlot_g1 <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[7]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "percent forest", dhat = bob_wtd_for_wtr_out, y_up = 0.6, season = "Winter")
  bob_wtd_for_wtr_overPlot_g4 <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[8]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "percent forest", dhat = bob_wtd_for_wtr_out, y_up = 0.6, season = "Winter")
  (bob_wtd_for_wtr_overlap_plot <- bob_wtd_for_wtr_overPlot_g1[[1]] + theme(legend.position = c(0.40, 0.895)) + bob_wtd_for_wtr_overPlot_g4[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_wtd_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtd_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  bob_wtd_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[4]][[9]], name1 = "Bobcat", name2 = "White-tailed deer", name3 = "percent forest", dhat = bob_wtd_for_sprg_out, y_up = 0.6, season = "Spring")
  (bob_wtd_for_sprg_overlap_plot <- bob_wtd_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + bob_wtd_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(bob_wtd_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtd_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Coyote - mule deer % forest  ####
  coy_md_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[5]][[1]], name1 = "Coyote", name2 = "Mule deer", name3 = "percent forest", dhat = coy_md_for_smr_out, y_up = 0.6, season = "Summer")
  (coy_md_for_smr_overlap_plot <- coy_md_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coy_md_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_md_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_md_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_md_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[5]][[2]], name1 = "Coyote", name2 = "Mule deer", name3 = "percent forest", dhat = coy_md_for_fall_out, y_up = 0.6, season = "Fall")
  (coy_md_for_fall_overlap_plot <- coy_md_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coy_md_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_md_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_md_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  # coy_md_for_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[5]][[3]], name1 = "Coyote", name2 = "Mule deer", name3 = "percent forest", dhat = coy_md_for_wtr_out, y_up = 0.65, season = "Winter")
  # (coy_md_for_wtr_overlap_plot <- coy_md_for_wtr_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coy_md_for_wtr_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  # ggsave(coy_md_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_md_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_md_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[5]][[4]], name1 = "Coyote", name2 = "Mule deer", name3 = "percent forest", dhat = coy_md_for_sprg_out, y_up = 0.6, season = "Spring")
  (coy_md_for_sprg_overlap_plot <- coy_md_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.30, 0.895)) + coy_md_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.30, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_md_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_md_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  ####  Coyote - white-tailed deer % forest  ####
  coy_wtd_for_smr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[5]][[5]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "percent forest", dhat = coy_wtd_for_smr_out, y_up = 0.6, season = "Summer")
  (coy_wtd_for_smr_overlap_plot <- coy_wtd_for_smr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coy_wtd_for_smr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_wtd_for_smr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtd_for_smr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtd_for_fall_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[5]][[6]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "percent forest", dhat = coy_wtd_for_fall_out, y_up = 0.6, season = "Fall")
  (coy_wtd_for_fall_overlap_plot <- coy_wtd_for_fall_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coy_wtd_for_fall_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_wtd_for_fall_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtd_for_fall.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtd_for_wtr_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[5]][[7]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "percent forest", dhat = coy_wtd_for_wtr_out, y_up = 0.6, season = "Winter")
  (coy_wtd_for_wtr_overlap_plot <- coy_wtd_for_wtr_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coy_wtd_for_wtr_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_wtd_for_wtr_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtd_for_wtr.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtd_for_sprg_overPlot_g <- overlap_pred_prey_plots(pred_prey_for_overlap[[5]][[8]], name1 = "Coyote", name2 = "White-tailed deer", name3 = "percent forest", dhat = coy_wtd_for_sprg_out, y_up = 0.6, season = "Spring")
  (coy_wtd_for_sprg_overlap_plot <- coy_wtd_for_sprg_overPlot_g[[1]] + theme(legend.position = c(0.40, 0.895)) + coy_wtd_for_sprg_overPlot_g[[2]] + theme(legend.position = c(0.40, 0.895)) & theme(text = element_text(size = 14)))
  ggsave(coy_wtd_for_sprg_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtd_for_sprg.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  # library(patchwork)
  # tst <- coug_md_tri_fall_overlap_plot + bear_md_tri_fall_overlap_plot + bob_md_tri_wtr_overlap_plot +
  #   coy_md_tri_smr_overlap_plot + coy_md_tri_wtr_overlap_plot + coug_md_for_smr_overlap_plot + 
  #   coug_md_for_fall_overlap_plot + coug_md_for_wtr_overlap_plot + wolf_md_for_fall_overlap_plot +
  #   bob_md_for_smr_overlap_plot + coy_md_for_wtr_overlap_plot + coy_md_for_sprg_overlap_plot
  # md_pred_patchwork1 <- coy_md_tri_smr_overlap_plot + bob_md_for_smr_overlap_plot + 
  #   coug_md_for_smr_overlap_plot + plot_layout(ncol = 1, widths = 1)
  # md_pred_patchwork2 <- bear_md_tri_fall_overlap_plot + coug_md_tri_fall_overlap_plot + 
  #   wolf_md_for_fall_overlap_plot + plot_layout(ncol = 1, widths = 1)
  # md_pred_patchwork3 <- bob_md_tri_wtr_overlap_plot + coy_md_tri_wtr_overlap_plot +
  #   coy_md_for_wtr_overlap_plot +plot_layout(ncol = 1, widths = 1)
  # ggsave("./Outputs/patchwork1.png", md_pred_patchwork1, width = 9, height = 15, units = "in")
  # ggsave("./Outputs/patchwork2.png", md_pred_patchwork2, width = 9, height = 15, units = "in")
  # ggsave("./Outputs/patchwork3.png", md_pred_patchwork3, width = 9, height = 15, units = "in")
  # 
  # 
  # md_pred_patchwork <- coy_md_tri_smr_overlap_plot + bear_md_tri_fall_overlap_plot + bob_md_tri_wtr_overlap_plot + 
  #                      bob_md_for_smr_overlap_plot + coug_md_tri_fall_overlap_plot + coy_md_tri_wtr_overlap_plot + 
  #                      coug_md_for_smr_overlap_plot + wolf_md_for_fall_overlap_plot + coy_md_for_wtr_overlap_plot +
  #   plot_layout(ncol = 3, widths = 1)
  # md_pred_patchwork
  # ggsave("./Outputs/patchwork.png", md_pred_patchwork, height = 15, units = "in")
  

  #'  ---------------------------------
  ####  Prey Activity Overlap Curves  ####
  #'  ---------------------------------
  #'  Plot activity curves for each prey species at cameras with low versus high 
  #'  levels of background risk
  overlap_singlespp_plots <- function(dat, name2, name3, dhat, y_up, season, x_start) {
    #'  Sample sizes for prey when background risk is low and high
    n2low <- dat$ndet_spp2.lowrisk; n2high <- dat$ndet_spp2.highrisk
    spp2low <- paste0("Low ", name3, " (n = ", n2low, ")"); spp2high <- paste0("High ", name3, " (n = ", n2high, ")")
    #'  Temporal overlap between predators and prey where background risk is low or high
    dhatmu <- dhat[1,4]; dhatl <- dhat[1,5]; dhatu<- dhat[1,6]
    #'  Density data for overlap plots
    overdensity <- dat[[6]]
    #'  Separate data sets based on whether background risk is low or high
    low <- overdensity[overdensity$BackgroundRisk == "Low",] %>% rename("densityA" = "density")
    high <- overdensity[overdensity$BackgroundRisk == "High",] %>% rename("densityB" = "density")
    overdensity <- full_join(low, high, by = c("x", "Species", "RiskType"))
    
    overlap <- ggplot(overdensity, aes(x, densityA, colour = BackgroundRisk.x)) +
      geom_line(lwd = 0.75) + 
      geom_line(aes(x, densityB, colour = BackgroundRisk.y), lwd = 0.75) +  
      geom_area(aes(y = pmin(densityA, densityB)),
                alpha = 0.3, color = NA) +
      scale_x_continuous(breaks = c(0, 1.57, 3.0, 4.71, 6.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = pi/2, linetype="dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype="dotted") +
      theme_bw() +
      theme(legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      ylim(0, y_up) +
      theme(legend.position = c(x_start, 0.895)) +
      labs(x = "Time of day", y = "Density", color = paste0("\u0394 = ", dhatmu, " (", dhatl, " - ", dhatu, ")"), title = paste0(name2, " temporal overlap, ", season)) +
      scale_colour_manual(breaks = c("Low", "High"), labels = c(spp2low, spp2high), values = c("red", "blue")) +
      theme(text = element_text(size = 14), legend.text = element_text(size = 14))
    plot(overlap)
    
    return(overlap)
  }
  ####  Mule deer summer  ####
  md_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[1]], name2 = "Mule deer", name3 = "TRI", dhat = md_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(md_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[2]], name2 = "Mule deer", name3 = "% Forest", dhat = md_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(md_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[3]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_smr_coug_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(md_smr_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_smr_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[4]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_smr_wolf_out, y_up = 0.6, season = "Summer", x_start = 0.25)
  ggsave(md_smr_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_smr_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[5]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_smr_bear_out, y_up = 0.6, season = "Summer", x_start = 0.30)
  ggsave(md_smr_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_smr_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_smr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[6]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_smr_bob_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(md_smr_bob_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_smr_bob.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_smr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[7]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_smr_coy_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(md_smr_coy_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_smr_coy.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Mule deer fall  ####
  md_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[8]], name2 = "Mule deer", name3 = "TRI", dhat = md_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(md_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[9]], name2 = "Mule deer", name3 = "% Forest", dhat = md_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(md_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[10]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_fall_coug_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(md_fall_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_fall_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_fall_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[11]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_fall_wolf_out, y_up = 0.6, season = "Fall", x_start = 0.25)
  ggsave(md_fall_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_fall_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[12]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_fall_bear_out, y_up = 0.6, season = "Fall", x_start = 0.30)
  ggsave(md_fall_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_fall_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_fall_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[13]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_fall_bob_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(md_fall_bob_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_fall_bob.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_fall_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[14]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_fall_coy_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(md_fall_coy_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_fall_coy.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Mule deer winter  ####
  md_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[15]], name2 = "Mule deer", name3 = "TRI", dhat = md_wtr_tri_out, y_up = 0.6, season = "Winter", x_start = 0.20)
  ggsave(md_wtr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_wtr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[16]], name2 = "Mule deer", name3 = "% Forest", dhat = md_wtr_for_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(md_wtr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_wtr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[17]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_wtr_coug_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(md_wtr_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_wtr_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[18]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_wtr_wolf_out, y_up = 0.6, season = "Winter", x_start = 0.25)
  ggsave(md_wtr_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_wtr_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_wtr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[19]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_wtr_bob_out, y_up = 0.6, season = "Winter", x_start = 0.38)
  ggsave(md_wtr_bob_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_wtr_bob.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_wtr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[20]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_wtr_coy_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(md_wtr_coy_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_wtr_coy.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Mule deer spring  ####
  md_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[21]], name2 = "Mule deer", name3 = "TRI", dhat = md_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(md_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[22]], name2 = "Mule deer", name3 = "% Forest", dhat = md_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(md_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[23]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_sprg_coug_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(md_sprg_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_sprg_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_sprg_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[24]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_sprg_wolf_out, y_up = 0.6, season = "Spring", x_start = 0.25)
  ggsave(md_sprg_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_sprg_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[25]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_sprg_bear_out, y_up = 0.6, season = "Spring", x_start = 0.30)
  ggsave(md_sprg_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_sprg_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_sprg_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[26]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_sprg_bob_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(md_sprg_bob_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_sprg_bob.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  md_sprg_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[27]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_sprg_coy_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(md_sprg_coy_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_md_sprg_coy.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ####  Elk summer  ####
  elk_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[1]], name2 = "Elk", name3 = "TRI", dhat = elk_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(elk_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[2]], name2 = "Elk", name3 = "% Forest", dhat = elk_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(elk_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[3]], name2 = "Elk", name3 = "cougar risk", dhat = elk_smr_coug_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(elk_smr_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_smr_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[4]], name2 = "Elk", name3 = "wolf risk", dhat = elk_smr_wolf_out, y_up = 0.6, season = "Summer", x_start = 0.25)
  ggsave(elk_smr_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_smr_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[5]], name2 = "Elk", name3 = "black bear risk", dhat = elk_smr_bear_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(elk_smr_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_smr_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Elk fall  ####
  elk_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[6]], name2 = "Elk", name3 = "TRI", dhat = elk_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(elk_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[7]], name2 = "Elk", name3 = "% Forest", dhat = elk_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(elk_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[8]], name2 = "Elk", name3 = "cougar risk", dhat = elk_fall_coug_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(elk_fall_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_fall_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[9]], name2 = "Elk", name3 = "black bear risk", dhat = elk_fall_bear_out, y_up = 0.6, season = "Fall", x_start = 0.30)
  ggsave(elk_fall_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_fall_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Elk winter  ####
  elk_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[10]], name2 = "Elk", name3 = "TRI", dhat = elk_wtr_tri_out, y_up = 0.6, season = "Winter", x_start = 0.20)
  ggsave(elk_wtr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_wtr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[11]], name2 = "Elk", name3 = "% Forest", dhat = elk_wtr_for_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(elk_wtr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_wtr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[12]], name2 = "Elk", name3 = "cougar risk", dhat = elk_wtr_coug_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(elk_wtr_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_wtr_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Elk spring  ####
  elk_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[13]], name2 = "Elk", name3 = "TRI", dhat = elk_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(elk_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[14]], name2 = "Elk", name3 = "% Forest", dhat = elk_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(elk_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[15]], name2 = "Elk", name3 = "cougar risk", dhat = elk_sprg_coug_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(elk_sprg_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_sprg_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  elk_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[16]], name2 = "Elk", name3 = "black bear risk", dhat = elk_sprg_bear_out, y_up = 0.6, season = "Spring", x_start = 0.30)
  ggsave(elk_sprg_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_elk_sprg_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ####  Moose summer  ####
  moose_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[1]], name2 = "Moose", name3 = "TRI", dhat = moose_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(moose_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[2]], name2 = "Moose", name3 = "% Forest", dhat = moose_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(moose_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[3]], name2 = "Moose", name3 = "cougar risk", dhat = moose_smr_coug_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(moose_smr_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_smr_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[4]], name2 = "Moose", name3 = "wolf risk", dhat = moose_smr_wolf_out, y_up = 0.6, season = "Summer", x_start = 0.25)
  ggsave(moose_smr_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_smr_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[5]], name2 = "Moose", name3 = "black bear risk", dhat = moose_smr_bear_out, y_up = 0.6, season = "Summer", x_start = 0.30)
  ggsave(moose_smr_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_smr_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Moose fall  ####
  moose_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[6]], name2 = "Moose", name3 = "TRI", dhat = moose_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(moose_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[7]], name2 = "Moose", name3 = "% Forest", dhat = moose_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(moose_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[8]], name2 = "Moose", name3 = "cougar risk", dhat = moose_fall_coug_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(moose_fall_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_fall_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_fall_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[9]], name2 = "Moose", name3 = "wolf risk", dhat = moose_fall_wolf_out, y_up = 0.6, season = "Fall", x_start = 0.25)
  ggsave(moose_fall_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_fall_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[10]], name2 = "Moose", name3 = "black bear risk", dhat = moose_fall_bear_out, y_up = 0.6, season = "Fall", x_start = 0.30)
  ggsave(moose_fall_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_fall_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Moose winter  ####
  moose_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[11]], name2 = "Moose", name3 = "TRI", dhat = moose_wtr_tri_out, y_up = 0.6, season = "Winter", x_start = 0.20)
  ggsave(moose_wtr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_wtr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[12]], name2 = "Moose", name3 = "% Forest", dhat = moose_wtr_for_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(moose_wtr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_wtr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[13]], name2 = "Moose", name3 = "cougar risk", dhat = moose_wtr_coug_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(moose_wtr_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_wtr_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[14]], name2 = "Moose", name3 = "wolf risk", dhat = moose_wtr_wolf_out, y_up = 0.6, season = "Winter", x_start = 0.25)
  ggsave(moose_wtr_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_wtr_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  Moose spring  ####
  moose_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[15]], name2 = "Moose", name3 = "TRI", dhat = moose_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(moose_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[16]], name2 = "Moose", name3 = "% Forest", dhat = moose_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(moose_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[17]], name2 = "Moose", name3 = "cougar risk", dhat = moose_sprg_coug_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(moose_sprg_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_sprg_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_sprg_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[18]], name2 = "Moose", name3 = "wolf risk", dhat = moose_sprg_wolf_out, y_up = 0.6, season = "Spring", x_start = 0.25)
  ggsave(moose_sprg_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_sprg_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  moose_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[19]], name2 = "Moose", name3 = "black bear risk", dhat = moose_sprg_bear_out, y_up = 0.6, season = "Spring", x_start = 0.30)
  ggsave(moose_sprg_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_moose_sprg_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ####  White-tailed deer summer  ####
  wtd_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[1]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(wtd_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[2]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(wtd_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[3]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_smr_coug_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(wtd_smr_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_smr_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[4]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_smr_wolf_out, y_up = 0.6, season = "Summer", x_start = 0.25)
  ggsave(wtd_smr_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_smr_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[5]], name2 = "White-tailed deer", name3 = "black bear risk", dhat = wtd_smr_bear_out, y_up = 0.6, season = "Summer", x_start = 0.30)
  ggsave(wtd_smr_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_smr_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_smr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[6]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_smr_bob_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(wtd_smr_bob_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_smr_bob.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_smr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[7]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_smr_coy_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(wtd_smr_coy_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_smr_coy.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  White-tailed deer fall  ####
  wtd_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[8]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(wtd_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[9]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(wtd_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[10]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_fall_coug_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(wtd_fall_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_fall_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_fall_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[11]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_fall_wolf_out, y_up = 0.6, season = "Fall", x_start = 0.25)
  ggsave(wtd_fall_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_fall_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[12]], name2 = "White-tailed deer", name3 = "black bear risk", dhat = wtd_fall_bear_out, y_up = 0.6, season = "Fall", x_start = 0.30)
  ggsave(wtd_fall_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_fall_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_fall_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[13]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_fall_bob_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(wtd_fall_bob_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_fall_bob.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_fall_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[14]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_fall_coy_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(wtd_fall_coy_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_fall_coy.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  White-tailed deer winter  ####
  wtd_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[15]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_wtr_tri_out, y_up = 0.6, season = "Winter", x_start = 0.20)
  ggsave(wtd_wtr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_wtr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[16]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_wtr_for_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(wtd_wtr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_wtr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[17]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_wtr_coug_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(wtd_wtr_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_wtr_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[18]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_wtr_wolf_out, y_up = 0.6, season = "Winter", x_start = 0.25)
  ggsave(wtd_wtr_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_wtr_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_wtr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[19]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_wtr_bob_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(wtd_wtr_bob_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_wtr_bob.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_wtr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[20]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_wtr_coy_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(wtd_wtr_coy_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_wtr_coy.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  ####  White-tailed deer spring  ####
  wtd_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[21]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(wtd_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[22]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(wtd_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[23]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_sprg_coug_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(wtd_sprg_coug_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_sprg_coug.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_sprg_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[24]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_sprg_wolf_out, y_up = 0.6, season = "Spring", x_start = 0.25)
  ggsave(wtd_sprg_wolf_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_sprg_wolf.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[25]], name2 = "White-tailed deer", name3 = "black bear risk", dhat = wtd_sprg_bear_out, y_up = 0.6, season = "Spring", x_start = 0.30)
  ggsave(wtd_sprg_bear_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_sprg_bear.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_sprg_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[26]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_sprg_bob_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(wtd_sprg_bob_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_sprg_bob.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wtd_sprg_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[27]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_sprg_coy_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(wtd_sprg_coy_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wtd_sprg_coy.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  ####  Black bear  ####
  bear_smr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[1]][[1]], name2 = "Black bear", name3 = "TRI", dhat = bear_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(bear_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bear_smr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[1]][[2]], name2 = "Black bear", name3 = "% Forest", dhat = bear_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(bear_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bear_fall_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[1]][[3]], name2 = "Black bear", name3 = "TRI", dhat = bear_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(bear_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bear_fall_for_overPlot <- overlap_singlespp_plots(pred_overlap[[1]][[4]], name2 = "Black bear", name3 = "% Forest", dhat = bear_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(bear_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bear_sprg_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[1]][[5]], name2 = "Black bear", name3 = "TRI", dhat = bear_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(bear_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bear_sprg_for_overPlot <- overlap_singlespp_plots(pred_overlap[[1]][[6]], name2 = "Black bear", name3 = "% Forest", dhat = bear_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(bear_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bear_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  

  ####  Bobcat  ####
  bob_smr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[2]][[1]], name2 = "Bobcat", name3 = "TRI", dhat = bob_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(bob_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_smr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[2]][[2]], name2 = "Bobcat", name3 = "% Forest", dhat = bob_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(bob_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_fall_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[2]][[3]], name2 = "Bobcat", name3 = "TRI", dhat = bob_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(bob_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_fall_for_overPlot <- overlap_singlespp_plots(pred_overlap[[2]][[4]], name2 = "Bobcat", name3 = "% Forest", dhat = bob_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(bob_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_wtr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[2]][[5]], name2 = "Bobcat", name3 = "TRI", dhat = bob_wtr_tri_out, y_up = 0.6, season = "Winter", x_start = 0.20)
  ggsave(bob_wtr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_wtr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[2]][[6]], name2 = "Bobcat", name3 = "% Forest", dhat = bob_wtr_for_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(bob_wtr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_wtr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_sprg_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[2]][[7]], name2 = "Bobcat", name3 = "TRI", dhat = bob_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(bob_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  bob_sprg_for_overPlot <- overlap_singlespp_plots(pred_overlap[[2]][[8]], name2 = "Bobcat", name3 = "% Forest", dhat = bob_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(bob_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_bob_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  

  ####  Cougar  ####
  coug_smr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[3]][[1]], name2 = "Cougar", name3 = "TRI", dhat = coug_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(coug_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coug_smr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[3]][[2]], name2 = "Cougar", name3 = "% Forest", dhat = coug_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(coug_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coug_fall_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[3]][[3]], name2 = "Cougar", name3 = "TRI", dhat = coug_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(coug_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coug_fall_for_overPlot <- overlap_singlespp_plots(pred_overlap[[3]][[4]], name2 = "Cougar", name3 = "% Forest", dhat = coug_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(coug_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coug_wtr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[3]][[5]], name2 = "Cougar", name3 = "TRI", dhat = coug_wtr_tri_out, y_up = 0.6, season = "Winter", x_start = 0.20)
  ggsave(coug_wtr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coug_wtr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[3]][[6]], name2 = "Cougar", name3 = "% Forest", dhat = coug_wtr_for_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(coug_wtr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_wtr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coug_sprg_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[3]][[7]], name2 = "Cougar", name3 = "TRI", dhat = coug_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(coug_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coug_sprg_for_overPlot <- overlap_singlespp_plots(pred_overlap[[3]][[8]], name2 = "Cougar", name3 = "% Forest", dhat = coug_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(coug_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coug_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')

  
  ####  Coyote  ####
  coy_smr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[4]][[1]], name2 = "Coyote", name3 = "TRI", dhat = coy_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(coy_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_smr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[4]][[2]], name2 = "Coyote", name3 = "% Forest", dhat = coy_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(coy_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_fall_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[4]][[3]], name2 = "Coyote", name3 = "TRI", dhat = coy_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(coy_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_fall_for_overPlot <- overlap_singlespp_plots(pred_overlap[[4]][[4]], name2 = "Coyote", name3 = "% Forest", dhat = coy_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(coy_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[4]][[5]], name2 = "Coyote", name3 = "TRI", dhat = coy_wtr_tri_out, y_up = 0.6, season = "Winter", x_start = 0.20)
  ggsave(coy_wtr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_wtr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[4]][[6]], name2 = "Coyote", name3 = "% Forest", dhat = coy_wtr_for_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(coy_wtr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_wtr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_sprg_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[4]][[7]], name2 = "Coyote", name3 = "TRI", dhat = coy_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(coy_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  coy_sprg_for_overPlot <- overlap_singlespp_plots(pred_overlap[[4]][[8]], name2 = "Coyote", name3 = "% Forest", dhat = coy_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(coy_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_coy_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')

  
  ####  Wolf  ####
  wolf_smr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[5]][[1]], name2 = "Wolf", name3 = "TRI", dhat = wolf_smr_tri_out, y_up = 0.6, season = "Summer", x_start = 0.20)
  ggsave(wolf_smr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_smr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_smr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[5]][[2]], name2 = "Wolf", name3 = "% Forest", dhat = wolf_smr_for_out, y_up = 0.6, season = "Summer", x_start = 0.28)
  ggsave(wolf_smr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_smr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_fall_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[5]][[3]], name2 = "Wolf", name3 = "TRI", dhat = wolf_fall_tri_out, y_up = 0.6, season = "Fall", x_start = 0.20)
  ggsave(wolf_fall_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_fall_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_fall_for_overPlot <- overlap_singlespp_plots(pred_overlap[[5]][[4]], name2 = "Wolf", name3 = "% Forest", dhat = wolf_fall_for_out, y_up = 0.6, season = "Fall", x_start = 0.28)
  ggsave(wolf_fall_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_fall_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_wtr_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[5]][[5]], name2 = "Wolf", name3 = "TRI", dhat = wolf_wtr_tri_out, y_up = 0.6, season = "Winter", x_start = 0.20)
  ggsave(wolf_wtr_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtr_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_wtr_for_overPlot <- overlap_singlespp_plots(pred_overlap[[5]][[6]], name2 = "Wolf", name3 = "% Forest", dhat = wolf_wtr_for_out, y_up = 0.6, season = "Winter", x_start = 0.28)
  ggsave(wolf_wtr_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_wtr_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_sprg_tri_overPlot <- overlap_singlespp_plots(pred_overlap[[5]][[7]], name2 = "Wolf", name3 = "TRI", dhat = wolf_sprg_tri_out, y_up = 0.6, season = "Spring", x_start = 0.20)
  ggsave(wolf_sprg_tri_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_sprg_tri.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  wolf_sprg_for_overPlot <- overlap_singlespp_plots(pred_overlap[[5]][[8]], name2 = "Wolf", name3 = "% Forest", dhat = wolf_sprg_for_out, y_up = 0.6, season = "Spring", x_start = 0.28)
  ggsave(wolf_sprg_for_overPlot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Overlap_Plot_wolf_sprg_for.tiff", width = 6, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  