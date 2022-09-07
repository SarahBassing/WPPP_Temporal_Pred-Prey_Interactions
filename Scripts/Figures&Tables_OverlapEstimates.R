  #'  =====================================================
  #'  Temporal overlap analysis: results tables and figures
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  August 2022
  #'  =====================================================
  #'  Script pulls in outputs from temporal overlap analyses conducted with the
  #'  Temporal_overlap.R script, then generates result tables & figures for publication.
  #'  =====================================================
  
  #'  Libraries
  library(ggplot2)
  library(khroma)
  library(patchwork)
  library(sp)
  library(raster)
  library(tidyverse)  

  #'  Load output from temporal overlap analysis
  load("./Outputs/Temporal Overlap/PredPrey_LowHi_Overlap_2022-08-25.RData")
  load("./Outputs/Temporal Overlap/PreyOnly_LowHi_Overlap_2022-08-25.RData")
  
  
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
    names(df)[names(df) == "risk"] <- paste0(risk_type, "_background_risk")
    df <- mutate(df, Dhat = as.numeric(Dhat),
                 l95 = as.numeric(l95),
                 u95 = as.numeric(u95))
    return(df)
  }
  #'  When using two different overlap estimators for same species pairing (e.g.,
  #'  when low vs high risk sample sizes requires different overlap estimators), 
  #'  need to filter to the appropriate estimator given sample size- 
  #'  dhat1 for when n < 50 observations, dhat4 for n > 50 observations
  #'  Note: overlap estimate when risk is low is always first row [1,],  
  #'  overlap estimate when risk is high is always second row [2,]
  ####  Cougar - Prey Output  ####
  #'  Cougar-mule deer
  coug_md_smr_out <- predprey_table(pred_prey_overlap[[1]][[1]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Complexity", season = "Summer")
  coug_md_fall_out <- predprey_table(pred_prey_overlap[[1]][[2]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Complexity", season = "Fall")
  coug_md_wtr_out <- predprey_table(pred_prey_overlap[[1]][[3]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Complexity", season = "Winter")
  #'  Low risk >50 cougars = dhat4; High risk <50 cougars = dhat1
  coug_md_sprg_out1 <- predprey_table(pred_prey_overlap[[1]][[4]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Complexity", season = "Spring")
  coug_md_sprg_out4 <- predprey_table(pred_prey_overlap[[1]][[5]], spp1 = "Cougar", spp2 = "Mule Deer", risk_type = "Complexity", season = "Spring")
  coug_md_sprg_out <- rbind(coug_md_sprg_out4[1,], coug_md_sprg_out1[2,])
  #'  Cougar-elk
  coug_elk_smr_out <- predprey_table(pred_prey_overlap[[1]][[6]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Complexity", season = "Summer")
  coug_elk_fall_out <- predprey_table(pred_prey_overlap[[1]][[7]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Complexity", season = "Fall")
  coug_elk_wtr_out <- predprey_table(pred_prey_overlap[[1]][[8]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Complexity", season = "Winter")
  coug_elk_sprg_out <- predprey_table(pred_prey_overlap[[1]][[9]], spp1 = "Cougar", spp2 = "Elk", risk_type = "Complexity", season = "Spring")
  #'  Cougar-moose
  coug_moose_smr_out <- predprey_table(pred_prey_overlap[[1]][[10]], spp1 = "Cougar", spp2 = "Moose", risk_type = "Complexity", season = "Summer")
  coug_moose_fall_out <- predprey_table(pred_prey_overlap[[1]][[11]], spp1 = "Cougar", spp2 = "Moose", risk_type = "Complexity", season = "Fall")
  coug_moose_wtr_out <- predprey_table(pred_prey_overlap[[1]][[12]], spp1 = "Cougar", spp2 = "Moose", risk_type = "Complexity", season = "Winter")
  coug_moose_sprg_out <- predprey_table(pred_prey_overlap[[1]][[13]], spp1 = "Cougar", spp2 = "Moose", risk_type = "Complexity", season = "Spring")
  #'  Cougar-white-tailed deer
  coug_wtd_smr_out <- predprey_table(pred_prey_overlap[[1]][[14]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Summer")
  #'  Low risk <50 cougars; High risk >50 cougars
  coug_wtd_fall_out1 <- predprey_table(pred_prey_overlap[[1]][[15]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Fall")
  coug_wtd_fall_out4 <- predprey_table(pred_prey_overlap[[1]][[16]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Fall")
  coug_wtd_fall_out <- rbind(coug_wtd_fall_out1[1,], coug_wtd_fall_out4[2,])
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_wtd_wtr_out1 <- predprey_table(pred_prey_overlap[[1]][[17]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Winter")
  coug_wtd_wtr_out4 <- predprey_table(pred_prey_overlap[[1]][[18]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Winter")
  coug_wtd_wtr_out <- rbind(coug_wtd_wtr_out4[1,], coug_wtd_wtr_out1[2,])
  coug_wtd_sprg_out <- predprey_table(pred_prey_overlap[[1]][[19]], spp1 = "Cougar", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Spring")
  #'  Cougar-prey results table
  coug_prey_out <- rbind(coug_md_smr_out, coug_md_fall_out, coug_md_wtr_out, coug_md_sprg_out,
                         coug_elk_smr_out, coug_elk_fall_out, coug_elk_wtr_out, coug_elk_sprg_out,
                         coug_moose_smr_out, coug_moose_fall_out, coug_moose_wtr_out, coug_moose_sprg_out,
                         coug_wtd_smr_out, coug_wtd_fall_out, coug_wtd_wtr_out, coug_wtd_sprg_out)
  
  ####  Wolf - Prey Output  ####
  #'  Wolf-mule deer
  wolf_md_smr_out <- predprey_table(pred_prey_overlap[[2]][[1]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "Complexity", season = "Summer")
  wolf_md_fall_out <- predprey_table(pred_prey_overlap[[2]][[2]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "Complexity", season = "Fall")
  wolf_md_sprg_out <- predprey_table(pred_prey_overlap[[2]][[3]], spp1 = "Wolf", spp2 = "Mule Deer", risk_type = "Complexity", season = "Spring")
  #'  Wolf-elk
  wolf_elk_smr_out <- predprey_table(pred_prey_overlap[[2]][[4]], spp1 = "Wolf", spp2 = "Elk", risk_type = "Complexity", season = "Summer")
  #'  Wolf-moose
  wolf_moose_smr_out <- predprey_table(pred_prey_overlap[[2]][[5]], spp1 = "Wolf", spp2 = "Moose", risk_type = "Complexity", season = "Summer")
  wolf_moose_fall_out <- predprey_table(pred_prey_overlap[[2]][[6]], spp1 = "Wolf", spp2 = "Moose", risk_type = "Complexity", season = "Fall")
  wolf_moose_wtr_out <- predprey_table(pred_prey_overlap[[2]][[7]], spp1 = "Wolf", spp2 = "Moose", risk_type = "Complexity", season = "Winter")
  #'  Wolf-white-tailed deer
  wolf_wtd_smr_out <- predprey_table(pred_prey_overlap[[2]][[8]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Summer")
  wolf_wtd_fall_out <- predprey_table(pred_prey_overlap[[2]][[9]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Fall")
  wolf_wtd_wtr_out <- predprey_table(pred_prey_overlap[[2]][[10]], spp1 = "Wolf", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Winter")
  #'  Wolf-prey results table
  wolf_prey_out <- rbind(wolf_md_smr_out, wolf_md_fall_out, wolf_md_sprg_out, wolf_elk_smr_out,
                         wolf_moose_smr_out, wolf_moose_fall_out, wolf_moose_wtr_out,
                         wolf_wtd_smr_out, wolf_wtd_fall_out, wolf_wtd_wtr_out)
  
  ####  Black bear - Prey Output  ####
  #'  Black bear-mule deer
  bear_md_smr_out <- predprey_table(pred_prey_overlap[[3]][[1]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "Complexity", season = "Summer")
  #'  Low risk >50 black bears; High risk <50 black bears
  bear_md_fall_out1 <- predprey_table(pred_prey_overlap[[3]][[2]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "Complexity", season = "Fall")
  bear_md_fall_out4 <- predprey_table(pred_prey_overlap[[3]][[3]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "Complexity", season = "Fall")
  bear_md_fall_out <- rbind(bear_md_fall_out4[1,], bear_md_fall_out1[2,])
  bear_md_sprg_out <- predprey_table(pred_prey_overlap[[3]][[4]], spp1 = "Black bear", spp2 = "Mule Deer", risk_type = "Complexity", season = "Spring")
  #'  Black bear-elk
  bear_elk_smr_out <- predprey_table(pred_prey_overlap[[3]][[5]], spp1 = "Black bear", spp2 = "Elk", risk_type = "Complexity", season = "Summer")
  bear_elk_fall_out <- predprey_table(pred_prey_overlap[[3]][[6]], spp1 = "Black bear", spp2 = "Elk", risk_type = "Complexity", season = "Fall")
  bear_elk_sprg_out <- predprey_table(pred_prey_overlap[[3]][[7]], spp1 = "Black bear", spp2 = "Elk", risk_type = "Complexity", season = "Spring")
  #'  Black bear-moose
  bear_moose_smr_out <- predprey_table(pred_prey_overlap[[3]][[8]], spp1 = "Black bear", spp2 = "Moose", risk_type = "Complexity", season = "Summer")
  bear_moose_fall_out <- predprey_table(pred_prey_overlap[[3]][[9]], spp1 = "Black bear", spp2 = "Moose", risk_type = "Complexity", season = "Fall")
  bear_moose_sprg_out <- predprey_table(pred_prey_overlap[[3]][[10]], spp1 = "Black bear", spp2 = "Moose", risk_type = "Complexity", season = "Spring")
  #'  Black bear-white-tailed deer
  bear_wtd_smr_out <- predprey_table(pred_prey_overlap[[3]][[11]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Summer")
  bear_wtd_fall_out <- predprey_table(pred_prey_overlap[[3]][[12]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Fall")
  #'  Low risk >50 black bears; High risk <50 black bears
  bear_wtd_sprg_out1 <- predprey_table(pred_prey_overlap[[3]][[13]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Spring")
  bear_wtd_sprg_out4 <- predprey_table(pred_prey_overlap[[3]][[14]], spp1 = "Black bear", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Spring")
  bear_wtd_sprg_out <- rbind(bear_wtd_sprg_out4[1,], bear_wtd_sprg_out1[2,])
  #'  Black bear-prey results table
  bear_prey_out <- rbind(bear_md_smr_out, bear_md_fall_out, bear_md_sprg_out,
                         bear_elk_smr_out, bear_elk_fall_out, bear_elk_sprg_out,
                         bear_moose_smr_out, bear_moose_fall_out, bear_moose_sprg_out,
                         bear_wtd_smr_out, bear_wtd_fall_out, bear_wtd_sprg_out)
  
  ####  Bobcat - Prey Overlap  ####
  #'  Bobcat-mule deer
  bob_md_smr_out <- predprey_table(pred_prey_overlap[[4]][[1]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "Complexity", season = "Summer")
  bob_md_fall_out <- predprey_table(pred_prey_overlap[[4]][[2]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "Complexity", season = "Fall")
  bob_md_sprg_out <- predprey_table(pred_prey_overlap[[4]][[3]], spp1 = "Bobcat", spp2 = "Mule Deer", risk_type = "Complexity", season = "Spring")
  #'  Bobcat-white-tailed deer
  bob_wtd_smr_out <- predprey_table(pred_prey_overlap[[4]][[4]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Summer")
  #'  Low risk <50 bobcat; High risk >50 bobcat
  bob_wtd_fall_out1 <- predprey_table(pred_prey_overlap[[4]][[5]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Fall")
  bob_wtd_fall_out4 <- predprey_table(pred_prey_overlap[[4]][[6]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Fall")
  bob_wtd_fall_out <- rbind(bob_wtd_fall_out1[1,], bob_wtd_fall_out4[2,])
  bob_wtd_wtr_out <- predprey_table(pred_prey_overlap[[4]][[7]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Winter")
  bob_wtd_sprg_out <- predprey_table(pred_prey_overlap[[4]][[8]], spp1 = "Bobcat", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Spring")
  #'  Bobcat-prey results table
  bob_prey_out <- rbind(bob_md_smr_out, bob_md_fall_out, bob_md_sprg_out,
                        bob_wtd_smr_out, bob_wtd_fall_out, bob_wtd_wtr_out, bob_wtd_sprg_out)
  
  ####  Coyote - Prey Overlap  ####
  #'  Coyote-mule deer
  coy_md_smr_out <- predprey_table(pred_prey_overlap[[5]][[1]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "Complexity", season = "Summer")
  coy_md_fall_out <- predprey_table(pred_prey_overlap[[5]][[2]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "Complexity", season = "Fall")
  coy_md_wtr_out <- predprey_table(pred_prey_overlap[[5]][[3]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "Complexity", season = "Winter")
  coy_md_sprg_out <- predprey_table(pred_prey_overlap[[5]][[4]], spp1 = "Coyote", spp2 = "Mule Deer", risk_type = "Complexity", season = "Spring")
  #'  Coyote-white-tailed deer
  coy_wtd_smr_out <- predprey_table(pred_prey_overlap[[5]][[5]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Summer")
  coy_wtd_fall_out <- predprey_table(pred_prey_overlap[[5]][[6]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Fall")
  coy_wtd_wtr_out <- predprey_table(pred_prey_overlap[[5]][[7]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Winter")
  coy_wtd_sprg_out <- predprey_table(pred_prey_overlap[[5]][[8]], spp1 = "Coyote", spp2 = "White-tailed Deer", risk_type = "Complexity", season = "Spring")
  #'  Coyote-prey results table
  coy_prey_out <- rbind(coy_md_smr_out, coy_md_fall_out, coy_md_wtr_out, coy_md_sprg_out,
                        coy_wtd_smr_out, coy_wtd_fall_out, coy_wtd_wtr_out, coy_wtd_sprg_out)
  
  pred_prey_out_tbl <- rbind(coug_prey_out, wolf_prey_out, bear_prey_out, bob_prey_out, coy_prey_out)
  # write.csv(pred_prey_out_tbl, file = paste0("./Outputs/Temporal Overlap/pred-prey_overlap_tbl_", Sys.Date(), ".csv"))

  #'  -------------------------------
  #'  Species-specific results tables
  #'  -------------------------------
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
    df <- mutate(df, Dhat = as.numeric(Dhat),
                 l95 = as.numeric(l95),
                 u95 = as.numeric(u95))
    return(df)
  }
  #'  Coefficient of overlap for each species at low vs high risk cameras
  ####  Mule deer Ouptut  ####
  #'  Summer
  md_smr_hab_out <- prey_table(prey_overlap[[1]][[1]], spp2 = "Mule deer", season = "Summer", risk_type = "Habitat complexity")
  md_smr_coug_out <- prey_table(prey_overlap[[1]][[2]], spp2 = "Mule deer", season = "Summer", risk_type = "Cougar detected")
  md_smr_wolf_out <- prey_table(prey_overlap[[1]][[3]], spp2 = "Mule deer", season = "Summer", risk_type = "Wolf detected")
  md_smr_bear_out <- prey_table(prey_overlap[[1]][[4]], spp2 = "Mule deer", season = "Summer", risk_type = "Black bear detected")
  md_smr_bob_out <- prey_table(prey_overlap[[1]][[5]], spp2 = "Mule deer", season = "Summer", risk_type = "Bobcat detected")
  md_smr_coy_out <- prey_table(prey_overlap[[1]][[6]], spp2 = "Mule deer", season = "Summer", risk_type = "Coyote detected")
  #'  Fall
  md_fall_hab_out <- prey_table(prey_overlap[[1]][[7]], spp2 = "Mule deer", season = "Fall", risk_type = "Habitat complexity")
  md_fall_coug_out <- prey_table(prey_overlap[[1]][[8]], spp2 = "Mule deer", season = "Fall", risk_type = "Cougar detected")
  md_fall_wolf_out <- prey_table(prey_overlap[[1]][[9]], spp2 = "Mule deer", season = "Fall", risk_type = "Wolf detected")
  md_fall_bear_out <- prey_table(prey_overlap[[1]][[10]], spp2 = "Mule deer", season = "Fall", risk_type = "Black bear detected")
  md_fall_bob_out <- prey_table(prey_overlap[[1]][[11]], spp2 = "Mule deer", season = "Fall", risk_type = "Bobcat detected")
  md_fall_coy_out <- prey_table(prey_overlap[[1]][[12]], spp2 = "Mule deer", season = "Fall", risk_type = "Coyote detected")
  #'  Winter
  md_wtr_hab_out <- prey_table(prey_overlap[[1]][[13]], spp2 = "Mule deer", season = "Winter", risk_type = "Habitat complexity")
  md_wtr_coug_out <- prey_table(prey_overlap[[1]][[14]], spp2 = "Mule deer", season = "Winter", risk_type = "Cougar detected")
  md_wtr_wolf_out <- prey_table(prey_overlap[[1]][[15]], spp2 = "Mule deer", season = "Winter", risk_type = "Wolf detected")
  md_wtr_bob_out <- prey_table(prey_overlap[[1]][[16]], spp2 = "Mule deer", season = "Winter", risk_type = "Bobcat detected")
  md_wtr_coy_out <- prey_table(prey_overlap[[1]][[17]], spp2 = "Mule deer", season = "Winter", risk_type = "Coyote detected")
  #'  Spring
  md_sprg_hab_out <- prey_table(prey_overlap[[1]][[18]], spp2 = "Mule deer", season = "Spring", risk_type = "Habitat complexity")
  md_sprg_coug_out <- prey_table(prey_overlap[[1]][[19]], spp2 = "Mule deer", season = "Spring", risk_type = "Cougar detected")
  md_sprg_wolf_out <- prey_table(prey_overlap[[1]][[20]], spp2 = "Mule deer", season = "Spring", risk_type = "Wolf detected")
  md_sprg_bear_out <- prey_table(prey_overlap[[1]][[21]], spp2 = "Mule deer", season = "Spring", risk_type = "Black bear detected")
  md_sprg_bob_out <- prey_table(prey_overlap[[1]][[22]], spp2 = "Mule deer", season = "Spring", risk_type = "Bobcat detected")
  md_sprg_coy_out <- prey_table(prey_overlap[[1]][[23]], spp2 = "Mule deer", season = "Spring", risk_type = "Coyote detected")
  #'  Save mule deer overlap results
  md_overlap <- rbind(md_smr_hab_out, md_smr_coug_out, md_smr_wolf_out, md_smr_bear_out, md_smr_bob_out, md_smr_coy_out,
                      md_fall_hab_out, md_fall_coug_out, md_fall_wolf_out, md_fall_bear_out, md_fall_bob_out, md_fall_coy_out,
                      md_wtr_hab_out, md_wtr_coug_out, md_wtr_wolf_out, md_wtr_bob_out, md_wtr_coy_out,
                      md_sprg_hab_out, md_sprg_coug_out, md_sprg_wolf_out, md_sprg_bear_out, md_sprg_bob_out, md_sprg_coy_out)
  
  ####  Elk Ouptut  ####
  #'  Summer
  elk_smr_hab_out <- prey_table(prey_overlap[[2]][[1]], spp2 = "Elk", season = "Summer", risk_type = "Habitat complexity")
  elk_smr_coug_out <- prey_table(prey_overlap[[2]][[2]], spp2 = "Elk", season = "Summer", risk_type = "Cougar detected")
  elk_smr_wolf_out <- prey_table(prey_overlap[[2]][[3]], spp2 = "Elk", season = "Summer", risk_type = "Wolf detected")
  elk_smr_bear_out <- prey_table(prey_overlap[[2]][[4]], spp2 = "Elk", season = "Summer", risk_type = "Black bear detected")
  #'  Fall
  elk_fall_hab_out <- prey_table(prey_overlap[[2]][[5]], spp2 = "Elk", season = "Fall", risk_type = "Habitat complexity")
  elk_fall_coug_out <- prey_table(prey_overlap[[2]][[6]], spp2 = "Elk", season = "Fall", risk_type = "Cougar detected")
  elk_fall_wolf_out <- prey_table(prey_overlap[[2]][[7]], spp2 = "Elk", season = "Fall", risk_type = "Wolf detected")
  elk_fall_bear_out <- prey_table(prey_overlap[[2]][[8]], spp2 = "Elk", season = "Fall", risk_type = "Black bear detected")
  #'  Winter
  elk_wtr_hab_out <- prey_table(prey_overlap[[2]][[9]], spp2 = "Elk", season = "Winter", risk_type = "Habitat complexity")
  elk_wtr_coug_out <- prey_table(prey_overlap[[2]][[10]], spp2 = "Elk", season = "Winter", risk_type = "Cougar detected")
  #'  Spring
  elk_sprg_hab_out <- prey_table(prey_overlap[[2]][[11]], spp2 = "Elk", season = "Spring", risk_type = "Habitat complexity")
  elk_sprg_coug_out <- prey_table(prey_overlap[[2]][[12]], spp2 = "Elk", season = "Spring", risk_type = "Cougar detected")
  elk_sprg_wolf_out <- prey_table(prey_overlap[[2]][[13]], spp2 = "Elk", season = "Spring", risk_type = "Wolf detected")
  elk_sprg_bear_out <- prey_table(prey_overlap[[2]][[14]], spp2 = "Elk", season = "Spring", risk_type = "Black bear detected")
  #'  Save elk overlap results
  elk_overlap <- rbind(elk_smr_hab_out, elk_smr_coug_out, elk_smr_wolf_out, elk_smr_bear_out,
                       elk_fall_hab_out, elk_fall_coug_out, elk_fall_wolf_out, elk_fall_bear_out,
                       elk_wtr_hab_out, elk_wtr_coug_out, 
                       elk_sprg_hab_out, elk_sprg_coug_out, elk_sprg_wolf_out, elk_sprg_bear_out)
  
  ####  Moose Ouptut  ####
  #'  Summer
  moose_smr_hab_out <- prey_table(prey_overlap[[3]][[1]], spp2 = "Moose", season = "Summer", risk_type = "Habitat complexity")
  moose_smr_coug_out <- prey_table(prey_overlap[[3]][[2]], spp2 = "Moose", season = "Summer", risk_type = "Cougar detected")
  moose_smr_wolf_out <- prey_table(prey_overlap[[3]][[3]], spp2 = "Moose", season = "Summer", risk_type = "Wolf detected")
  moose_smr_bear_out <- prey_table(prey_overlap[[3]][[4]], spp2 = "Moose", season = "Summer", risk_type = "Black bear detected")
  #'  Fall
  moose_fall_hab_out <- prey_table(prey_overlap[[3]][[5]], spp2 = "Moose", season = "Fall", risk_type = "Habitat complexity")
  moose_fall_coug_out <- prey_table(prey_overlap[[3]][[6]], spp2 = "Moose", season = "Fall", risk_type = "Cougar detected")
  moose_fall_wolf_out <- prey_table(prey_overlap[[3]][[7]], spp2 = "Moose", season = "Fall", risk_type = "Wolf detected")
  moose_fall_bear_out <- prey_table(prey_overlap[[3]][[8]], spp2 = "Moose", season = "Fall", risk_type = "Black bear detected")
  #'  Winter
  moose_wtr_hab_out <- prey_table(prey_overlap[[3]][[9]], spp2 = "Moose", season = "Winter", risk_type = "Habitat complexity")
  moose_wtr_coug_out <- prey_table(prey_overlap[[3]][[10]], spp2 = "Moose", season = "Winter", risk_type = "Cougar detected")
  moose_wtr_wolf_out <- prey_table(prey_overlap[[3]][[11]], spp2 = "Moose", season = "Winter", risk_type = "Wolf detected")
  #'  Spring
  moose_sprg_hab_out <- prey_table(prey_overlap[[3]][[12]], spp2 = "Moose", season = "Spring", risk_type = "Habitat complexity")
  moose_sprg_coug_out <- prey_table(prey_overlap[[3]][[13]], spp2 = "Moose", season = "Spring", risk_type = "Cougar detected")
  moose_sprg_wolf_out <- prey_table(prey_overlap[[3]][[14]], spp2 = "Moose", season = "Spring", risk_type = "Wolf detected")
  moose_sprg_bear_out <- prey_table(prey_overlap[[3]][[15]], spp2 = "Moose", season = "Spring", risk_type = "Black bear detected")
  #'  Save elk overlap results
  moose_overlap <- rbind(moose_smr_hab_out, moose_smr_coug_out, moose_smr_wolf_out, moose_smr_bear_out,
                         moose_fall_hab_out, moose_fall_coug_out, moose_fall_wolf_out, moose_fall_bear_out,
                         moose_wtr_hab_out, moose_wtr_coug_out,  moose_wtr_wolf_out, 
                         moose_sprg_hab_out, moose_sprg_coug_out, moose_sprg_wolf_out, moose_sprg_bear_out)
  
  ####  White-tailed deer Ouptut  ####
  #'  Summer
  wtd_smr_hab_out <- prey_table(prey_overlap[[1]][[1]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Habitat complexity")
  wtd_smr_coug_out <- prey_table(prey_overlap[[1]][[2]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Cougar detected")
  wtd_smr_wolf_out <- prey_table(prey_overlap[[1]][[3]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Wolf detected")
  wtd_smr_bear_out <- prey_table(prey_overlap[[1]][[4]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Black bear detected")
  wtd_smr_bob_out <- prey_table(prey_overlap[[1]][[5]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Bobcat detected")
  wtd_smr_coy_out <- prey_table(prey_overlap[[1]][[6]], spp2 = "White-tailed deer", season = "Summer", risk_type = "Coyote detected")
  #'  Fall
  wtd_fall_hab_out <- prey_table(prey_overlap[[1]][[7]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Habitat complexity")
  wtd_fall_coug_out <- prey_table(prey_overlap[[1]][[8]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Cougar detected")
  wtd_fall_wolf_out <- prey_table(prey_overlap[[1]][[9]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Wolf detected")
  wtd_fall_bear_out <- prey_table(prey_overlap[[1]][[10]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Black bear detected")
  wtd_fall_bob_out <- prey_table(prey_overlap[[1]][[11]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Bobcat detected")
  wtd_fall_coy_out <- prey_table(prey_overlap[[1]][[12]], spp2 = "White-tailed deer", season = "Fall", risk_type = "Coyote detected")
  #'  Winter
  wtd_wtr_hab_out <- prey_table(prey_overlap[[1]][[13]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Habitat complexity")
  wtd_wtr_coug_out <- prey_table(prey_overlap[[1]][[14]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Cougar detected")
  wtd_wtr_wolf_out <- prey_table(prey_overlap[[1]][[15]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Wolf detected")
  wtd_wtr_bob_out <- prey_table(prey_overlap[[1]][[16]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Bobcat detected")
  wtd_wtr_coy_out <- prey_table(prey_overlap[[1]][[17]], spp2 = "White-tailed deer", season = "Winter", risk_type = "Coyote detected")
  #'  Spring
  wtd_sprg_hab_out <- prey_table(prey_overlap[[1]][[18]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Habitat complexity")
  wtd_sprg_coug_out <- prey_table(prey_overlap[[1]][[19]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Cougar detected")
  wtd_sprg_wolf_out <- prey_table(prey_overlap[[1]][[20]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Wolf detected")
  wtd_sprg_bear_out <- prey_table(prey_overlap[[1]][[21]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Black bear detected")
  wtd_sprg_bob_out <- prey_table(prey_overlap[[1]][[22]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Bobcat detected")
  wtd_sprg_coy_out <- prey_table(prey_overlap[[1]][[23]], spp2 = "White-tailed deer", season = "Spring", risk_type = "Coyote detected")
  #'  Save mule deer overlap results
  wtd_overlap <- rbind(wtd_smr_hab_out, wtd_smr_coug_out, wtd_smr_wolf_out, wtd_smr_bear_out, wtd_smr_bob_out, wtd_smr_coy_out,
                       wtd_fall_hab_out, wtd_fall_coug_out, wtd_fall_wolf_out, wtd_fall_bear_out, wtd_fall_bob_out, wtd_fall_coy_out,
                       wtd_wtr_hab_out, wtd_wtr_coug_out, wtd_wtr_wolf_out, wtd_wtr_bob_out, wtd_wtr_coy_out,
                       wtd_sprg_hab_out, wtd_sprg_coug_out, wtd_sprg_wolf_out, wtd_sprg_bear_out, wtd_sprg_bob_out, wtd_sprg_coy_out)
  
  prey_out_tbl <- rbind(md_overlap, elk_overlap, moose_overlap, wtd_overlap)
  # write.csv(prey_out_tbl, file = paste0("./Outputs/Temporal Overlap/prey_overlap_tbl_", Sys.Date(), ".csv"))
  
  #'  ---------------------------------
  ####  Coefficient Comparison Plots  ####
  #'  ---------------------------------
  #'  Compare coefficient of overlapping estimates for each species-pairing or
  #'  individual species at cameras with low versus high background risk
  
  #'  Load data in table format
  pred_prey_overlap_tbl <- read.csv("./Outputs/Temporal Overlap/pred-prey_overlap_tbl_2022-08-26.csv")
  prey_overlap_tbl <- read.csv("./Outputs/Temporal Overlap/prey_overlap_tbl_2022-08-26.csv")
  
  #'  Set seasonal factor levels
  pred_prey_overlap_tbl$season <- factor(pred_prey_overlap_tbl$season, levels = c("Summer", "Fall", "Winter", "Spring"))
  prey_overlap_tbl$Risk <- factor(prey_overlap_tbl$Risk, levels = c("Habitat complexity", "Black bear detected", "Bobcat detected", "Cougar detected", "Coyote detected", "Wolf detected"))
  
  #'  Effect of habitat complexity on seasonal predator-prey temporal overlap
  pred_prey_smr_coeff_plot <- ggplot(pred_prey_overlap_tbl[pred_prey_overlap_tbl$season == "Summer",], aes(x = predator, y = Dhat, group = Complexity_background_risk)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = Complexity_background_risk), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Background predation risk")) + 
    ggtitle("Effect of habitat complexity on summer predator-prey diel activity patterns") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~prey, scales = "free", space = "free")
  pred_prey_smr_coeff_plot
  
  pred_prey_fall_coeff_plot <- ggplot(pred_prey_overlap_tbl[pred_prey_overlap_tbl$season == "Fall",], aes(x = predator, y = Dhat, group = Complexity_background_risk)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = Complexity_background_risk), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Background predation risk")) + 
    ggtitle("Effect of habitat complexity on fall predator-prey diel activity patterns") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~prey, scales = "free", space = "free")
  pred_prey_fall_coeff_plot
  
  pred_prey_wtr_coeff_plot <- ggplot(pred_prey_overlap_tbl[pred_prey_overlap_tbl$season == "Winter",], aes(x = predator, y = Dhat, group = Complexity_background_risk)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = Complexity_background_risk), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Background predation risk")) + 
    ggtitle("Effect of habitat complexity on winter predator-prey diel activity patterns") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~prey, scales = "free", space = "free")
  pred_prey_wtr_coeff_plot
  
  pred_prey_sprg_coeff_plot <- ggplot(pred_prey_overlap_tbl[pred_prey_overlap_tbl$season == "Spring",], aes(x = predator, y = Dhat, group = Complexity_background_risk)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = predator, shape = Complexity_background_risk), size = 2.75, position = position_dodge(width = 0.4)) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none", shape = guide_legend(title = "Background predation risk")) + 
    ggtitle("Effect of habitat complexity on spring predator-prey diel activity patterns") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~prey, scales = "free", space = "free")
  pred_prey_sprg_coeff_plot
  
    
  #'  Effect of habitat complexity on prey temporal overlap
  prey_smr_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Season == "Summer",], aes(x = Risk, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
    geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Effect of background risk on summer predator-prey diel activity patterns") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~Species, scales = "free", space = "free")
  prey_smr_coeff_plot
  
  prey_fall_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Season == "Fall",], aes(x = Risk, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
    geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Effect of background risk on fall predator-prey diel activity patterns") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~Species, scales = "free", space = "free")
  prey_fall_coeff_plot
  
  prey_wtr_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Season == "Winter",], aes(x = Risk, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
    geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Effect of background risk on winter predator-prey diel activity patterns") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~Species, scales = "free", space = "free")
  prey_wtr_coeff_plot
  
  prey_sprg_coeff_plot <- ggplot(prey_overlap_tbl[prey_overlap_tbl$Season == "Spring",], aes(x = Risk, y = Dhat)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Risk), width = 0.3) +
    geom_point(stat = "identity", aes(col = Risk), size = 2.75) +
    scale_color_bright() + 
    ylim(0,1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,-10,-10,0)) +
    guides(color = "none") + 
    ggtitle("Effect of background risk on spring predator-prey diel activity patterns") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~Species, scales = "free", space = "free")
  prey_sprg_coeff_plot
  
  #' Save 'em - predator-prey overlap coefficients
  ggsave(pred_prey_smr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_smr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_prey_fall_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_fall_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_prey_wtr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_wtr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(pred_prey_wtr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_sprg_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  #' Save 'em - prey overlap coefficients
  ggsave(prey_smr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_prey_smr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(prey_fall_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_fall_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(prey_wtr_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_wtr_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  ggsave(prey_sprg_coeff_plot, filename = "./Outputs/Temporal Overlap/Figures/Coeff_pred_prey_sprg_bckgrd_risk.tiff", width = 7, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  #'  ----------------------------
  ####  Activity Overlap Curves  ####
  #'  ----------------------------
  #'  Plot activity curves for each species-pairing or individual species at 
  #'  cameras with low versus high levels of background risk
  overlap_pred_prey_plots <- function(dat, name1, name2, name3, dhat, y_up) {
    #'  Sample sizes for predators[1] and prey[2] where background risk is low or high
    n1low <- dat[[7]]; n1high <- dat[[9]]
    n2low <- dat[[8]]; n2high <- dat[[10]]
    spp1low <- paste0(name1, ", n = ", n1low); spp1high <- paste0(name1, " (n = ", n1high, ")")
    spp2low <- paste0(name2, ", n = ", n2low); spp2high <- paste0(name2, " (n = ", n2high, ")")
    #'  Temporal overlap between predators and prey where background risk is low or high
    dhatlow <- dhat[1,5]; dhatlowl <- dhat[1,6]; dhatlowu<- dhat[1,7]
    dhathigh <- dhat[2,5]; dhathighl <- dhat[2,6]; dhathighu<- dhat[2,7]
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
      geom_line(aes(x, y, colour =  Species.z), linetype = "dashed", lwd = 0.75) +  
      scale_x_continuous(breaks = c(0, 1.57, 3.0, 4.71, 6.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = pi/2, linetype="dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype="dotted") +
      theme_bw() +
      theme(legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      ylim(0, y_up) +
      labs(x = "Time of day", y = "Density", color = paste0("\u0394 = ", dhatp, " (", dhatpl, " - ", dhatpu, ")"), title = paste0(name3, " low")) + 
      scale_color_manual(labels = c(name3, spp1p, spp2p), values = c("black", "red", "blue")) 
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
      labs(x = "Time of day", y = "Density", color = paste0("\u0394 = ", dhata, " (", dhatal, " - ", dhatau, ")"), title = paste0(name3, " high")) + 
      scale_color_manual(labels = c(spp1a, spp2a), values = c("red", "blue"))  
    plot(overlap_high)
    
    plots <- list(overlap_low, overlap_high)
    return(plots)
  }
  ####  Cattle Activity Overlap Plots  ####
  #'  Keep track of list positions when dhat1 and dhat4 are being combined
  #'  Dhat1 for sample sizes <50, Dhat4 for sample sizes >50, fig [[1]] = present, fig [[2]] = absent
  coug_md_overPlot_g <- overlap_pred_prey_plots(pred_prey_overlap[[1]][[1]], name1 = "Cougar", name2 = "Mule deer", name3 = "Habitat complexity", dhat = coug_md_graze_out, y_up = 0.6)
  (coug_md_graze_overlap_plot <- coug_md_overPlot_g[[2]] + theme(legend.position = c(0.24, 0.92)) + coug_md_overPlot_g[[1]] + theme(legend.position = c(0.24, 0.895)))
  ggsave(coug_md_graze_overlap_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap_Plot_coug_md_cattle.tiff", width = 9, height = 6, dpi = 600, units = "in", device='tiff')
  
  
  
  
  
  
