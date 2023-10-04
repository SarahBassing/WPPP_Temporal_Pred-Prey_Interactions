  #'  =====================================================
  #'  Temporal overlap plots: Species-specific activity curves
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  August 2022
  #'  =====================================================
  #'  Script pulls in outputs from temporal overlap analyses, formats results for 
  #'  easier plotting, then creates figures depicting species-specific seasonal 
  #'  activity curves at sites with high vs low levels of background risk. 
  #'  =====================================================
  
  #'  Libraries
  library(ggplot2)
  library(khroma)
  library(grid)
  library(patchwork)
  library(tidyverse)  
  
  #'  Load output from temporal overlap analysis
  # load("./Outputs/Temporal Overlap/PredPrey_TRI_Overlap_2022-09-19.RData")
  # load("./Outputs/Temporal Overlap/PredPrey_PercForest_Overlap_2022-09-19.RData")
  load("./Outputs/Temporal Overlap/PreyOnly_TRI_Forest_Pred_Overlap_2022-09-26.RData")
  # load("./Outputs/Temporal Overlap/PredOnly_TRI_Forest_Overlap_2022-10-15.RData")
  load("./Outputs/Temporal OVerlap/Prey_Overall_Activity_2023-09-13.RData")
  
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
  
  
  #'  ---------------------------------
  ####  Seasonal Prey Activity Curves  ####
  #'  ---------------------------------
  #'  Plot overall activity curves for each prey species per season
  singlespp_activity_df <- function(dat, spp) {
    dat[[1]]$Season <- "Summer"
    dat[[2]]$Season <- "Fall"
    dat[[3]]$Season <- "Winter"
    dat[[4]]$Season <- "Spring"
    allactivity <- rbind(dat[[1]], dat[[2]], dat[[3]], dat[[4]])
    allactivity$Species <- spp
    colnames(allactivity) <- c("Time", "Density", "Season", "Species")
    allactivity <- allactivity %>%
      filter(Time >= 0 & Time <= 24.0) %>%
      mutate(Season = factor(Season, levels = c("Summer", "Fall", "Winter", "Spring")))
    return(allactivity)
  }
  ####  Mule deer summer  ####
  md_activity <- singlespp_activity_df(prey_activity[[1]], spp = "Mule deer")
  elk_activity <- singlespp_activity_df(prey_activity[[2]], spp = "Elk")
  moose_activity <- singlespp_activity_df(prey_activity[[3]], spp = "Moose")
  wtd_activity <- singlespp_activity_df(prey_activity[[4]], spp = "White-tailed deer")
  spp_activity_df <- rbind(md_activity, elk_activity, moose_activity, wtd_activity)
  
  #'  Plot facetted figure for all species together
  spp_activity_plot <- ggplot(spp_activity_df, aes(Time, Density, colour = Season)) + #, linetype = Season
    geom_line(lwd = 0.75) + 
    scale_color_manual(values = c("#009E73", "#F0E442", "#56B4E9", "#CC79A7")) +
    # scale_linetype_manual(values = c("solid", "dashed", "dotted", "longdash")) +
    scale_x_continuous(breaks = c(0, 6.0, 12.0, 18.0, 24.0),
                       labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
    facet_wrap(.~Species, scales = "free_y", ncol = 4) +
    geom_vline(xintercept = 6.0, linetype="dotted") +
    geom_vline(xintercept = 18.0, linetype="dotted") +
    theme_bw() +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(colour = NA, fill = NA)) +
    theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Time of day", y = "Density") +
    theme(text = element_text(size = 14), legend.text = element_text(size = 14)) 
  plot(spp_activity_plot)  
  # ggsave(spp_activity_plot, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Ungulate_Overall_Activity_Plots.tiff", width = 15, height = 4, dpi = 600, units = "in", device='tiff')
  
  #'  Keep each plot separate 
  seasonal_activity_plot <- function(activity_df, spp, x_lab, y_lab, x_text, x_tick) {
    all_seasons <- ggplot(activity_df, aes(Time, Density, colour = Season)) + #, linetype = Season
      geom_line(lwd = 0.75) + 
      scale_color_manual(values = c("#009E73", "#F0E442", "#56B4E9", "#CC79A7")) +
      scale_x_continuous(breaks = c(0, 6.0, 12.0, 18.0, 24.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = 6.0, linetype="dotted") +
      geom_vline(xintercept = 18.0, linetype="dotted") +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.border = element_blank(), 
            axis.text.x = x_text, #element_text(angle = 45, vjust = 1, hjust=1),
            axis.ticks.x = x_tick,
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      ggtitle(spp) +
      labs(x = x_lab, y = y_lab) +
      theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  #face = "bold.italic"
            text = element_text(size = 14), legend.text = element_text(size = 14)) 
  }
  (md_all_seasons <- seasonal_activity_plot(md_activity, spp = "Mule deer", x_lab = NULL, y_lab = NULL, x_text = element_blank(), x_tick = element_blank()))
  (elk_all_seasons <- seasonal_activity_plot(elk_activity, spp = "Elk", x_lab = NULL, y_lab = "Activity density", x_text = element_blank(), x_tick = element_blank()))
  (moose_all_seasons <- seasonal_activity_plot(moose_activity, spp = "Moose", x_lab = NULL, y_lab = NULL, x_text = element_blank(), x_tick = element_blank()))
  (wtd_all_seasons <- seasonal_activity_plot(wtd_activity, spp = "White-tailed deer", x_lab = NULL, y_lab = NULL, x_text = element_blank(), x_tick = element_blank()))
 
  
  #'  --------------------------------
  ####  Prey Activity Overlap Curves  ####
  #'  --------------------------------
  #'  Plot activity curves for each prey species at cameras with low versus high 
  #'  levels of background risk
  #'  Seasonal fill colors for summer, fall, winter, spring --> c("#009E73", "#F0E442", "#56B4E9", "#CC79A7")
  overlap_singlespp_plots <- function(dat, name2, name3, dhat, x_lab, y_lab, season, x_start, x_text, x_tick, legend_spot, legend_title, seasoncol, linecol1, linecol2) { #y_up
    #'  Create risk level labels
    # spp2low <- paste0("Low ", name3); spp2high <- paste0("High ", name3)
    spp2low <- "Low"; spp2high <- "High"
    #'  Overlap coefficient comparing activity patterns where background risk is low vs high
    dhatmu <- dhat[1,4]; dhatl <- dhat[1,5]; dhatu<- dhat[1,6]
    #'  Overlap coefficient text to be superimposed on plot
    grob <- grobTree(textGrob(paste0("\u0394 = ", dhatmu, " (", dhatl, " - ", dhatu, ")"), x = x_start, y = 0.895, hjust = 0,
                              gp = gpar(col = "black", fontsize = 14, fontface = "italic")))
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
                alpha = 0.25, fill = seasoncol, color = NA) +
      scale_x_continuous(breaks = c(0, 1.57, 3.0, 4.71, 6.0),
                         labels = c('Midnight', 'Dawn', 'Noon', 'Dusk', 'Midnight')) +
      geom_vline(xintercept = pi/2, linetype="dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype="dotted") +
      theme_minimal() +
      ylim(0, 0.6) + #y_up
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.border = element_blank(), 
            axis.text.x = x_text, 
            axis.ticks.x = x_tick,
            text = element_text(size = 14), 
            legend.position = legend_spot,
            legend.text = element_text(size = 14),
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      scale_colour_manual(breaks = c("Low", "High"), labels = c(spp2low, spp2high), values = c(linecol1, linecol2), name = legend_title) +
      labs(x = x_lab, y = y_lab) + 
      annotation_custom(grob)
    plot(overlap)
    
    return(overlap)
  }
  ####  Mule deer summer  ####
  md_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[1]], name2 = "Mule deer", name3 = "TRI", dhat = md_smr_tri_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.20, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  md_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[2]], name2 = "Mule deer", name3 = "% Forest", dhat = md_smr_for_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  md_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[3]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_smr_coug_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  md_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[4]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_smr_wolf_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.25, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  md_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[5]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_smr_bear_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.30, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  md_smr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[6]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_smr_bob_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  md_smr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[7]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_smr_coy_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  
  ####  Mule deer fall  ####
  md_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[8]], name2 = "Mule deer", name3 = "TRI", dhat = md_fall_tri_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.20, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[9]], name2 = "Mule deer", name3 = "% Forest", dhat = md_fall_for_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[10]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_fall_coug_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[11]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_fall_wolf_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.25, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[12]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_fall_bear_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.30, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[13]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_fall_bob_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[14]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_fall_coy_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  
  ####  Mule deer winter  ####
  md_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[15]], name2 = "Mule deer", name3 = "TRI", dhat = md_wtr_tri_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.20, x_text = element_blank(), x_tick = element_blank(), legend_spot = "bottom", legend_title = "Winter background risk level", seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  md_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[16]], name2 = "Mule deer", name3 = "% Forest", dhat = md_wtr_for_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  md_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[17]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_wtr_coug_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  md_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[18]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_wtr_wolf_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.25, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  md_wtr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[19]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_wtr_bob_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.38, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  md_wtr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[20]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_wtr_coy_out, x_lab = "Time of day", y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_text(angle = 45, vjust = 1, hjust=1), x_tick = element_line(color = "black"), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  
  ####  Mule deer spring  ####
  md_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[21]], name2 = "Mule deer", name3 = "TRI", dhat = md_sprg_tri_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.20, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[22]], name2 = "Mule deer", name3 = "% Forest", dhat = md_sprg_for_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[23]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_sprg_coug_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[24]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_sprg_wolf_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.25, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[25]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_sprg_bear_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.30, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[26]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_sprg_bob_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[27]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_sprg_coy_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  
  
  ####  Elk summer  ####
  elk_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[1]], name2 = "Elk", name3 = "TRI", dhat = elk_smr_tri_out, y_up = 0.6, y_lab = "Activity density", season = "Summer", x_start = 0.20, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  elk_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[2]], name2 = "Elk", name3 = "% Forest", dhat = elk_smr_for_out, y_up = 0.6, y_lab = "Activity density", season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  elk_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[3]], name2 = "Elk", name3 = "cougar risk", dhat = elk_smr_coug_out, y_up = 0.6, y_lab = "Activity density", season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  elk_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[4]], name2 = "Elk", name3 = "wolf risk", dhat = elk_smr_wolf_out, y_up = 0.6, y_lab = "Activity density", season = "Summer", x_start = 0.25, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  elk_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[5]], name2 = "Elk", name3 = "black bear risk", dhat = elk_smr_bear_out, y_up = 0.6, y_lab = "Activity density", season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  
  ####  Elk fall  ####
  elk_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[6]], name2 = "Elk", name3 = "TRI", dhat = elk_fall_tri_out, y_up = 0.6, y_lab = "Activity density", season = "Fall", x_start = 0.20, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  elk_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[7]], name2 = "Elk", name3 = "% Forest", dhat = elk_fall_for_out, y_up = 0.6, y_lab = "Activity density", season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  elk_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[8]], name2 = "Elk", name3 = "cougar risk", dhat = elk_fall_coug_out, y_up = 0.6, y_lab = "Activity density", season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  elk_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[9]], name2 = "Elk", name3 = "black bear risk", dhat = elk_fall_bear_out, y_up = 0.6, y_lab = "Activity density", season = "Fall", x_start = 0.30, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  
  ####  Elk winter  ####
  elk_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[10]], name2 = "Elk", name3 = "TRI", dhat = elk_wtr_tri_out, x_lab = NULL, y_lab = "Activity density", season = "Winter", x_start = 0.20, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  elk_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[11]], name2 = "Elk", name3 = "% Forest", dhat = elk_wtr_for_out, x_lab = NULL, y_lab = "Activity density", season = "Winter", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  elk_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[12]], name2 = "Elk", name3 = "cougar risk", dhat = elk_wtr_coug_out, x_lab = "Time of day", y_lab = "Activity density", season = "Winter", x_start = 0.28, x_text = element_text(angle = 45, vjust = 1, hjust=1), x_tick = element_line(color = "black"), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  
  ####  Elk spring  ####
  elk_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[13]], name2 = "Elk", name3 = "TRI", dhat = elk_sprg_tri_out, y_up = 0.6, y_lab = "Activity density", season = "Spring", x_start = 0.20, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  elk_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[14]], name2 = "Elk", name3 = "% Forest", dhat = elk_sprg_for_out, y_up = 0.6, y_lab = "Activity density", season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  elk_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[15]], name2 = "Elk", name3 = "cougar risk", dhat = elk_sprg_coug_out, y_up = 0.6, y_lab = "Activity density", season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  elk_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[16]], name2 = "Elk", name3 = "black bear risk", dhat = elk_sprg_bear_out, y_up = 0.6, y_lab = "Activity density", season = "Spring", x_start = 0.30, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  
  
  ####  Moose summer  ####
  moose_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[1]], name2 = "Moose", name3 = "TRI", dhat = moose_smr_tri_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.20, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  moose_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[2]], name2 = "Moose", name3 = "% Forest", dhat = moose_smr_for_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  moose_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[3]], name2 = "Moose", name3 = "cougar risk", dhat = moose_smr_coug_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  moose_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[4]], name2 = "Moose", name3 = "wolf risk", dhat = moose_smr_wolf_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.25, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  moose_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[5]], name2 = "Moose", name3 = "black bear risk", dhat = moose_smr_bear_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.30, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  
  ####  Moose fall  ####
  moose_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[6]], name2 = "Moose", name3 = "TRI", dhat = moose_fall_tri_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.20, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  moose_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[7]], name2 = "Moose", name3 = "% Forest", dhat = moose_fall_for_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  moose_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[8]], name2 = "Moose", name3 = "cougar risk", dhat = moose_fall_coug_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  moose_fall_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[9]], name2 = "Moose", name3 = "wolf risk", dhat = moose_fall_wolf_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.25, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  moose_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[10]], name2 = "Moose", name3 = "black bear risk", dhat = moose_fall_bear_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.30, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
 
  ####  Moose winter  ####
  moose_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[11]], name2 = "Moose", name3 = "TRI", dhat = moose_wtr_tri_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.20, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  moose_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[12]], name2 = "Moose", name3 = "% Forest", dhat = moose_wtr_for_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  moose_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[13]], name2 = "Moose", name3 = "cougar risk", dhat = moose_wtr_coug_out, x_lab = "Time of day", y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_text(angle = 45, vjust = 1, hjust=1), x_tick = element_line(color = "black"), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  moose_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[14]], name2 = "Moose", name3 = "wolf risk", dhat = moose_wtr_wolf_out, x_lab = "Time of day", y_lab = NULL, season = "Winter", x_start = 0.25, x_text = element_text(angle = 45, vjust = 1, hjust=1), x_tick = element_line(color = "black"), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")

  ####  Moose spring  ####
  moose_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[15]], name2 = "Moose", name3 = "TRI", dhat = moose_sprg_tri_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.20, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  moose_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[16]], name2 = "Moose", name3 = "% Forest", dhat = moose_sprg_for_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  moose_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[17]], name2 = "Moose", name3 = "cougar risk", dhat = moose_sprg_coug_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  moose_sprg_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[18]], name2 = "Moose", name3 = "wolf risk", dhat = moose_sprg_wolf_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.25, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  moose_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[19]], name2 = "Moose", name3 = "black bear risk", dhat = moose_sprg_bear_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.30, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
 
  ####  White-tailed deer summer  ####
  wtd_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[1]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_smr_tri_out, y_up = 0.6, y_lab = NULL, y_lab = NULL, season = "Summer", x_start = 0.20, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  wtd_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[2]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_smr_for_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  wtd_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[3]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_smr_coug_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  wtd_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[4]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_smr_wolf_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.25, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  wtd_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[5]], name2 = "White-tailed deer", name3 = "black bear risk", dhat = wtd_smr_bear_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.30, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  wtd_smr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[6]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_smr_bob_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  wtd_smr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[7]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_smr_coy_out, y_up = 0.6, y_lab = NULL, season = "Summer", x_start = 0.28, seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55")
  
  ####  White-tailed deer fall  ####
  wtd_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[8]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_fall_tri_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.20, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  wtd_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[9]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_fall_for_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  wtd_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[10]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_fall_coug_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  wtd_fall_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[11]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_fall_wolf_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.25, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  wtd_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[12]], name2 = "White-tailed deer", name3 = "black bear risk", dhat = wtd_fall_bear_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.30, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  wtd_fall_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[13]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_fall_bob_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  wtd_fall_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[14]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_fall_coy_out, y_up = 0.6, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  
  ####  White-tailed deer winter  ####
  wtd_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[15]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_wtr_tri_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.20, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  wtd_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[16]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_wtr_for_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  wtd_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[17]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_wtr_coug_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  wtd_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[18]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_wtr_wolf_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.25, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  wtd_wtr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[19]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_wtr_bob_out, x_lab = NULL, y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  wtd_wtr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[20]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_wtr_coy_out, x_lab = "Time of day", y_lab = NULL, season = "Winter", x_start = 0.28, x_text = element_text(angle = 45, vjust = 1, hjust=1), x_tick = element_line(color = "black"), legend_spot = "none", legend_title = NULL, seasoncol = "#56B4E9", linecol1 = "#b4d9f4", linecol2 = "#99a3ac")
  
  ####  White-tailed deer spring  ####
  wtd_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[21]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_sprg_tri_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.20, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[22]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_sprg_for_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[23]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_sprg_coug_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[24]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_sprg_wolf_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.25, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[25]], name2 = "White-tailed deer", name3 = "black bear risk", dhat = wtd_sprg_bear_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.30, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[26]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_sprg_bob_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[27]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_sprg_coy_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  

  #'  -------------------------
  ####  Patchwork mega figure  ####
  #'  -------------------------
  #'  Patchwork all individual figures together to create one massive figure
  fig1 <- elk_all_seasons + moose_all_seasons + md_all_seasons + wtd_all_seasons +
    elk_wtr_tri_overPlot + moose_wtr_tri_overPlot + md_wtr_tri_overPlot + wtd_wtr_tri_overPlot +
    elk_wtr_coug_overPlot + moose_wtr_coug_overPlot + md_wtr_coug_overPlot + plot_spacer() +
    plot_spacer() + plot_spacer() + md_wtr_wolf_overPlot + wtd_wtr_wolf_overPlot +
    plot_spacer() + plot_spacer() + md_wtr_coy_overPlot + wtd_wtr_coy_overPlot +
    plot_layout(ncol = 4, nrow = 5, guides = 'collect') & theme(legend.position = "bottom")
  plot(fig1)
  ggsave(fig1, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Figure1_prey_winter_activity_overlap.tiff", width = 17, height = 15, dpi = 600, units = "in", device='tiff')
  