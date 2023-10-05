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
  library(rphylopic)
  library(grid)
  library(png)
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
            axis.text.x = x_text, 
            axis.ticks.x = x_tick,
            text = element_text(size = 18), 
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      ggtitle(spp) +
      labs(x = x_lab, y = y_lab) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))  #face = "bold.italic"
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
  overlap_singlespp_plots <- function(dat, name2, name3, dhat, axis_theme, x_start, season, legend_spot, legend_title, season_colors) { #x_lab, y_lab, x_start, x_text, x_tick, axis_col, 
    #'  Define shading and line colors
    seasoncol <- season_colors[1]
    linecol1 <- season_colors[2]
    linecol2 <- season_colors[3]
    
    #'  Define axis theme
    x_lab <- axis_theme[[1]]
    y_lab <- axis_theme[[2]]
    x_text <- axis_theme[[3]]
    x_tick <- axis_theme[[4]]
    y_text <- axis_theme[[5]]
    y_tick <- axis_theme[[6]]
    axis_col <- axis_theme[[7]]
    
    #'  Create risk level labels
    spp2low <- "Low"; spp2high <- "High"
    #'  Overlap coefficient comparing activity patterns where background risk is low vs high
    dhatmu <- dhat[1,4]; dhatl <- dhat[1,5]; dhatu<- dhat[1,6]
    #'  Overlap coefficient text to be superimposed on plot
    grob <- grobTree(textGrob(paste0("\u0394 = ", dhatmu, " (", dhatl, " - ", dhatu, ")"), x = x_start, y = 0.895, hjust = 0,
                              gp = gpar(col = "black", fontsize = 16, fontface = "italic")))
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
      geom_vline(xintercept = pi/2, linetype = "dotted") +
      geom_vline(xintercept = (3*pi)/2, linetype = "dotted") +
      theme_minimal() +
      ylim(0, 0.6) + 
      theme_bw() +
      theme(axis.line = element_line(colour = axis_col),
            panel.border = element_blank(), 
            axis.text.x = x_text, 
            axis.ticks.x = x_tick,
            axis.text.y = y_text, 
            axis.ticks.y = y_tick,
            text = element_text(size = 18), 
            legend.position = legend_spot,
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = NA, fill = NA)) +
      scale_colour_manual(breaks = c("Low", "High"), labels = c(spp2low, spp2high), values = c(linecol1, linecol2), name = legend_title) +
      labs(x = x_lab, y = y_lab) + 
      annotation_custom(grob)
    plot(overlap)
    
    return(overlap)
  }
  #'  Define seasonal color schemes
  summer_colors <- c("#009E73", "#74bfa0", "#4d5c55")
  fall_colors <- c("#F0E442", "#f5e867", "#716b27")
  winter_colors <- c("#56B4E9", "#b4d9f4", "#99a3ac")
  spring_colors <- c("#CC79A7", "#e8bcd2", "#ab9da4")
  gray_scale <- c("#5d5d5d", "#a2a2a2", "#414141")
  
  #'  Define different axis themes
  axis_none <- list(x_lab = NULL, y_lab = NULL, x_text = element_blank(), x_tick = element_blank(), y_text = element_blank(), y_tick = element_blank(), axis_col = "black")
  axis_x <- list(x_lab = "Time of day", y_lab = NULL, x_text = element_text(angle = 45, vjust = 1, hjust=1), x_tick = element_line(color = "black"), y_text = element_blank(), y_tick = element_blank(), axis_col = "black")
  axis_y <- list(x_lab = NULL, y_lab = "Activity density", x_text = element_blank(), x_tick = element_blank(), y_text = element_text(), y_tick = element_line(color = "black"), axis_col = "black")
  axis_xy <- list(x_lab = "Time of day", y_lab = "Activity density", x_text = element_text(angle = 45, vjust = 1, hjust=1), x_tick = element_line(color = "black"), y_text = element_text(), y_tick = element_line(color = "black"), axis_col = "black")  
  
  ####  Mule deer summer  ####
  md_smr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[1]], name2 = "Mule deer", name3 = "TRI", dhat = md_smr_tri_out, x_lab = NULL, y_lab = NULL, season = "Summer", x_start = 0.20, x_text = element_blank(), x_tick = element_blank(), axis_col = "black", legend_spot = "none", legend_title = "Summer background risk level", season_colors = summer_colors) #seasoncol = "#009E73", linecol1 = "#74bfa0", linecol2 = "#4d5c55"
  md_smr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[2]], name2 = "Mule deer", name3 = "% Forest", dhat = md_smr_for_out, x_lab = NULL, y_lab = NULL, season = "Summer", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, season_colors = summer_colors)
  md_smr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[3]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_smr_coug_out, x_lab = NULL, y_lab = NULL, season = "Summer", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, season_colors = summer_colors)
  md_smr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[4]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_smr_wolf_out, x_lab = NULL, y_lab = NULL, season = "Summer", x_start = 0.25, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, season_colors = summer_colors)
  md_smr_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[5]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_smr_bear_out, x_lab = NULL, y_lab = NULL, season = "Summer", x_start = 0.30, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, season_colors = summer_colors)
  md_smr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[6]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_smr_bob_out, x_lab = NULL, y_lab = NULL, season = "Summer", x_start = 0.28, x_text = element_blank(), x_tick = element_blank(), legend_spot = "none", legend_title = NULL, season_colors = summer_colors)
  md_smr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[7]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_smr_coy_out, x_lab = "Time of day", y_lab = NULL, season = "Summer", x_start = 0.28, x_text = element_text(angle = 45, vjust = 1, hjust=1), x_tick = element_line(color = "black"), legend_spot = "none", legend_title = NULL, season_colors = summer_colors)
  
  ####  Mule deer fall  ####
  md_fall_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[8]], name2 = "Mule deer", name3 = "TRI", dhat = md_fall_tri_out, x_lab = NULL, y_lab = NULL, season = "Fall", x_start = 0.20, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[9]], name2 = "Mule deer", name3 = "% Forest", dhat = md_fall_for_out, x_lab = NULL, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[10]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_fall_coug_out, x_lab = NULL, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[11]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_fall_wolf_out, x_lab = NULL, y_lab = NULL, season = "Fall", x_start = 0.25, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[12]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_fall_bear_out, x_lab = NULL, y_lab = NULL, season = "Fall", x_start = 0.30, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[13]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_fall_bob_out, x_lab = NULL, y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  md_fall_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[14]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_fall_coy_out, x_lab = "Time of day", y_lab = NULL, season = "Fall", x_start = 0.28, seasoncol = "#F0E442", linecol1 = "#f5e867", linecol2 = "#716b27")
  
  ####  Mule deer winter  ####
  md_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[15]], name2 = "Mule deer", name3 = "TRI", dhat = md_wtr_tri_out, season = "Winter", x_start = 0.20, axis_theme = axis_y, legend_spot = "right", legend_title = "Winter background \nrisk level", season_colors = winter_colors)
  md_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[16]], name2 = "Mule deer", name3 = "% Forest", dhat = md_wtr_for_out, season = "Winter", x_start = 0.28, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  md_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[17]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_wtr_coug_out, season = "Winter", x_start = 0.28, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  md_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[18]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_wtr_wolf_out, season = "Winter", x_start = 0.25, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  md_wtr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[19]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_wtr_bob_out, season = "Winter", x_start = 0.38, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  md_wtr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[20]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_wtr_coy_out, season = "Winter", x_start = 0.28, axis_theme = axis_x, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  
  ####  Mule deer spring  ####
  md_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[21]], name2 = "Mule deer", name3 = "TRI", dhat = md_sprg_tri_out, x_lab = NULL, y_lab = NULL, season = "Spring", x_start = 0.20, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[22]], name2 = "Mule deer", name3 = "% Forest", dhat = md_sprg_for_out, x_lab = NULL, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[23]], name2 = "Mule deer", name3 = "cougar risk", dhat = md_sprg_coug_out, x_lab = NULL, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[24]], name2 = "Mule deer", name3 = "wolf risk", dhat = md_sprg_wolf_out, x_lab = NULL, y_lab = NULL, season = "Spring", x_start = 0.25, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[25]], name2 = "Mule deer", name3 = "black bear risk", dhat = md_sprg_bear_out, x_lab = NULL, y_lab = NULL, season = "Spring", x_start = 0.30, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[26]], name2 = "Mule deer", name3 = "bobcat risk", dhat = md_sprg_bob_out, x_lab = NULL, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  md_sprg_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[1]][[27]], name2 = "Mule deer", name3 = "coyote risk", dhat = md_sprg_coy_out, x_lab = "Time of day", y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  
  
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
  elk_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[10]], name2 = "Elk", name3 = "TRI", dhat = elk_wtr_tri_out, season = "Winter", x_start = 0.20, axis_theme = axis_y, legend_spot = "right", legend_title = "Winter background \nrisk level", season_colors = winter_colors)
  elk_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[11]], name2 = "Elk", name3 = "% Forest", dhat = elk_wtr_for_out, season = "Winter", x_start = 0.28, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  elk_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[2]][[12]], name2 = "Elk", name3 = "cougar risk", dhat = elk_wtr_coug_out, season = "Winter", x_start = 0.28, axis_theme = axis_xy, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  
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
  moose_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[11]], name2 = "Moose", name3 = "TRI", dhat = moose_wtr_tri_out, season = "Winter", x_start = 0.20, axis_theme = axis_y, legend_spot = "right", legend_title = "Winter background \nrisk level", season_colors = winter_colors)
  moose_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[12]], name2 = "Moose", name3 = "% Forest", dhat = moose_wtr_for_out, season = "Winter", x_start = 0.28, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  moose_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[13]], name2 = "Moose", name3 = "cougar risk", dhat = moose_wtr_coug_out, season = "Winter", x_start = 0.28, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  moose_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[3]][[14]], name2 = "Moose", name3 = "wolf risk", dhat = moose_wtr_wolf_out, season = "Winter", x_start = 0.25, axis_theme = axis_xy, legend_spot = "none", legend_title = NULL, season_colors = gray_scale)

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
  wtd_wtr_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[15]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_wtr_tri_out, season = "Winter", x_start = 0.20, axis_theme = axis_y, legend_spot = "right", legend_title = "Winter background \nrisk level", season_colors = winter_colors)
  wtd_wtr_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[16]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_wtr_for_out, season = "Winter", x_start = 0.28, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  wtd_wtr_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[17]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_wtr_coug_out, season = "Winter", x_start = 0.28, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  wtd_wtr_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[18]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_wtr_wolf_out, season = "Winter", x_start = 0.25, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  wtd_wtr_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[19]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_wtr_bob_out, season = "Winter", x_start = 0.28, axis_theme = axis_none, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  wtd_wtr_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[20]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_wtr_coy_out, season = "Winter", x_start = 0.28, axis_theme = axis_x, legend_spot = "none", legend_title = NULL, season_colors = winter_colors)
  
  ####  White-tailed deer spring  ####
  wtd_sprg_tri_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[21]], name2 = "White-tailed deer", name3 = "TRI", dhat = wtd_sprg_tri_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.20, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_for_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[22]], name2 = "White-tailed deer", name3 = "% Forest", dhat = wtd_sprg_for_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_coug_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[23]], name2 = "White-tailed deer", name3 = "cougar risk", dhat = wtd_sprg_coug_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_wolf_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[24]], name2 = "White-tailed deer", name3 = "wolf risk", dhat = wtd_sprg_wolf_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.25, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_bear_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[25]], name2 = "White-tailed deer", name3 = "black bear risk", dhat = wtd_sprg_bear_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.30, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_bob_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[26]], name2 = "White-tailed deer", name3 = "bobcat risk", dhat = wtd_sprg_bob_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  wtd_sprg_coy_overPlot <- overlap_singlespp_plots(prey_overlap[[4]][[27]], name2 = "White-tailed deer", name3 = "coyote risk", dhat = wtd_sprg_coy_out, y_up = 0.6, y_lab = NULL, season = "Spring", x_start = 0.28, seasoncol = "#CC79A7", linecol1 = "#e8bcd2", linecol2 = "#ab9da4")
  

  #'  ------------------------
  ####  Species silhouettes  ####
  #'  ------------------------
  #'  Silhouettes for each species from PhyloPic in two different formats (PNG & rastergrob)
  #'  Cougar, mule deer, and white-tailed deer silhouettes created by the talented 
  #'  Gabriela Palomo-Munoz and uploaded to http://phylopic.org/
  cougurlGB <- "https://images.phylopic.org/images/cbe2a3c9-2c11-4f36-a51f-8a6c8de6a420/raster/1024x489.png?v=17cfd978c92.png"
  cougimgGB <- readPNG(getURLContent(cougurlGB), native = T)
  couggrid <- rasterGrob(cougimgGB, interpolate = TRUE)
  wolfurl <- "https://images.phylopic.org/images/8cad2b22-30d3-4cbd-86a3-a6d2d004b201/raster/1024x797.png?v=16ff245cc7c.png"
  wolfimg <- readPNG(getURLContent(wolfurl), native = T)
  wolfgrid <- rasterGrob(wolfimg, interpolate = TRUE)
  boburl <- "https://images.phylopic.org/images/ab6cfd4f-aef7-40fa-b5a5-1b79b7d112aa/raster/1024x740.png?v=17a30c18af1.png"
  bobimg <- readPNG(getURLContent(boburl), native = T)
  bobgrid <- rasterGrob(bobimg, interpolate = TRUE)
  coyurl <- "https://images.phylopic.org/images/5a0398e3-a455-4ca6-ba86-cf3f1b25977a/raster/1024x894.png?v=16fe8749858.png"
  coyimg <- readPNG(getURLContent(coyurl), native = T) 
  coygrid <- rasterGrob(coyimg, interpolate = TRUE)
  coyurlGB <- "https://images.phylopic.org/images/e6a2fa4b-85df-43b4-989c-34a65ba7eee3/raster/1024x911.png?v=17f2638df97.png"
  coyimgGB <- readPNG(getURLContent(coyurlGB), native = T)
  coygridGB <- rasterGrob(coyimg, interpolate = TRUE)
  bearurl <- "https://images.phylopic.org/images/369a7880-4798-41bf-851a-ec5da17fafa3/raster/1024x753.png?v=178afd80706.png"
  bearimg <- readPNG(getURLContent(bearurl), native = T)
  beargrid <- rasterGrob(bearimg, interpolate = TRUE)
  elkfurl <- "https://images.phylopic.org/images/97f83f5e-9afe-4ce8-812e-337f506ca841/raster/1024x1005.png?v=1402ea30c27.png"
  elkfimg <- readPNG(getURLContent(elkfurl), native = T)
  elkgrid <- rasterGrob(elkfimg, interpolate = TRUE)
  elkmurl <- "https://images.phylopic.org/images/72f2f997-e474-4caf-bbd5-72fc8dbcc40d/raster/866x1024.png?v=1356f9ea6de.png"
  elkmimg <- readPNG(getURLContent(elkmurl), native = T)
  elkmgrid <- rasterGrob(elkmimg, interpolate = TRUE)
  wtdurlGB1 <- "https://images.phylopic.org/images/8569838c-c725-4772-b0a3-b5eb04baaada/raster/1024x850.png?v=17cfdbaf920.png"
  wtdimgGB1 <- readPNG(getURLContent(wtdurlGB1), native = T)
  wtdgrid <- rasterGrob(wtdimgGB1, interpolate = TRUE)
  wtdurlGB2 <- "https://images.phylopic.org/images/6038e80c-398d-47b2-9a69-2b9edf436f64/raster/1023x1024.png?v=17cfdb9f8b6.png"
  wtdimgGB2 <- readPNG(getURLContent(wtdurlGB2), native = T)
  wtdgridGB2 <- rasterGrob(wtdimgGB2, interpolate = TRUE)
  
  
  #'  -------------------------
  ####  Patchwork mega figure  ####
  #'  -------------------------
  #'  Patchwork all individual figures together to create one massive figure
  fig1 <- elk_all_seasons + moose_all_seasons + md_all_seasons + wtd_all_seasons + 
    elk_wtr_tri_overPlot + moose_wtr_tri_overPlot + md_wtr_tri_overPlot + wtd_wtr_tri_overPlot + 
    elk_wtr_coug_overPlot + moose_wtr_coug_overPlot + md_wtr_coug_overPlot + plot_spacer() + 
    plot_spacer() + moose_wtr_wolf_overPlot + md_wtr_wolf_overPlot + wtd_wtr_wolf_overPlot + 
    guide_area() + plot_spacer() + md_wtr_coy_overPlot + wtd_wtr_coy_overPlot + theme(legend.position = "none") +
    plot_layout(ncol = 4, nrow = 5, guides = 'collect') + plot_annotation(tag_levels = 'a')
  plot(fig1)
  ggsave(fig1, filename = "./Outputs/Temporal Overlap/Figures/Overlap Plots/Figure1_prey_winter_activity_overlap.tiff", width = 17, height = 18, dpi = 600, units = "in", device='tiff')
  