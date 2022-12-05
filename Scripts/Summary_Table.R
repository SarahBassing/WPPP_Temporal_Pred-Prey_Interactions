  #'  =================================
  #'  Result tables for publication
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  =================================
  
  #'  Load libraries
  library(camtrapR)
  library(chron)
  library(lubridate)
  library(tidyverse)
  
  #'  Read in data
  #'  Camera stations: DON'T format the dates in the cam_stations file yet!
  cam_stationsYr1 <- read.csv("./Data/All_Camera_Stations_18-19_updated_1.21.21.csv")
  cam_stationsYr2 <- read.csv("./Data/All_Camera_Stations_19-20.csv")
  cam_stationsYr3 <- read.csv("./Data/All_Camera_Stations_20-21.csv")
  cam_stations <- rbind(cam_stationsYr1, cam_stationsYr2, cam_stationsYr3)
  #'  Seasonal detection events
  seasonal_detections <- read.csv("./Outputs/seasonal_detections_2022-09-19.csv") %>%
    mutate(Year = ifelse(Date < "2019-06-01", "Year1", "Year2"),
           Year = ifelse(Date > "2020-05-31", "Year3", Year))
  
  #'  Camera operation table
  camop_problem <- cameraOperation(CTtable = cam_stations,
                                   stationCol = "CameraLocation",
                                   setupCol = "Setup_date",
                                   retrievalCol = "Retrieval_date",
                                   hasProblems = TRUE,
                                   dateFormat = "%m/%d/%Y", # Define date format here!
                                   writecsv = FALSE) 
  
  probs <- as.data.frame(camop_problem)
  
  #'  Detection history
  DH <- function(images, spp, start_date) {
    det_hist <- detectionHistory(recordTable = images,
                                 camOp = camop_problem,
                                 stationCol = "CameraLocation",
                                 speciesCol = "Species",
                                 recordDateTimeCol = "DateTime",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 species = spp,
                                 occasionLength = 7,
                                 day1 = start_date, 
                                 # datesAsOccasionNames = TRUE,
                                 # occasionStartTime = 12, # starts at noon
                                 timeZone = "America/Los_Angeles",
                                 output = "binary",
                                 includeEffort = TRUE,
                                 scaleEffort = FALSE,
                                 # writecsv = TRUE,
                                 outDir = "./Data/Detection_Histories")
    
    return(det_hist)
  }
  ####  BOBCATS  ####
  bob_yr1 <- seasonal_detections[seasonal_detections$Species == "Bobcat" & seasonal_detections$Year == "Year1",]
  bob_yr2 <- seasonal_detections[seasonal_detections$Species == "Bobcat" & seasonal_detections$Year == "Year2",]
  bob_yr3 <- seasonal_detections[seasonal_detections$Species == "Bobcat" & seasonal_detections$Year == "Year3",]
  #'  Summer season
  bob_smr1 <- DH(bob_yr1, "Bobcat", "2018-06-01")
  DH_bob_smr1 <- bob_smr1[[1]][1:125,1:17]   
  bob_smr2 <- DH(bob_yr2, "Bobcat", "2019-06-01")
  DH_bob_smr2 <- bob_smr2[[1]][126:242,1:17] 
  bob_smr3 <- DH(bob_yr3, "Bobcat", "2020-06-01")
  DH_bob_smr3 <- bob_smr3[[1]][243:361,1:17] 
  DH_bob_smr1820 <- rbind(DH_bob_smr1, DH_bob_smr2, DH_bob_smr3)
  #' #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  #' #'  All detections are missing at sites: 10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336
  #' DH_bob_smr1820 <- DH_bob_smr1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  
  #'  Fall season
  bob_fll1 <- DH(bob_yr1, "Bobcat", "2018-10-01")
  DH_bob_fll1 <- bob_fll1[[1]][1:125,1:9]   
  bob_fll2 <- DH(bob_yr2, "Bobcat", "2019-10-01")
  DH_bob_fll2 <- bob_fll2[[1]][126:242,1:9] 
  bob_fll3 <- DH(bob_yr3, "Bobcat", "2020-10-01")
  DH_bob_fll3 <- bob_fll3[[1]][243:361,1:9] 
  DH_bob_fll1820 <- rbind(DH_bob_fll1, DH_bob_fll2, DH_bob_fll3)
  #' #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  #' #'  All detections are missing at sites: 10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336
  #' DH_bob_fll1820 <- DH_bob_fll1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  
  #'  Winter season
  bob_wtr1 <- DH(bob_yr1, "Bobcat", "2018-12-01")
  DH_bob_wtr1 <- bob_wtr1[[1]][1:125,1:17]   
  bob_wtr2 <- DH(bob_yr2, "Bobcat", "2019-12-01")
  DH_bob_wtr2 <- bob_wtr2[[1]][126:242,1:17] 
  bob_wtr3 <- DH(bob_yr3, "Bobcat", "2020-12-01")
  DH_bob_wtr3 <- bob_wtr3[[1]][243:361,1:17] 
  DH_bob_wtr1820 <- rbind(DH_bob_wtr1, DH_bob_wtr2, DH_bob_wtr3)
  #' #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  #' #'  All detections are missing at sites: 10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336
  #' DH_bob_wtr1820 <- DH_bob_wtr1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  
  #'  Spring season
  bob_sprg1 <- DH(bob_yr1, "Bobcat", "2018-10-01")
  DH_bob_sprg1 <- bob_sprg1[[1]][1:125,1:9]   
  bob_sprg2 <- DH(bob_yr2, "Bobcat", "2019-10-01")
  DH_bob_sprg2 <- bob_sprg2[[1]][126:242,1:9] 
  bob_sprg3 <- DH(bob_yr3, "Bobcat", "2020-10-01")
  DH_bob_sprg3 <- bob_sprg3[[1]][243:361,1:9] 
  DH_bob_sprg1820 <- rbind(DH_bob_sprg1, DH_bob_sprg2, DH_bob_sprg3)
  #' #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  #' #'  All detections are missing at sites: 10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336
  #' DH_bob_sprg1820 <- DH_bob_sprg1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  
  #'  Sampling effort
  Effort_smr18 <- bob_smr1[[2]][1:125,1:17]
  Effort_smr19 <- bob_smr2[[2]][126:242,1:17]
  Effort_smr20 <- bob_smr3[[2]][243:361,1:17]
  Effort_fll18 <- bob_fll1[[2]][1:125,1:9]
  Effort_fll19 <- bob_fll2[[2]][126:242,1:9]
  Effort_fll20 <- bob_fll3[[2]][243:361,1:9]
  Effort_wtr18 <- bob_wtr1[[2]][1:125,1:17]
  Effort_wtr19 <- bob_wtr2[[2]][126:242,1:17]
  Effort_wtr20 <- bob_wtr3[[2]][243:361,1:17]
  Effort_sprg18 <- bob_sprg1[[2]][1:125,1:9]
  Effort_sprg19 <- bob_sprg2[[2]][126:242,1:9]
  Effort_sprg20 <- bob_sprg3[[2]][243:361,1:9]
  
  Effort_smr1820 <- rbind(Effort_smr18, Effort_smr19, Effort_smr20)
  Effort_fll1820 <- rbind(Effort_fll18, Effort_fll19, Effort_fll20)
  Effort_wtr1820 <- rbind(Effort_wtr18, Effort_wtr19, Effort_wtr20)
  Effort_sprg1820 <- rbind(Effort_sprg18, Effort_sprg19, Effort_sprg20)
  
  Effort_Yr1 <- cbind(Effort_smr18, Effort_fll18, Effort_wtr18, Effort_sprg18)
  Effort_Yr2 <- cbind(Effort_smr19, Effort_fll19 , Effort_wtr19, Effort_sprg19)
  Effort_Yr3 <- cbind(Effort_smr20, Effort_fll20, Effort_wtr20, Effort_sprg20)
  
  #'  Summary stats on trap nights (active sites only)
  #'  Total number of trap nights per season (across all years)
  (trapnights_smr <- sum(Effort_smr1820, na.rm = T)) #29494.5
  (trapnights_fll <- sum(Effort_fll1820, na.rm = T))   #20729
  (trapnights_wtr <- sum(Effort_wtr1820, na.rm = T)) #38089
  (trapnights_sprg <- sum(Effort_sprg1820, na.rm = T))   #20729
  trapnights <- c(trapnights_smr, trapnights_fll, trapnights_wtr, trapnights_sprg)
  print(mu_trapnight <- mean(trapnights)) #27260.38
  print(se_trapnights <- sd(trapnights)/length(trapnights)) #2079.504
  
  #'  Total number of trap nights per year (across all seasons)
  (trapnights_yr1 <- sum(Effort_Yr1, na.rm = T)) #36524
  (trapnights_yr2 <- sum(Effort_Yr2, na.rm = T)) #36382.5
  (trapnights_yr3 <- sum(Effort_Yr3, na.rm = T)) #36135
  trapnights_yr <- c(trapnights_yr1, trapnights_yr2, trapnights_yr3)
  print(mu_trapnight_yr <- mean(trapnights_yr)) #36347.17
  print(se_trapnights_yr <- sd(trapnights_yr)/length(trapnights_yr)) #65.63077
  
  #'  Remove rows of all NAs (inactive camera stations)
  eff_smr1820 <- Effort_smr1820[rowSums(is.na(Effort_smr1820)) != ncol(Effort_smr1820), ]
  eff_fll1820 <- Effort_fll1820[rowSums(is.na(Effort_fll1820)) != ncol(Effort_fll1820), ]
  eff_wtr1820 <- Effort_wtr1820[rowSums(is.na(Effort_wtr1820)) != ncol(Effort_wtr1820), ]
  eff_sprg1820 <- Effort_sprg1820[rowSums(is.na(Effort_sprg1820)) != ncol(Effort_sprg1820), ]
  #'  Number of active cameras per season
  print(ncams_smr <- nrow(eff_smr1820)) #350, ~116.66/year
  print(ncams_fll <- nrow(eff_fll1820))  #342, ~114/year
  print(ncams_wtr <- nrow(eff_wtr1820)) #333, ~111/year
  print(ncams_sprg <- nrow(eff_sprg1820))  #342, ~114/year
  #'  Mean number of active cameras per season
  ncams <- c(ncams_smr, ncams_fll, ncams_wtr, ncams_sprg)
  print(mean(ncams))
  
  
  
  ####  Photo-capture Summary Tables  ####
  #'  --------------------------------
  #'  Number of independent detections per season from camera traps
  skinny_dets <- seasonal_detections %>%
    mutate(Category = ifelse(Species == "Bobcat" | Species == "Black Bear" | 
                               Species == "Cougar" | Species == "Coyote" | 
                               Species == "Lynx" | Species == "Wolf", "Predator", "Other"),
           Category = ifelse(Species == "Elk" | Species == "Mule Deer" | 
                               Species == "Moose" | Species == "White-tailed Deer", "Prey", Category),
           Season = factor(Season, levels = c("Summer", "Fall", "Winter", "Spring"))) %>%
    filter(Category != "Other") %>%
    filter(Species != "Lynx") 
  ndet <- skinny_dets %>%
    group_by(Season, Species) %>%
    summarise(n = n()) %>%
    ungroup() 
  
  #'  Average number of detections per species across seasons
  summary_dets <- group_by(ndet, Species) %>% 
    summarize(mu_dets = mean(n), sd = sd(n), se_dets = sd(n)/sqrt(n())) %>% 
    ungroup()
  print(mean_dets <- mean(summary_dets$mu_dets)) #1031.361
  print(se_dets <- (sd(summary_dets$mu_dets)/(sqrt(nrow(summary_dets))))) #474.609
  colnames(summary_dets) <- c("Species", "Mean detections", "SD", "SE")

  #'  Proportion of cameras where each species was detected
  perc_cams <- skinny_dets %>%
    dplyr::select(CameraLocation, Species, Season) %>% 
    group_by(Species, Season, CameraLocation) %>%
    filter(row_number(CameraLocation) == 1) %>%
    ungroup() %>%
    group_by(Species, Season) %>%
    summarise(ncams = n()) %>%
    ungroup() %>%
    mutate(
      propcams = ifelse(Season == "Summer", ncams/ncams_smr, ncams/ncams_fll),
      propcams = ifelse(Season == "Winter", ncams/ncams_wtr, propcams),
      propcams = ifelse(Season == "Spring", ncams/ncams_sprg, propcams),
      propcams = round(propcams, 2)) %>%
    dplyr::select(-ncams)
  colnames(perc_cams) <- c("Species", "Season", "Proportion of cameras")
  
  #'  Create and save summery table of detection data for publication
  detection_summary_data <- full_join(ndet, perc_cams, by = c("Species", "Season")) %>%
    rename("Independent detections (n)" = "n")
  write.csv(detection_summary_data, "./Outputs/summary_table_detections.csv")

  
  
  ####  Time-between-detections Summary Tables  ####
  #'  -------------------------------------------
  tbd_pred.prey <- read.csv("./Outputs/tbd_pred.prey_NoOutliers.csv")
  tbd_ungulate <- read.csv("./Outputs/tbd_ungulate_NoOutliers.csv")
  tbd_conspif <- read.csv("./Outputs/tbd_conspif_NoOutliers.csv")
  
  tbd_pp <- tbd_pred.prey %>%
      dplyr::select(c(CameraLocation, Species, PredatorID, Study_Area, Season, Year, 
                      backgroundRisk_TRI, backgroundRisk_For, tbd_min_round, tbd_min)) %>%
      filter(PredatorID != "Lynx") %>%
      #'  Change units of time
      mutate(Pair = paste0(Species, " - ", PredatorID),
             Season = factor(Season, levels = c("Summer", "Fall", "Winter", "Spring")))
  nrow(tbd_pp) # 2677 after outliers were removed (2759 including outliers)
  tbd_ung <- tbd_ungulate %>%
    dplyr::select(c(CameraLocation, Species, PreviousDet, Study_Area, Season, Year, 
                    backgroundRisk_TRI, backgroundRisk_For, tbd_min_round, tbd_min)) %>%
    #'  Change units of time
    mutate(Pair = paste0(Species, " - ", PreviousDet),
           Season = factor(Season, levels = c("Summer", "Fall", "Winter", "Spring")))
  nrow(tbd_ung) # 1747 after outliers were removed (1809 including outliers)
  tbd_con <- tbd_conspif %>%
    dplyr::select(c(CameraLocation, Species, Study_Area, Season, Year, 
                    backgroundRisk_TRI, backgroundRisk_For, tbd_min_round, tbd_min)) %>%
    #'  Change units of time
    mutate(Pair = paste0(Species, " - ", Species),
           Season = factor(Season, levels = c("Summer", "Fall", "Winter", "Spring")))
  nrow(tbd_con) # 22895 (23,614 including outliers)
  
  #'  Summarize number of observations per species pair and season
  ntbd_pp <- tbd_pp %>%
    group_by(Season, Pair) %>%
    summarise(n = n()) %>%
    ungroup()
  mu_season_tbd_pp <- ntbd_pp %>%
    group_by(Season) %>%
    summarise(mu = mean(n), se = (sd(n))/sqrt(length(n))) %>%
    ungroup()
  #'  Average observations per species
  mu_spp_tbd_pp <- tbd_pp %>%
    group_by(Species) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    summarise(mu = mean(n), se = (sd(n))/sqrt(length(n)))
  
  ntbd_ung <- tbd_ung %>%
    group_by(Season, Pair) %>%
    summarise(n = n()) %>%
    ungroup()
  mu_season_tbd_ung <- ntbd_ung %>%
    group_by(Season) %>%
    summarise(mu = mean(n), se = (sd(n))/sqrt(length(n))) %>%
    ungroup()
  #'  Average observations per species
  mu_spp_tbd_ung <- tbd_ung %>%
    group_by(Species) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    summarise(mu = mean(n), se = (sd(n))/sqrt(length(n)))
  
  ntbd_con <- tbd_con %>%
    group_by(Season, Pair) %>%
    summarise(n = n()) %>%
    ungroup()
  mu_season_tbd_con <- ntbd_con %>%
    group_by(Season) %>%
    summarise(mu = mean(n), se = (sd(n))/sqrt(length(n))) %>%
    ungroup()
  #'  Average observations per species
  mu_spp_tbd_con <- tbd_con %>%
    group_by(Species) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    summarise(mu = mean(n), se = (sd(n))/sqrt(length(n)))
  
  #'  Split by number of observations by season
  ntbd_pp_smr <- ntbd_pp %>% filter(Season == "Summer") %>% 
    mutate(Analysis = "Predator-Prey") %>%
    rename("Summer (n)" = "n") %>% dplyr::select(-Season)
  ntbd_pp_fall <- ntbd_pp %>% filter(Season == "Fall") %>%  
    mutate(Analysis = "Predator-Prey") %>%
    rename("Fall (n)" = "n") %>% dplyr::select(-Season)
  ntbd_pp_wtr <- ntbd_pp %>% filter(Season == "Winter") %>%  
    mutate(Analysis = "Predator-Prey") %>%
    rename("Winter (n)" = "n") %>% dplyr::select(-Season)
  ntbd_pp_sprg <- ntbd_pp %>% filter(Season == "Spring") %>%  
    mutate(Analysis = "Predator-Prey") %>%
    rename("Spring (n)" = "n") %>% dplyr::select(-Season)
  
  ntbd_ung_smr <- ntbd_ung %>% filter(Season == "Summer") %>% 
    mutate(Analysis = "Ungulate") %>%
    rename("Summer (n)" = "n") %>% dplyr::select(-Season)
  ntbd_ung_fall <- ntbd_ung %>% filter(Season == "Fall") %>%  
    mutate(Analysis = "Ungulate") %>%
    rename("Fall (n)" = "n") %>% dplyr::select(-Season)
  ntbd_ung_wtr <- ntbd_ung %>% filter(Season == "Winter") %>%  
    mutate(Analysis = "Ungulate") %>%
    rename("Winter (n)" = "n") %>% dplyr::select(-Season)
  ntbd_ung_sprg <- ntbd_ung %>% filter(Season == "Spring") %>%  
    mutate(Analysis = "Ungulate") %>%
    rename("Spring (n)" = "n") %>% dplyr::select(-Season)
  
  ntbd_con_smr <- ntbd_con %>% filter(Season == "Summer") %>%  
    mutate(Analysis = "Conspecific") %>%
    rename("Summer (n)" = "n") %>% dplyr::select(-Season)
  ntbd_con_fall <- ntbd_con %>% filter(Season == "Fall") %>%   
    mutate(Analysis = "Conspecific") %>%
    rename("Fall (n)" = "n") %>% dplyr::select(-Season)
  ntbd_con_wtr <- ntbd_con %>% filter(Season == "Winter") %>%   
    mutate(Analysis = "Conspecific") %>%
    rename("Winter (n)" = "n") %>% dplyr::select(-Season)
  ntbd_con_sprg <- ntbd_con %>% filter(Season == "Spring") %>%   
    mutate(Analysis = "Conspecific") %>%
    rename("Spring (n)" = "n") %>% dplyr::select(-Season)
  
  #'  Merge seasonal predator-prey and ungulate data together
  ntbd_smr <- rbind(ntbd_pp_smr, ntbd_ung_smr) #ntbd_con_smr 
  ntbd_fall <- rbind(ntbd_pp_fall, ntbd_ung_fall) #ntbd_con_fall 
  ntbd_wtr <- rbind(ntbd_pp_wtr, ntbd_ung_wtr) #ntbd_con_wtr 
  ntbd_sprg <- rbind(ntbd_pp_sprg, ntbd_ung_sprg) #ntbd_con_sprg 
  
  #'  Create table in wide format that combines all season data
  ntbd_wide <- full_join(ntbd_smr, ntbd_fall, by = c("Analysis", "Pair")) %>%
    full_join(ntbd_wtr, by = c("Analysis", "Pair")) %>%
    full_join(ntbd_sprg, by = c("Analysis", "Pair")) %>%
    relocate("Analysis", .after = "Pair") %>%
    arrange(Pair) %>%
    rename("Species pair" = "Pair") 

  write.csv(ntbd_wide, "./Outputs/summary_table_tbd_pairs.csv")
    
  
  #'  Summarize tbd 
  #'  Range of times between detections of predators and prey
  range(tbd_pp$tbd_min_round)
  #'  Mean number of minutes (SE) between detections of each species pair & season
  (tbd_mins_pp <- tbd_pp %>%
      group_by(Season, Pair) %>%
      summarise(mu = mean(tbd_min_round), se = (sd(tbd_min_round))/sqrt(length(tbd_min_round))) %>%
      ungroup())
  #'  Mean number of minutes (SE) between detections across all predator-prey pairs
  (tbd_mu_pp <- tbd_pp %>%
      summarise(mu = mean(tbd_min_round), se = (sd(tbd_min_round))/sqrt(length(tbd_min_round)))) 
  
  #'  Range of times between detections of different ungulates
  range(tbd_ung$tbd_min_round)
  #'  Mean number of minutes (SE) between detections of different ungulates for each species & season
  (tbd_mins_ung <- tbd_ung %>%
      group_by(Season, Pair) %>%
      summarise(mu = mean(tbd_min_round), se = (sd(tbd_min_round))/sqrt(length(tbd_min_round))) %>%
      ungroup())
  #'  Mean number of minutes (SE) between detections of different ungulates across all ungulate species
  (tbd_mu_ung <- tbd_ung %>%
      summarise(mu = mean(tbd_min_round), se = (sd(tbd_min_round))/sqrt(length(tbd_min_round)))) 
  
  #'  Range of times between detections of conspecifics
  range(tbd_con$tbd_min_round)
  #'  Mean number of minutes (SE) between detections of conspecifics for each species & season
  (tbd_mins_con <- tbd_con %>%
      group_by(Season, Pair) %>%
      summarise(mu = mean(tbd_min_round), se = (sd(tbd_min_round))/sqrt(length(tbd_min_round))) %>%
      ungroup())
  #'  Mean number of minutes (SE) between detections of conspecifics across all ungulate species
  (tbd_mu_con <- tbd_con %>%
      summarise(mu = mean(tbd_min_round), se = (sd(tbd_min_round))/sqrt(length(tbd_min_round)))) 
  
  
  
  
  
  
  