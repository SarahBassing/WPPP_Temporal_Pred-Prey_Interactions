  #'  =======================================
  #'  Time-between-detections metric
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing, Cameron Ho, Hunter Hicks
  #'  July 2022
  #'  =======================================
  #'  Calculate times-between-detections of predators followed by prey detected 
  #'  at the same camera site, of different ungulate species detected at the same
  #'  camera site, of individuals of the same species (conspecifics) detected at
  #'  the same camera site, and of prey followed by predators detected at the same
  #'  camera site. GLMMs will be fit to these measures of TBD to estimate the 
  #'  latency of site-use by ungulates in response to predation risk.
  #'  ================================
  
  #'  Load packages
  library(data.table)
  library(lubridate)
  library(chron)
  library(sp)
  library(raster)
  library(tidyverse)
  
  #'  Read in and format data
  megadata <- read.csv("./Data/full_camdata18-21_2022-08-19.csv") %>%
    dplyr::select("File", "DateTime", "Date", "Time", "CameraLocation", 
                  "Camera_Lat", "Camera_Long", "Animal", "Human", "Vehicle", 
                  "Species", "HumanActivity", "Count", "Monitoring") %>%
    filter(!grepl("Moultrie", CameraLocation)) %>%
    #  Need to have something in the Species column for each detection
    mutate(
      Species = ifelse(Human == "TRUE" | Human == "true", "Human", Species),
      Species = ifelse(Vehicle == "TRUE" | Vehicle == "true", "Vehicle", Species),
      Species = ifelse(Species == "", "NA", Species),
      HumanActivity = ifelse(HumanActivity == "", "NA", HumanActivity),
      Monitoring = ifelse(Monitoring == "Closed road", "Dirt road", Monitoring),
      Monitoring = ifelse(Monitoring == "Decommissioned road", "Dirt road", Monitoring),
      Monitoring = ifelse(Monitoring == "Game trail", "Trail", Monitoring)) %>%
    #  Remove rows where no detection occurred but snuck into this data set somehow
    filter(!(Animal == "FALSE" & Human == "FALSE" & Vehicle == "FALSE") | (Animal == "false" & Human == "false" & Vehicle == "false")) %>%
    #'  Remove observations that are still NA
    filter(!is.na(Species)) %>%
    #'  Format date and time data
    mutate(
      DateTime = as.POSIXct(DateTime,
                            format="%Y-%m-%d %H:%M:%S",tz="America/Los_Angeles"),
      Date = as.Date(Date, format = "%Y-%m-%d"),
      Time = chron(times = Time)
    ) %>%
    mutate(Category = ifelse(Species == "Bobcat" | Species == "Black Bear" | 
                               Species == "Cougar" | Species == "Coyote" | 
                               Species == "Lynx" | Species == "Wolf", "Predator", "Other"),
           Category = ifelse(Species == "Elk" | Species == "Mule Deer" | 
                               Species == "Moose" | Species == "White-tailed Deer", "Prey", Category)) 
  
  ungulate_mega <- megadata %>%
    #'  Remove other species - not needed and sometimes random birds or vehicles
    #'  show up in sequence of same species but breaks up the sequence so code 
    #'  below thinks they are independent events
    filter(Category != "Other")
  
  ####  Extract independent detections for wildlife species  ####
  #'  -------------------------------------------------------
  #'  Create a column identifying whether each image is an "independent" event
  #'  If camera site is diff from previous row then give unique value. If not then...
  #'  If species detected is diff from previous row at same site then give unique value. If not then...
  #'  If DateTime is >30 min from previous DateTime at same site for same species then give unique value. If not then...
  #'  Capture value is the same as that in the previous row.
  det_events <- function(dets) {
    dat <- arrange(dets, CameraLocation, DateTime)
    caps <- c()
    caps[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps[i] = i
      else (if (dat$Species[i-1] != dat$Species[i]) caps[i] = i
            else (if (difftime(dat$DateTime[i], dat$DateTime[i-1], units = c("mins")) > 30) caps[i] = i
                  else caps[i] = caps[i-1]))
    }
    
    caps <- as.factor(caps)
    
    #'  Add new column to larger data set
    capdata <- cbind(as.data.frame(dat), caps)
    
    #'  Filter data to the first image from each unique detection event
    firstprey <- capdata[capdata$Category == "Prey",] %>%
      group_by(caps) %>%
      slice(1L) %>%
      ungroup() %>%
      mutate(Det_type = "first")
    firstpred <- capdata[capdata$Category == "Predator",] %>%
      group_by(caps) %>%
      slice(1L) %>%
      ungroup() %>%
      mutate(Det_type = "first")
    
    #'  Filter data to the last image of each unique detection event
    lastpredator <- capdata[capdata$Category == "Predator",] %>% 
      group_by(caps) %>% 
      slice_tail() %>%
      ungroup() %>%
      mutate(Det_type = "last")
    lastprey <- capdata[capdata$Category == "Prey",] %>%
      group_by(caps) %>%
      slice_tail() %>%
      ungroup() %>%
      mutate(Det_type = "last")
    lastother <- capdata[capdata$Category == "Other",] %>% 
      group_by(caps) %>% 
      slice_tail() %>%
      ungroup() %>%
      mutate(Det_type = "last")
    
    data_list <- list(firstprey, firstpred, lastpredator, lastprey, lastother)
    names(data_list) <- c("firstprey", "firstpred", "lastpredator", "lastprey", "lastother")
    return(data_list)
  }
  predprey_caps <- det_events(megadata)
  prey_caps <- det_events(ungulate_mega)

  
  #'  Merge data based on last image of each predator/other detection and first 
  #'  image of each prey detection
  lastPred_firstPrey <- rbind(predprey_caps$lastpredator, predprey_caps$firstprey, predprey_caps$lastother)
  #'  Merge data based on last image of each prey/other detection and first image
  #'  of each predatory detection
  lastPrey_firstPred <- rbind(predprey_caps$lastprey, predprey_caps$firstpred, predprey_caps$lastother)
  #'  Merge last and first images of each ungulate detection
  back2back_prey <- rbind(prey_caps$lastprey, prey_caps$firstprey) %>%
    arrange(CameraLocation, DateTime, caps, Det_type)
  #'  Merge all together to be used for last ungulate vs first ungulate of diff spp
  lastUng_firstUng <- rbind(predprey_caps$lastprey, predprey_caps$lastpredator, predprey_caps$firstprey, predprey_caps$lastother)
  
  #'  Detections with only one image get duplicated b/c its the first & last image
  #'  Identify duplicate images (ignoring last column where they differ) and filter
  #'  to just one image per detection event
  dups <- back2back_prey %>%
    group_by_at(vars(-Det_type)) %>% 
    filter(n() > 1) %>%
    filter(Det_type == "last")
  #'  Remove the duplicate images from the larger data set
  back2back_prey <- anti_join(back2back_prey, dups)
  
  #'  Merge and filter to specific date ranges of interest
  #'  Want to include detections of ALL specie and humans to account for non-focal
  #'  species that are detected between species of interest, but doesn't matter if
  #'  its the first or last image of their unique detection events here
  #'  Filter data to specific date ranges
  seasonal_filter <- function(dets, yr) {
    summer <- dets %>%
      filter(Date > paste0(yr, "-05-31")) %>%
      filter(Date < paste0(yr, "-10-01")) %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "Category", "caps", "Det_type") %>%
      arrange(CameraLocation, DateTime)
    fall <- dets %>%
      filter(Date > paste0(yr, "-09-30")) %>%
      filter(Date < paste0(yr, "-12-01")) %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "Category", "caps", "Det_type") %>%
      arrange(CameraLocation, DateTime)
    winter <- dets %>%
      filter(Date > paste0(yr, "-11-30")) %>%
      filter(Date < paste0(yr+1, "-04-01")) %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "Category", "caps", "Det_type") %>%
      arrange(CameraLocation, DateTime)
    spring <- dets %>%
      filter(Date > paste0(yr+1, "-03-29")) %>%
      filter(Date < paste0(yr+1, "-06-01")) %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species", "Category", "caps", "Det_type") %>%
      arrange(CameraLocation, DateTime)
    det_list <- list(summer, fall, winter, spring)
    return(det_list)
  }
  #'  Create list of seasonal detection data per year
  lastfirst18 <- seasonal_filter(lastPred_firstPrey, yr = 2018)
  lastfirst19 <- seasonal_filter(lastPred_firstPrey, yr = 2019)
  lastfirst20 <- seasonal_filter(lastPred_firstPrey, yr = 2020)
  
  lastfirstpreypred18 <- seasonal_filter(lastPrey_firstPred, yr = 2018)
  lastfirstpreypred19 <- seasonal_filter(lastPrey_firstPred, yr = 2019)
  lastfirstpreypred20 <- seasonal_filter(lastPrey_firstPred, yr = 2020)
  
  lastfirstung18 <- seasonal_filter(lastUng_firstUng, yr = 2018)
  lastfirstung19 <- seasonal_filter(lastUng_firstUng, yr = 2019)
  lastfirstung20 <- seasonal_filter(lastUng_firstUng, yr = 2020)
  
  back2back18 <- seasonal_filter(back2back_prey, yr = 2018)
  back2back19 <- seasonal_filter(back2back_prey, yr = 2019)
  back2back20 <- seasonal_filter(back2back_prey, yr = 2020)
  
  #'  Join all detection data per season, across years
  #'  Last image of predator followed by first image of prey
  lpfprey_smr <- rbind(lastfirst18[[1]], lastfirst19[[1]], lastfirst20[[1]])
  lpfprey_fall <- rbind(lastfirst18[[2]], lastfirst19[[2]], lastfirst20[[2]])
  lpfprey_wtr <- rbind(lastfirst18[[3]], lastfirst19[[3]], lastfirst20[[3]])
  lpfprey_sprg <- rbind(lastfirst18[[4]], lastfirst19[[4]], lastfirst20[[4]])
  
  #'  Last image of prey followed by first image of predator
  lpfpred_smr <- rbind(lastfirstpreypred18[[1]], lastfirstpreypred19[[1]], lastfirstpreypred20[[1]])
  lpfpred_fall <- rbind(lastfirstpreypred18[[2]], lastfirstpreypred19[[2]], lastfirstpreypred20[[2]])
  lpfpred_wtr <- rbind(lastfirstpreypred18[[3]], lastfirstpreypred19[[3]], lastfirstpreypred20[[3]])
  lpfpred_sprg <- rbind(lastfirstpreypred18[[4]], lastfirstpreypred19[[4]], lastfirstpreypred20[[4]])
  
  #'  Last image of one ungulate species followed by first image of different ungulate species
  lufu_smr <- rbind(lastfirstung18[[1]], lastfirstung19[[1]], lastfirstung20[[1]])
  lufu_fall <- rbind(lastfirstung18[[2]], lastfirstung19[[2]], lastfirstung20[[2]])
  lufu_wtr <- rbind(lastfirstung18[[3]], lastfirstung19[[3]], lastfirstung20[[3]])
  lufu_sprg <- rbind(lastfirstung18[[4]], lastfirstung19[[4]], lastfirstung20[[4]])
  
  #'  Last image of ungulate in an independent detection followed by first image 
  #'  of same ungulate species in a next independent detection
  b2b_smr <- rbind(back2back18[[1]], back2back19[[1]], back2back20[[1]])
  b2b_fall <- rbind(back2back18[[2]], back2back19[[2]], back2back20[[2]])
  b2b_wtr <- rbind(back2back18[[3]], back2back19[[3]], back2back20[[3]])
  b2b_sprg <- rbind(back2back18[[4]], back2back19[[4]], back2back20[[4]])
  
  ####  Filter predator-prey detection data  ####
  #'  ----------------------------------------
  #'  Group multiple detection events of same category (but of different species) 
  #'  when they occur sequentially, then reduce to a single observation (e.g., 
  #'  we only care about the LAST of the last predator detections in a series of 
  #'  predator detections).
  thin_dat <- function(dets) {
    dat <- arrange(dets, CameraLocation, DateTime) %>%
      dplyr::select(-c(caps, Det_type))
    caps_new <- c()
    caps_new[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps_new[i] = i
      else(if (dat$Category[i-1] != dat$Category[i]) caps_new[i] = i
           else caps_new[i] = caps_new[i-1])
    }
    
    caps_new <- as.factor(caps_new)
    
    #'  Add new column to larger data set
    capdata <- cbind(as.data.frame(dat), caps_new)
    
    #'  Remove all extra detections when multiple detections of same category occur in a row
    firstpreyspp <- capdata[capdata$Category == "Prey",] %>%
      group_by(caps_new) %>% 
      slice(1L) %>%
      ungroup()
    lasteverythingelse <- capdata[capdata$Category != "Prey",] %>%
      group_by(caps_new) %>% 
      slice_tail() %>%
      ungroup()
    #'  Combine into final data set
    dets <- rbind(firstpreyspp, lasteverythingelse) %>%
      arrange(CameraLocation, DateTime)
    return(dets)
  }
  #'  Last predator followed by first prey images by season
  lpfprey_smr_thin <- thin_dat(lpfprey_smr)
  lpfprey_fall_thin <- thin_dat(lpfprey_fall)
  lpfprey_wtr_thin <- thin_dat(lpfprey_wtr)
  lpfprey_sprg_thin <- thin_dat(lpfprey_sprg)
  
  #'  For first predator detection
  thin_dat <- function(dets) {
    dat <- arrange(dets, CameraLocation, DateTime) %>%
      dplyr::select(-c(caps, Det_type))
    caps_new <- c()
    caps_new[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps_new[i] = i
      else(if (dat$Category[i-1] != dat$Category[i]) caps_new[i] = i
           else caps_new[i] = caps_new[i-1])
    }
    
    caps_new <- as.factor(caps_new)
    
    #'  Add new column to larger data set
    capdata <- cbind(as.data.frame(dat), caps_new)
    
    #'  Remove all extra detections when multiple detections of same category occur in a row
    firstpredspp <- capdata[capdata$Category == "Predator",] %>%
      group_by(caps_new) %>% 
      slice(1L) %>%
      ungroup()
    lasteverythingelse <- capdata[capdata$Category != "Predator",] %>%
      group_by(caps_new) %>% 
      slice_tail() %>%
      ungroup()
    #'  Combine into final data set
    dets <- rbind(firstpredspp, lasteverythingelse) %>%
      arrange(CameraLocation, DateTime)
    return(dets)
  }
  #'  Last prey followed by first predator images by season
  lpfpred_smr_thin <- thin_dat(lpfpred_smr)
  lpfpred_fall_thin <- thin_dat(lpfpred_fall)
  lpfpred_wtr_thin <- thin_dat(lpfpred_wtr)
  lpfpred_sprg_thin <- thin_dat(lpfpred_sprg)
  
  #'  Function to reduce detections to just series of a spp1 detection followed
  #'  by a spp2 detection (e.g., predator then prey, but no other species detected 
  #'  in between) at each camera site
  spppair_dat <- function(dets, spp1, spp2) {
    #'  Assign same ID to all detection events from the same camera
    dat <- arrange(dets, CameraLocation, DateTime)
    cam <- c()
    cam[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) cam[i] = i
      else cam[i] = cam[i-1]
    }
    
    #'  Identify images where a spp1 is detected right before a spp2
    spp_caps1 <- c()
    spp_caps1[1] <- "N"
    for (i in 2:nrow(dat)){
      if (dat$Category[i-1] == spp1 & dat$Category[i] == spp2) spp_caps1[i-1] = "Y"
      else spp_caps1[i-1] = "N"
      #'  Flag situations where last image at a camera site is of a predator - 
      #'  this confuses time-btwn-detection calculations because script ignores 
      #'  changes in location when looping through each row
      if (spp_caps1[i-1] == "Y" & dat$CameraLocation[i-1] != dat$CameraLocation[i]) spp_caps1[i-1] = "FLAG"
    }
    #'  Add "N" to very end of vector so the length matches number of rows in dets
    spp_caps1[nrow(dat)] <- "N"
    
    #'  Identify images where a spp2 is detected right after a spp1
    spp_caps2 <- c()
    spp_caps2[1] <- "N"
    for (i in 2:nrow(dat)){
      if (dat$Category[i-1] == spp1 & dat$Category[i] == spp2) spp_caps2[i] = "Y"
      else spp_caps2[i] = "N"
      #'  Flag situations where first image at a camera is of a prey - this 
      #'  confuses time-btwn-detection calculations because script ignores changes 
      #'  in location when looping through each row
      if (spp_caps2[i] == "Y" & dat$CameraLocation[i] != dat$CameraLocation[i-1]) spp_caps2[i] = "FLAG"
    }

    #'  Identify the species pair for easier filtering later on
    spp1spp2 <- c()
    spp1spp2[1] <- NA
    for(i in 2:nrow(dat)){
      if (dat$Category[i-1] == spp1 & dat$Category[i] == spp2) spp1spp2[i] = paste0(dat$Species[i-1], "_", dat$Species[i])
      else spp1spp2[i] =  NA
    }
    
    #'  Add new columns to larger data set and filter
    spp_new1 <- as.factor(spp_caps1)
    spp_new2 <- as.factor(spp_caps2)
    spp_pair <- as.factor(spp1spp2)
    capdata <- cbind(as.data.frame(dat), cam, spp_new1, spp_new2, spp_pair) %>%
      #'  Filter detection events to just situations where spp2 follows spp1 at
      #'  the same camera site (drop any "N"s and "FLAG" rows)
      filter(spp_new1 == "Y" | spp_new2 == "Y") %>%
      #'  Add columns designating predator hunting modes (ambush vs coursing) &  
      #'  trophic level of each predator (apex vs meso). Note, considering
      #'  black bears to be ambush predators for now because can't really run
      #'  down prey over long distances.
      mutate(HuntingMode = ifelse(Species == "Coyote" | Species == "Wolf", "Coursing", "Ambush"),
             HuntingMode = ifelse(grepl("Bobcat", spp_pair), "Ambush", HuntingMode),
             HuntingMode = ifelse(grepl("Cougar", spp_pair), "Ambush", HuntingMode),
             HuntingMode = ifelse(grepl("Lynx", spp_pair), "Ambush", HuntingMode),
             TrophicLevel = ifelse(Species == "Black Bear" | Species == "Cougar" | Species == "Wolf", "Apex", "Meso"),
             TrophicLevel = ifelse(grepl("Black Bear", spp_pair), "Apex", TrophicLevel),
             TrophicLevel = ifelse(grepl("Cougar", spp_pair), "Apex", TrophicLevel),
             TrophicLevel = ifelse(grepl("Wolf", spp_pair), "Apex", TrophicLevel),
             PredatorID = gsub("_.*","", spp_pair))
    return(capdata)
  }
  resp2pred_smr <- spppair_dat(lpfprey_smr_thin, spp1 = "Predator", spp2 = "Prey")
  resp2pred_fall <- spppair_dat(lpfprey_fall_thin, spp1 = "Predator", spp2 = "Prey")
  resp2pred_wtr <- spppair_dat(lpfprey_wtr_thin, spp1 = "Predator", spp2 = "Prey")
  resp2pred_sprg <- spppair_dat(lpfprey_sprg_thin, spp1 = "Predator", spp2 = "Prey")
  
  resp2prey_smr <- spppair_dat(lpfpred_smr_thin, spp1 = "Prey", spp2 = "Predator")
  resp2prey_fall <- spppair_dat(lpfpred_fall_thin, spp1 = "Prey", spp2 = "Predator")
  resp2prey_wtr <- spppair_dat(lpfpred_wtr_thin, spp1 = "Prey", spp2 = "Predator")
  resp2prey_sprg <- spppair_dat(lpfpred_sprg_thin, spp1 = "Prey", spp2 = "Predator")
  
  ####  Calculate times between detection events of predator-prey  ####
  #'  -------------------------------------------------------------
  #'  Function to calculate time between detection events of two focal species
  #'  Data structured so only last image of spp1 and first image of spp2 per
  #'  detection event are included in data frame.
  tbd <- function(detection_data, spp1, unittime) {
    #'  Create empty vector to be filled
    detection_data$TimeSinceLastDet <- c()
    #'  Fill first element of the vector to get it started
    detection_data$TimeSinceLastDet[1] <- 0
    #'  Loop through each row to calculate elapsed time since previous detection
    for (i in 2:nrow(detection_data)){
      #'  If previous detection was spp1, set time to 0
      if (detection_data$Category[i-1] == spp1) detection_data$TimeSinceLastDet[i] = 0
      #'  If current detection is spp2 and follows detection of spp1, calculate
      #'  the difference in time from previous detection to current detection
      if (detection_data$Category[i] != spp1) detection_data$TimeSinceLastDet[i] = difftime(detection_data$DateTime[i], detection_data$DateTime[i-1], units = unittime)
    }
    #'  Retain only prey observations (don't need the actual predator detections)
    # detection_data <- filter(detection_data, Category == "Prey")
    return(detection_data)
  }
  #'  Calculate time between detections for different pairs of species of interest
  #'  spp1 should be the species detected first, unittime is the unit of time 
  #'  to make calculations in (options are: "sec", "min", "hour", "day")
  #'  Note: there should be NO negative values! If there are negative values this
  #'  means the script is calculating times between detections across camera sites
  tbd_pred.prey_smr <- tbd(resp2pred_smr, spp1 = "Predator", unittime = "min")
  tbd_pred.prey_fall <- tbd(resp2pred_fall, spp1 = "Predator", unittime = "min")
  tbd_pred.prey_wtr <- tbd(resp2pred_wtr, spp1 = "Predator", unittime = "min")
  tbd_pred.prey_sprg <- tbd(resp2pred_sprg, spp1 = "Predator", unittime = "min")
  
  #'  Calculate time between detections for different pairings of prey species
  #'  followed by predator species
  tbd_prey.pred_smr <- tbd(resp2prey_smr, spp1 = "Prey", unittime = "min")
  tbd_prey.pred_fall <- tbd(resp2prey_fall, spp1 = "Prey", unittime = "min")
  tbd_prey.pred_wtr <- tbd(resp2prey_wtr, spp1 = "Prey", unittime = "min")
  tbd_prey.pred_sprg <- tbd(resp2prey_sprg, spp1 = "Prey", unittime = "min")
  
  
  ####  Calculate times between detection events of different ungulate species  ####
  #'  --------------------------------------------------------------------------
  #'  Function to calculate time between detection events of ungulate species
  #'  Data structured so only first and last image per detection event are included.
  tbd_ungulate <- function(dat, unittime) {
    detection_data <- arrange(dat, CameraLocation, DateTime, File)
    
    #'  Identify the species pair for easier filtering later on
    spp1spp2 <- c()
    spp1spp2[1] <- NA
    for(i in 2:nrow(detection_data)){
      if (detection_data$Species[i-1] != detection_data$Species[i]) spp1spp2[i] = paste0(detection_data$Species[i-1], "_", detection_data$Species[i])
      else spp1spp2[i] =  NA
    }
    spp_pair <- as.factor(spp1spp2)
    detection_data <- cbind(detection_data, spp_pair) %>%
      mutate(PreviousDet = gsub("_.*","", spp_pair))
    
    #'  Create empty vector to be filled
    detection_data$TimeSinceLastDet <- c()
    #'  Fill first element of the vector to get it started
    detection_data$TimeSinceLastDet[1] <- 0
    #'  Loop through each row to calculate elapsed time since previous detection
    for (i in 2:nrow(detection_data)){
      #'  If previous detection was different species, calculate
      #'  the difference in time from previous detection to current detection
      if (detection_data$Species[i-1] != detection_data$Species[i]) detection_data$TimeSinceLastDet[i] = difftime(detection_data$DateTime[i], detection_data$DateTime[i-1], units = unittime)
      #'  If previous detection is same species as current detection, set to 0 
      if (detection_data$Species[i-1] == detection_data$Species[i]) detection_data$TimeSinceLastDet[i] = 0
      #'  If previous detection is not an ungulate, set to 0 
      if (detection_data$Category[i-1] != "Prey") detection_data$TimeSinceLastDet[i] = 0
      #'  If current detection is not an ungulate, set to 0 
      if (detection_data$Category[i] != "Prey") detection_data$TimeSinceLastDet[i] = 0
      #'  If previous detection was from a different camera site, set to 0
      if (detection_data$CameraLocation[i-1] != detection_data$CameraLocation[i]) detection_data$TimeSinceLastDet[i] = 0
    }
    
    #'  Filter to just back-to-back ungulate detections
    ungulates_only <- detection_data %>%
      filter(Category == "Prey") %>%
      filter(TimeSinceLastDet > 0) %>%
      dplyr::select(-c(caps, Det_type))
    return(ungulates_only)
  }
  #'  Calculate time between detections for different pairs of ungulate species.
  #'  unittime is unit of time to make calculations in ("sec", "min", "hour", "day")
  #'  Note: there should be NO negative values! If there are negative values this
  #'  means the script is calculating times between detections across camera sites
  tbd_lufu_smr <- tbd_ungulate(lufu_smr, unittime = "min")
  tbd_lufu_fall <- tbd_ungulate(lufu_fall, unittime = "min")
  tbd_lufu_wtr <- tbd_ungulate(lufu_wtr, unittime = "min")
  tbd_lufu_sprg <- tbd_ungulate(lufu_sprg, unittime = "min")
  
  ####  Calculate times between detection events of conspecifics  ####
  #'  ------------------------------------------------------------
  #'  Function to calculate time between detection events of conspecifics
  #'  Data structured so only first and last image per detection event are included.
  tbd <- function(dat, unittime) {
    #'  Group multiple detection events of same species when they occur sequentially
    #'  at same camera site
    dat <- arrange(dat, CameraLocation, DateTime, File)
    caps_new <- c()
    caps_new[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps_new[i] = i
      else(if (dat$Species[i-1] != dat$Species[i]) caps_new[i] = i
           else caps_new[i] = caps_new[i-1])
    }
    #'  Add new column to larger data set
    caps_new <- as.factor(caps_new)
    detection_data <- cbind(as.data.frame(dat), caps_new)
    detection_data <- arrange(detection_data, CameraLocation, DateTime, caps, Det_type)
    
    #'  Create empty vector to be filled
    detection_data$TimeSinceLastDet <- c()
    #'  Fill first element of the vector to get it started
    detection_data$TimeSinceLastDet[1] <- 0
    #'  Loop through each row to calculate elapsed time since previous detection
    for (i in 2:nrow(detection_data)){
      #'  If previous detection was different species, set time to 0
      if (detection_data$Species[i-1] != detection_data$Species[i]) detection_data$TimeSinceLastDet[i] = 0
      #'  If previous detection is same species as current detection, calculate
      #'  the difference in time from previous detection to current detection
      if (detection_data$Species[i-1] == detection_data$Species[i]) detection_data$TimeSinceLastDet[i] = difftime(detection_data$DateTime[i], detection_data$DateTime[i-1], units = unittime)
      #'  If previous detection was from a different camera site, set to 0
      if (detection_data$CameraLocation[i-1] != detection_data$CameraLocation[i]) detection_data$TimeSinceLastDet[i] = 0
    }
    #'  Retain only Det_type = first from each cap event (don't want the tbd from
    #'  the first to last image of the same detection event, only the tbd from
    #'  the last image of a detection to the first image of the next detection)
    detection_data <- filter(detection_data, Det_type == "first") %>%
      #'  Filter out times = 0 (first detection at a camera site) - smallest length
      #'  of time possible btwn detections is 30min due to definition of "independent
      #'  detection events" for conspecificss
      filter(TimeSinceLastDet > 0) %>%
      dplyr::select(-c(caps, Det_type, caps_new))
    return(detection_data)
  }
  #'  Calculate time between detections for different pairs of species of interest
  #'  spp1 should be the species detected first, unittime is the unit of time 
  #'  to make calculations in (options are: "sec", "min", "hour", "day")
  #'  Note: there should be NO negative values! If there are negative values this
  #'  means the script is calculating times between detections across camera sites
  tbd_conspif_smr <- tbd(b2b_smr, unittime = "min")
  tbd_conspif_fall <- tbd(b2b_fall, unittime = "min")
  tbd_conspif_wtr <- tbd(b2b_wtr, unittime = "min")
  tbd_conspif_sprg <- tbd(b2b_sprg, unittime = "min")
  
  ####  Add covariate data to tbd data  ####
  #'  ----------------------------------
  #'  Add habitat complexity and other site-level covariates to each observation
  #'  FYI: Complexity_index1 = TRI_250m * (PercForest*100)
  #'  "Low" HCI = Complexity_index1 values < mean(Complexity_index1)
  #'  "High" HCI = Complexity_index1 values => mean(Complexity_index1)
  stations_data <- read.csv("./Data/cam_stations_hab_complex_data.csv") %>%
    dplyr::select(-"X")
  
  #'  Join time-between-detection data with site-level covariates
  tbd_pred.prey_smr <- left_join(tbd_pred.prey_smr, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Summer")
  tbd_pred.prey_fall <- left_join(tbd_pred.prey_fall, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Fall")
  tbd_pred.prey_wtr <- left_join(tbd_pred.prey_wtr, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Winter")
  tbd_pred.prey_sprg <- left_join(tbd_pred.prey_sprg, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Spring")
  
  tbd_prey.pred_smr <- left_join(tbd_prey.pred_smr, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Summer")
  tbd_prey.pred_fall <- left_join(tbd_prey.pred_fall, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Fall")
  tbd_prey.pred_wtr <- left_join(tbd_prey.pred_wtr, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Winter")
  tbd_prey.pred_sprg <- left_join(tbd_prey.pred_sprg, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Spring")
  
  tbd_ungulate_smr <- left_join(tbd_lufu_smr, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Summer")
  tbd_ungulate_fall <- left_join(tbd_lufu_fall, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Fall")
  tbd_ungulate_wtr <- left_join(tbd_lufu_wtr, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Winter")
  tbd_ungulate_sprg <- left_join(tbd_lufu_sprg, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Spring")
  
  tbd_conspif_smr <- left_join(tbd_conspif_smr, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Summer")
  tbd_conspif_fall <- left_join(tbd_conspif_fall, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Fall")
  tbd_conspif_wtr <- left_join(tbd_conspif_wtr, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Winter")
  tbd_conspif_sprg <- left_join(tbd_conspif_sprg, stations_data, by = "CameraLocation") %>%
    mutate(Season = "Spring")
  
  #'  Merge into one large data set
  tbd_pred.prey <- rbind(tbd_pred.prey_smr, tbd_pred.prey_fall, tbd_pred.prey_wtr, tbd_pred.prey_sprg) %>%
    arrange(CameraLocation, DateTime) %>%
    dplyr::select(-c(File, Category, caps_new, cam, spp_new1, spp_new2)) %>%
    relocate(Year, .before = DateTime) %>%
    relocate(Study_Area, .before = Year) %>%
    relocate(Season, .before = DateTime) %>%
    relocate(TimeSinceLastDet, .after = backgroundRisk_For)
  
  tbd_prey.pred <- rbind(tbd_prey.pred_smr, tbd_prey.pred_fall, tbd_prey.pred_wtr, tbd_prey.pred_sprg) %>%
    arrange(CameraLocation, DateTime) %>%
    dplyr::select(-c(File, Category, caps_new, cam, spp_new1, spp_new2)) %>%
    relocate(Year, .before = DateTime) %>%
    relocate(Study_Area, .before = Year) %>%
    relocate(Season, .before = DateTime) %>%
    relocate(TimeSinceLastDet, .after = backgroundRisk_For)
  
  tbd_ungulate <- rbind(tbd_ungulate_smr, tbd_ungulate_fall, tbd_ungulate_wtr, tbd_ungulate_sprg) %>%
    arrange(CameraLocation, DateTime) %>%
    dplyr::select(-c(File, Category)) %>%
    relocate(Year, .before = DateTime) %>%
    relocate(Study_Area, .before = Year) %>%
    relocate(Season, .before = DateTime) %>%
    relocate(TimeSinceLastDet, .after = backgroundRisk_For)
  
  tbd_conspif <- rbind(tbd_conspif_smr, tbd_conspif_fall, tbd_conspif_wtr, tbd_conspif_sprg) %>%
    arrange(CameraLocation, DateTime) %>%
    dplyr::select(-c(File, Category)) %>%
    relocate(Year, .before = DateTime) %>%
    relocate(Study_Area, .before = Year) %>%
    relocate(Season, .before = DateTime) %>%
    relocate(TimeSinceLastDet, .after = backgroundRisk_For)
  
  #'  Double check there are no negative times-between-detection
  summary(tbd_pred.prey$TimeSinceLastDet)
  summary(tbd_prey.pred$TimeSinceLastDet)
  summary(tbd_ungulate$TimeSinceLastDet)
  summary(tbd_conspif$TimeSinceLastDet)
  
  #'  Conspecific data there should be NO tbd < 30min because 30min minimum was
  #'  required to identify independent detection events of the same species.
  #'  The few instances where tbd < 30min in this data set arise b/c another spp
  #'  (usually coyote) shows up within sequential ungulate images. Without 2nd
  #'  spp being detected in the mix all images would be considered 1 detection
  #'  event. Removing these instances b/c they are wrong based on definition of
  #'  independent detection event & because they belong (and are) in the pred-prey
  #'  tbd data set, not the conspecific data set.
  tbd_conspif <- filter(tbd_conspif, TimeSinceLastDet >= 30.00) #5.00 for 5-min interval
  
  #'  SAVE!
  write.csv(tbd_pred.prey, file = paste0("./Outputs/tbd_pred.prey_", Sys.Date(), ".csv"))
  save(tbd_pred.prey, file = paste0("./Outputs/tbd_pred.prey_", Sys.Date(), ".RData"))
  
  write.csv(tbd_prey.pred, file = paste0("./Outputs/tbd_prey.pred_", Sys.Date(), ".csv"))
  save(tbd_prey.pred, file = paste0("./Outputs/tbd_prey.pred_", Sys.Date(), ".RData"))
  
  write.csv(tbd_ungulate, file = paste0("./Outputs/tbd_ungulate_", Sys.Date(), ".csv"))
  save(tbd_ungulate, file = paste0("./Outputs/tbd_ungulate_", Sys.Date(), ".RData"))
  
  write.csv(tbd_conspif, file = paste0("./Outputs/tbd_conspif_", Sys.Date(), ".csv"))
  save(tbd_conspif, file = paste0("./Outputs/tbd_conspif_", Sys.Date(), ".RData"))

  
  
  #'  Next up: ExponentialGLMM_Latency_Analysis scripts (the exact script depends 
  #'  on which data set is being analyzed - Pred-Prey, Prey-Pred, Ungulate, Conspecific)
  
  
    
  
  
  