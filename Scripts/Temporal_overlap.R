  #'  ==================================================================
  #'  Temporal overlap patterns in response to background predation risk
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  June 2022
  #'  ==================================================================
  #'  Script to estimate temporal overlap of prey species in repsonse to varying
  #'  levels of background predation risk in a multi-predator system. 
  #'  ==================================================================
  
  #'  Load libraries
  library(lubridate)
  library(chron)
  library(overlap)
  library(circular)
  library(ggplot2)
  library(khroma)
  library(patchwork)
  library(sp)
  library(raster)
  library(tidyverse)
  
  #'  Load and format detection data
  megadata <- read.csv("./Data/full_camdata18-21_2022-08-19.csv") %>%  
    dplyr::select("File", "DateTime", "Date", "Time", "CameraLocation", 
                  "Camera_Lat", "Camera_Long", "Animal", "Human", "Vehicle", 
                  "Species", "HumanActivity", "Count", "Land_Mgnt", "Land_Owner") %>%
    filter(!grepl("Moultrie", CameraLocation)) %>%
    #  Need to have something in the Species column for each detection
    mutate(
      Species = ifelse(Human == "TRUE" | Human == "true", "Human", Species),
      Species = ifelse(Vehicle == "TRUE" | Vehicle == "true", "Vehicle", Species),
      Species = ifelse(Species == "", "NA", Species),
      HumanActivity = ifelse(HumanActivity == "", "NA", HumanActivity),
      #'  Identify if cameras were on public (1) or private (0) land
      Public1 = ifelse(Land_Mgnt != "Private", 1, 0),
      #'  Even though timberland is private, it's generally open to public recreation
      #'  so considering it public for these purposes
      Public1 = ifelse(Land_Mgnt == "Private" & Land_Owner == "Private timber", 1, Public1)
    ) %>%
    #  Remove rows where no detection occurred but snuck into this data set somehow
    filter(!(Animal == "FALSE" & Human == "FALSE" & Vehicle == "FALSE") | (Animal == "false" & Human == "false" & Vehicle == "false")) %>%
    #'  Remove observations of humans
    filter(!is.na(Species)) %>%
    mutate(
      DateTime = as.POSIXct(DateTime,
                            format="%Y-%m-%d %H:%M:%S",tz="America/Los_Angeles"),
      Date = as.Date(Date, format = "%Y-%m-%d"),
      Time = chron(times = Time),
      radTime = ((hour(DateTime)*60 + minute(DateTime) + second(DateTime)/(60))/1440)*2*pi,
    )
  stations_data <- read.csv("./Data/cam_stations_data.csv") %>%
    dplyr::select(c("CameraLocation", "Year", "Study_Area", "Elev", "Slope", "TRI", 
                    "PercForest", "PercForestMix2", "Canopy_Cov", "Latitude", "Longitude")) %>%
    mutate(Canopy_Cov = ifelse(is.na(Canopy_Cov), 22, Canopy_Cov)) # fill in missing values based on mean canopy cover (22%)
        
  #'  Make camera coordinates spatial
  cam_locs <- SpatialPoints(cbind(stations_data$Longitude, stations_data$Latitude), 
                      proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  #'  Read in mean TRI raster and extract values at each camera site
  mean_tri <- raster("./Shapefiles/mean_tri_250m.tif")
  cam_tri <- raster::extract(mean_tri, cam_locs)
  stations_data$TRI_250m <- cam_tri
  
  #'  Identity which predators were detected at each camera site (e.g., predation risk)
  bear_dets <- stations_data$CameraLocation %in% megadata$CameraLocation[megadata$Species == "Black Bear"]
  bear_det <- as.data.frame(bear_dets) %>% transmute(bear_det = ifelse(bear_dets == TRUE, 1, 0))
  bob_dets <- stations_data$CameraLocation %in% megadata$CameraLocation[megadata$Species == "Bobcat"]
  bob_det <- as.data.frame(bob_dets) %>% transmute(bob_det = ifelse(bob_dets == TRUE, 1, 0))
  coug_dets <- stations_data$CameraLocation %in% megadata$CameraLocation[megadata$Species == "Cougar"]
  coug_det <- as.data.frame(coug_dets) %>% transmute(coug_det = ifelse(coug_dets == TRUE, 1, 0))
  coy_dets <- stations_data$CameraLocation %in% megadata$CameraLocation[megadata$Species == "Coyote"]
  coy_det <- as.data.frame(coy_dets) %>% transmute(coy_det = ifelse(coy_dets == TRUE, 1, 0))
  lynx_dets <- stations_data$CameraLocation %in% megadata$CameraLocation[megadata$Species == "Lynx"]
  lynx_det <- as.data.frame(lynx_dets) %>% transmute(lynx_det = ifelse(lynx_dets == TRUE, 1, 0))
  wolf_dets <- stations_data$CameraLocation %in% megadata$CameraLocation[megadata$Species == "Wolf"]
  wolf_det <- as.data.frame(wolf_dets) %>% transmute(wolf_det = ifelse(wolf_dets == TRUE, 1, 0))
  
  #'  Add covariates representing background predation risk to
  stations_data <- stations_data %>%
    #'  Create habitat complexity index by multiplying TRI by forest habitat covs
    #'  index 1) 250m spatial scale; index 2) 30m/camera site spatial scale
    mutate(Complexity_index1 = TRI_250m * (PercForest*100), # multiply % forest (in decimals) by 100 so on same scale as canopy cover
           Complexity_index2 = TRI * Canopy_Cov)
  #'  Add naive occupancy for each predator and tally number of predators spp detected per site
  stations_data <- cbind(stations_data, bear_det, bob_det, coug_det, coy_det, lynx_det, wolf_det) %>%
    rowwise() %>%mutate(total_preds = sum(c_across(15:20)))
  
  #'  Summary stats on background predation risk
  summary(stations_data)
  hist(stations_data$Complexity_index1, breaks = 20)
  abline(v = median(stations_data$Complexity_index1), col = "blue")
  abline(v = mean(stations_data$Complexity_index1), col = "red")
  hist(stations_data$Complexity_index2, breaks = 20)
  abline(v = median(stations_data$Complexity_index2), col = "blue")
  abline(v = mean(stations_data$Complexity_index2), col = "red")
  #'  Using the mean as a cutoff to differentiate high vs low predation risk:
  #'  low background risk < mean complexity index; high background risk >= mean complexity risk
  
  #'  Add categorical variable designating level of background risk
  stations_data <- mutate(stations_data, backgroundRisk = ifelse(Complexity_index1 < mean(stations_data$Complexity_index1), "Low", "High"))
  
  #'  Is naive predator occupancy associated with habitat complexity?
  bear_index_cor <- glm(bear_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(bear_index_cor)
  bob_index_cor <- glm(bob_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(bob_index_cor)
  coug_index_cor <- glm(coug_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(coug_index_cor)
  coy_index_cor <- glm(coy_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(coy_index_cor)
  lynx_index_cor <- glm(lynx_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(lynx_index_cor)
  wolf_index_cor <- glm(wolf_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(wolf_index_cor)
  
  #'  Join detection data with camera station data
  megadata <- right_join(megadata, stations_data, by = "CameraLocation")
  #'  Make each observation spatial so times can be adjusted by suntime at their precise location
  xy <- SpatialPoints(cbind(megadata$Camera_Long, megadata$Camera_Lat),
                      proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  #'  Convert radians to sun times to account for seasonal changes in sunlight 
  #'  Converts a vector of clock times to "sun times", by mapping sunrise to π/2 
  #'  and sunset to 3π/2. Sunrise & sunset times are determined based on the dates 
  #'  and locations provided
  sunTime <- sunTime(megadata$radTime, Dates = megadata$DateTime, Coords = xy)
  head(sunTime)
  
  #'  Add to larger data set
  megadata$sunTime <- sunTime
  
  ####  Extract independent detections for wildlife species  ####
  #'  Create a column identifying whether each image is an "independent" event
  #'  If camera site is diff from previous row then give unique value. If not then...
  #'  If species detected is diff from previous row at same site then give unique value. If not then...
  #'  If DateTime is >30 min from previous DateTime at same site for same species then give unique value. If not then...
  #'  Capture value is the same as that in the previous row.
  dat <- arrange(megadata, CameraLocation, DateTime)
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
  
  #'  Retain only the first image from each unique detection event
  detections <- capdata %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  # detections_OK <- detections %>%
  #   filter(grepl("OK", CameraLocation))
  
  #' #'  Retain only the last image from each unique detection event
  #' last_det <- capdata %>%
  #'   group_by(caps) %>%
  #'   slice_tail(n = 1) %>%
  #'   ungroup()
  #' last_det_OK <- last_det %>%
  #'   filter(grepl("OK", CameraLocation))
  
  #'  Filter data to specific date ranges
  seasonal_filter <- function(dets, yr) {
    summer <- dets %>%
      filter(Date > paste0(yr, "-05-31")) %>%
      filter(Date < paste0(yr, "-10-01")) %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species")
    fall <- dets %>%
      filter(Date > paste0(yr, "-09-30")) %>%
      filter(Date < paste0(yr, "-12-01")) %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species")
    winter <- dets %>%
      filter(Date > paste0(yr, "-11-30")) %>%
      filter(Date < paste0(yr+1, "-04-01")) %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species")
    spring <- dets %>%
      filter(Date > paste0(yr+1, "-03-29")) %>%
      filter(Date < paste0(yr+1, "-06-01")) %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species")
    det_list <- list(summer, fall, winter, spring)
    return(det_list)
  }
  #'  Create list of seasonal detection data per year
  dets18 <- seasonal_filter(detections, yr = 2018)
  dets19 <- seasonal_filter(detections, yr = 2019)
  dets20 <- seasonal_filter(detections, yr = 2020)
  #'  Join all detection data per season, across years
  dets_smr <- rbind(dets18[[1]], dets19[[1]], dets20[[1]])
  dets_fall <- rbind(dets18[[2]], dets19[[2]], dets20[[2]])
  dets_wtr <- rbind(dets18[[3]], dets19[[3]], dets20[[3]])
  dets_sprg <- rbind(dets18[[4]], dets19[[4]], dets20[[4]])

  #'  ------------------------------------------
  ####  Predator-prey temporal overlap analysis ####
  #'  ------------------------------------------
  #'  Function to estimate temporal overlap between predators (spp1) and prey (spp2)
  #'  at camera sites located in areas with high or low background predation risk.
  pred_prey_overlap <- function(spp1, spp2, spp3, name1, name2, name3, nboot, dhat) { #i
    #'  Create logical vectors (T/F) indicating whether spp1 & spp2 were detected 
    #'  at the same site and reduce detection events to just those cameras --> These 
    #'  species need to spatially overlap for any temporal overlap to be meaningful
    both.present <- spp1$CameraLocation %in% spp2$CameraLocation
    spp1_dat <- cbind(spp1, both.present) 
    spp1_dat <- spp1_dat[spp1_dat$both.present == T,]
    both.present <- spp2$CameraLocation %in% spp1$CameraLocation
    spp2_dat <- cbind(spp2, both.present)
    spp2_dat <- spp2_dat[spp2_dat$both.present == T,]
    #'  Double check the same number of cameras are included in each
    length(unique(spp1_dat$CameraLocation)); length(unique(spp2_dat$CameraLocation))
    
    #'  Create a logical vector (T/F) indicating whether background risk was high 
    #'  or low at each site where spp1 was detected (based on were both were detected)
    risk <- spp1_dat$CameraLocation %in% stations_data$CameraLocation[stations_data$backgroundRisk == "Low"]
    spp1_dat <- cbind(spp1_dat, risk) %>% mutate(risk = ifelse(risk == TRUE, "Low", "High"))
    spp1_lowrisk <- filter(spp1_dat, spp1_dat$risk == "Low")
    spp1_highrisk <- filter(spp1_dat, spp1_dat$risk == "High")
    
    #'  Create a logical vector (T/F) indicating whether background risk was high 
    #'  or low at each site where spp2 was detected (based on were both were detected)
    risk <- spp2_dat$CameraLocation %in% stations_data$CameraLocation[stations_data$backgroundRisk == "Low"]
    spp2_dat <- cbind(spp2_dat, risk) %>% mutate(risk = ifelse(risk == TRUE, "Low", "High"))
    spp2_lowrisk <- filter(spp2_dat, spp2_dat$risk == "Low")
    spp2_highrisk <- filter(spp2_dat, spp2_dat$risk == "High")
    
    #'  Review sample size per species- smaller sample will determine which coefficient
    #'  of overlap estimator to use (Meredith & Ridout 2017 suggest using delta1 
    #'  for small samples [<50 detection events] and delta4 for larger samples 
    #'  [>50 detection events]).
    ndet_spp1_lowrisk <- nrow(spp1_lowrisk)
    ndet_spp2_lowrisk <- nrow(spp2_lowrisk)
    ndet_spp1_highrisk <- nrow(spp1_highrisk)
    ndet_spp2_highrisk <- nrow(spp2_highrisk)
    print(ndet_spp1_lowrisk); print(ndet_spp2_lowrisk)
    print(ndet_spp1_highrisk); print(ndet_spp2_highrisk)
    
    #'  Visualize general temporal activity with density plots
    densityPlot(spp1$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name1, " Daily Activity"))
    densityPlot(spp2$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name2, " Daily Activity"))
    
    #'  Visualize temporal overlap
    #'  Overlap when background risk is low
    saveOverlap_lowrisk <- overlapPlot(spp1_lowrisk$sunTime, spp2_lowrisk$sunTime, rug = T, 
                                       xscale = NA, xcenter = "noon", linet = c(1, 1), linec = c("red", "blue"), 
                                       linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " (red) & ", name2, " (blue) \ndiel activity when background risk is low")) 
    
    #'  Wrangle density data from wide to long format
    DensityA_low <- saveOverlap_lowrisk[,1:2] %>%
      mutate(Species = "Predator",
             BackgroundRisk = "Low")
    DensityB_low <- saveOverlap_lowrisk[,c(1,3)] %>%
      mutate(Species = "Prey",
             BackgroundRisk = "Low")
    overlap_lowrisk <- full_join(DensityA_low, DensityB_low, by = c("x", "BackgroundRisk")) 
    
    #'  Overlap when background risk is high
    saveOverlap_highrisk <- overlapPlot(spp1_highrisk$sunTime, spp2_highrisk$sunTime, rug = T, 
                                      xscale = NA, xcenter = "noon", linet = c(1, 1), linec = c("red", "blue"), 
                                      linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " (red) & ", name2, " (blue) \ndiel activity when background risk is high")) 
    
    #'  Wrangle from wide to long format
    DensityA_high <- saveOverlap_highrisk[,1:2] %>%
      mutate(Species = "Predator",
             BackgroundRisk = "High")
    DensityB_high <- saveOverlap_highrisk[,c(1,3)] %>%
      mutate(Species = "Prey",
             BackgroundRisk = "High")
    overlap_highrisk <- full_join(DensityA_high, DensityB_high, by = c("x", "BackgroundRisk")) 
    
    #'  Bind into single long data set of density estimates to make custom overlap plots
    plotdata <- rbind(overlap_lowrisk, overlap_highrisk)
    
    #'  Calculate coefficient of overlap
    dhats_spp1.spp2_low <- overlapEst(A = spp1_lowrisk$sunTime,
                                      B = spp2_lowrisk$sunTime, type = dhat)
    dhats_spp1.spp2_high <- overlapEst(A = spp1_highrisk$sunTime,
                                       B = spp2_highrisk$sunTime, type = dhat)
    
    #'  Bootstrap to estimate standard errors
    #'  FYI: smooth = TRUE is default and allows bootstrap to randomly sample from
    #'  a distribution of times that have a wider range than the original sample
    #'  (see pg. 5 in Overlap package vignette for details)
    spp12.lowrisk.boot <- bootstrap(spp1_lowrisk$sunTime, spp2_lowrisk$sunTime,
                                    nboot, smooth = TRUE, type = dhat)
    spp12.highrisk.boot <- bootstrap(spp1_highrisk$sunTime, spp2_highrisk$sunTime,
                                     nboot, smooth = TRUE, type = dhat)
    #'  Bootstrap mean will be a little different then detla coefficient due to
    #'  bootstrap bias (BSmean - delta) that needs to be accounted for in 95% CIs
    BSmean.lowrisk <- mean(spp12.lowrisk.boot)
    BSmean.highrisk <- mean(spp12.highrisk.boot)
    
    #'  Bootstrap 95% Confidence Intervals
    #'  norm0 uses the standard deviation of bootstrap results to calculate CI (delta +/- 1.96*SDboot)
    #'  basic0 takes the 2.5% and 97.5% percentiles and adjusts based on BS bias (percentile - BSbias)
    #'  If sampling distribution is normal, norm0 and basic0 should be similar;
    #'  if sampling distribution is skewed (i.e., if delta is close to 0 or 1) then
    #'  basic0 is the better estimator - using basic0 b/c should be good in either situation
    #'  Using bootCIlogit instead of bootCI so that bias corrections are done on
    #'  the logit scale, then backtransformed. Without this, 95% CIs can fall
    #'  outside (0, 1) interval. See Overlap vignette for more details.
    lowrisk_CI <- bootCIlogit(dhats_spp1.spp2_low, spp12.lowrisk.boot) #[i]
    highrisk_CI <- bootCIlogit(dhats_spp1.spp2_high, spp12.highrisk.boot) #[i]
    
    #'  Print results
    #'  Effect of low background risk
    print("Overlap coefficients when background risk is low"); print(dhats_spp1.spp2_low)
    print("Bootstrap mean"); print(BSmean.lowrisk)
    print("Bootstrap 95% CI"); print(lowrisk_CI)
    
    #'  Effect of high background risk
    print("Overlap coefficients when background risk is high"); print(dhats_spp1.spp2_high)
    print("Bootstrap mean"); print(BSmean.highrisk)
    print("Bootstrap 95% CI"); print(highrisk_CI)
    
    #'  Save as a giant list
    overlap_list <- list(dhats_spp1.spp2_low, dhats_spp1.spp2_high,
                         spp12.lowrisk.boot, spp12.highrisk.boot,
                         lowrisk_CI, highrisk_CI, ndet_spp1_lowrisk,
                         ndet_spp2_lowrisk, ndet_spp1_highrisk, ndet_spp2_highrisk,
                         plotdata)
    names(overlap_list) <- c("dhat_lowrisk", "dhat_highrisk", "dhat_lowrisk_boot",
                             "dhat_highrisk_boot", "lowrisk_CI", "highrisk_CI",
                             "ndet_spp1_lowrisk", "ndet_spp2_lowrisk",
                             "ndet_spp1_highrisk", "ndet_spp2_highrisk",
                             "overlap.plot.data")
    
    return(overlap_list)
  }
  ####  Predator-Prey Overlap Grazing Season  ####
  #'  Estimate temporal overlap between predators and prey when cattle are/aren't detected
  #'  Focusing on only OK study area since big difference in number of cameras 
  #'  with cattle in NE vs OK, pooling across study areas can bias results
  nboot <- 100
  ####  Cougar - Mule deer  ####
  coug_md_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Cougar"), 
                                          spp2 = filter(dets_smr, Species == "Mule Deer"), 
                                          name1 = "Cougar", name2 = "Mule Deer", 
                                          nboot = nboot, dhat = "Dhat4")
  coug_md_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Cougar"), 
                                        spp2 = filter(dets_fall, Species == "Mule Deer"), 
                                        name1 = "Cougar", name2 = "Mule Deer", 
                                        nboot = nboot, dhat = "Dhat1")
  coug_md_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Cougar"), 
                                        spp2 = filter(dets_wtr, Species == "Mule Deer"), 
                                        name1 = "Cougar", name2 = "Mule Deer", 
                                        nboot = nboot, dhat = "Dhat1")
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_md_sprg_over1 <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Cougar"), 
                                        spp2 = filter(dets_sprg, Species == "Mule Deer"), 
                                        name1 = "Cougar", name2 = "Mule Deer", 
                                        nboot = nboot, dhat = "Dhat1")
  coug_md_sprg_over4 <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Cougar"), 
                                         spp2 = filter(dets_sprg, Species == "Mule Deer"), 
                                         name1 = "Cougar", name2 = "Mule Deer", 
                                         nboot = nboot, dhat = "Dhat4")
  ####  Cougar - Elk  ####
  coug_elk_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Cougar"), 
                                        spp2 = filter(dets_smr, Species == "Elk"), 
                                        name1 = "Cougar", name2 = "Elk", 
                                        nboot = nboot, dhat = "Dhat4")
  coug_elk_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Cougar"), 
                                         spp2 = filter(dets_fall, Species == "Elk"), 
                                         name1 = "Cougar", name2 = "Elk", 
                                         nboot = nboot, dhat = "Dhat1")
  coug_elk_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Cougar"), 
                                        spp2 = filter(dets_wtr, Species == "Elk"), 
                                        name1 = "Cougar", name2 = "Elk", 
                                        nboot = nboot, dhat = "Dhat1")
  coug_elk_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Cougar"), 
                                          spp2 = filter(dets_sprg, Species == "Elk"), 
                                          name1 = "Cougar", name2 = "Elk", 
                                          nboot = nboot, dhat = "Dhat1")
  ####  Cougar - Moose  ####
  coug_moose_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Cougar"), 
                                         spp2 = filter(dets_smr, Species == "Moose"), 
                                         name1 = "Cougar", name2 = "Moose", 
                                         nboot = nboot, dhat = "Dhat4")
  coug_moose_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Cougar"), 
                                          spp2 = filter(dets_fall, Species == "Moose"), 
                                          name1 = "Cougar", name2 = "Moose", 
                                          nboot = nboot, dhat = "Dhat1")
  coug_moose_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Cougar"), 
                                         spp2 = filter(dets_wtr, Species == "Moose"), 
                                         name1 = "Cougar", name2 = "Moose", 
                                         nboot = nboot, dhat = "Dhat1")
  coug_moose_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Cougar"), 
                                           spp2 = filter(dets_sprg, Species == "Moose"), 
                                           name1 = "Cougar", name2 = "Moose", 
                                           nboot = nboot, dhat = "Dhat1")
  ####  Cougar - White-tailed deer  ####
  coug_wtd_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Cougar"), 
                                           spp2 = filter(dets_smr, Species == "White-tailed Deer"), 
                                           name1 = "Cougar", name2 = "wtd", 
                                           nboot = nboot, dhat = "Dhat4")
  #'  Low risk <50 cougars; High risk >50 cougars
  coug_wtd_fall_over1 <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Cougar"), 
                                            spp2 = filter(dets_fall, Species == "White-tailed Deer"), 
                                            name1 = "Cougar", name2 = "wtd", 
                                            nboot = nboot, dhat = "Dhat1")
  coug_wtd_fall_over4 <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Cougar"), 
                                          spp2 = filter(dets_fall, Species == "White-tailed Deer"), 
                                          name1 = "Cougar", name2 = "wtd", 
                                          nboot = nboot, dhat = "Dhat4")
  #'  Low risk >50 cougars; High risk <50 cougars
  coug_wtd_wtr_over1 <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Cougar"), 
                                           spp2 = filter(dets_wtr, Species == "White-tailed Deer"), 
                                           name1 = "Cougar", name2 = "wtd", 
                                           nboot = nboot, dhat = "Dhat1")
  coug_wtd_wtr_over4 <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Cougar"), 
                                         spp2 = filter(dets_wtr, Species == "White-tailed Deer"), 
                                         name1 = "Cougar", name2 = "wtd", 
                                         nboot = nboot, dhat = "Dhat4")
  coug_wtd_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Cougar"), 
                                             spp2 = filter(dets_sprg, Species == "White-tailed Deer"), 
                                             name1 = "Cougar", name2 = "wtd", 
                                             nboot = nboot, dhat = "Dhat1")
  ####  Wolf - Mule deer  ####
  wolf_md_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Wolf"), 
                                        spp2 = filter(dets_smr, Species == "Mule Deer"), 
                                        name1 = "Wolf", name2 = "Mule Deer", 
                                        nboot = nboot, dhat = "Dhat1")
  wolf_md_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Wolf"), 
                                         spp2 = filter(dets_fall, Species == "Mule Deer"), 
                                         name1 = "Wolf", name2 = "Mule Deer", 
                                         nboot = nboot, dhat = "Dhat1")
  # wolf_md_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Wolf"), 
  #                                       spp2 = filter(dets_wtr, Species == "Mule Deer"), 
  #                                       name1 = "Wolf", name2 = "Mule Deer", 
  #                                       nboot = nboot, dhat = "Dhat1")
  wolf_md_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Wolf"), 
                                          spp2 = filter(dets_sprg, Species == "Mule Deer"), 
                                          name1 = "Wolf", name2 = "Mule Deer", 
                                          nboot = nboot, dhat = "Dhat1")
  ####  Wolf - Elk  ####
  wolf_elk_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Wolf"), 
                                        spp2 = filter(dets_smr, Species == "Elk"), 
                                        name1 = "Wolf", name2 = "Elk", 
                                        nboot = nboot, dhat = "Dhat1")
  # wolf_elk_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Wolf"), 
  #                                        spp2 = filter(dets_fall, Species == "Elk"), 
  #                                        name1 = "Wolf", name2 = "Elk", 
  #                                        nboot = nboot, dhat = "Dhat1")
  # wolf_elk_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Wolf"),
  #                                       spp2 = filter(dets_wtr, Species == "Elk"),
  #                                       name1 = "Wolf", name2 = "Elk",
  #                                       nboot = nboot, dhat = "Dhat1")
  # wolf_elk_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Wolf"), 
  #                                        spp2 = filter(dets_sprg, Species == "Elk"), 
  #                                        name1 = "Wolf", name2 = "Elk", 
  #                                        nboot = nboot, dhat = "Dhat1")
  ####  Wolf - Moose  ####
  wolf_moose_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Wolf"), 
                                        spp2 = filter(dets_smr, Species == "Moose"), 
                                        name1 = "Wolf", name2 = "Moose", 
                                        nboot = nboot, dhat = "Dhat1")
  wolf_moose_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Wolf"), 
                                         spp2 = filter(dets_fall, Species == "Moose"), 
                                         name1 = "Wolf", name2 = "Moose", 
                                         nboot = nboot, dhat = "Dhat1")
  wolf_moose_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Wolf"),
                                        spp2 = filter(dets_wtr, Species == "Moose"),
                                        name1 = "Wolf", name2 = "Moose",
                                        nboot = nboot, dhat = "Dhat1")
  # wolf_moose_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Wolf"), 
  #                                        spp2 = filter(dets_sprg, Species == "Moose"), 
  #                                        name1 = "Wolf", name2 = "Moose", 
  #                                        nboot = nboot, dhat = "Dhat1")
  ####  Wolf - White-tailed deer ####
  wolf_wtd_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Wolf"), 
                                           spp2 = filter(dets_smr, Species == "White-tailed Deer"), 
                                           name1 = "Wolf", name2 = "White-tailed Deer", 
                                           nboot = nboot, dhat = "Dhat1")
  wolf_wtd_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Wolf"), 
                                            spp2 = filter(dets_fall, Species == "White-tailed Deer"), 
                                            name1 = "Wolf", name2 = "White-tailed Deer", 
                                            nboot = nboot, dhat = "Dhat1")
  wolf_wtd_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Wolf"),
                                           spp2 = filter(dets_wtr, Species == "White-tailed Deer"),
                                           name1 = "Wolf", name2 = "White-tailed Deer",
                                           nboot = nboot, dhat = "Dhat1")
  # wolf_wtd_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Wolf"),
  #                                        spp2 = filter(dets_sprg, Species == "White-tailed Deer"),
  #                                        name1 = "Wolf", name2 = "White-tailed Deer",
  #                                        nboot = nboot, dhat = "Dhat1")
  ####  Black bear - Mule Deer  ####
  bear_md_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Black Bear"), 
                                        spp2 = filter(dets_smr, Species == "Mule Deer"), 
                                        name1 = "Black Bear", name2 = "Mule Deer", 
                                        nboot = nboot, dhat = "Dhat4")
  #'  Low risk >50 black bears; High risk <50 black bears
  bear_md_fall_over1 <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Black Bear"), 
                                         spp2 = filter(dets_fall, Species == "Mule Deer"), 
                                         name1 = "Black Bear", name2 = "Mule Deer", 
                                         nboot = nboot, dhat = "Dhat1")
  bear_md_fall_over4 <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Black Bear"), 
                                         spp2 = filter(dets_fall, Species == "Mule Deer"), 
                                         name1 = "Black Bear", name2 = "Mule Deer", 
                                         nboot = nboot, dhat = "Dhat4")
  # bear_md_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Black Bear"),
  #                                       spp2 = filter(dets_wtr, Species == "Mule Deer"),
  #                                       name1 = "Black Bear", name2 = "Mule Deer",
  #                                       nboot = nboot, dhat = "Dhat1")
  bear_md_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Black Bear"), 
                                         spp2 = filter(dets_sprg, Species == "Mule Deer"), 
                                         name1 = "Black Bear", name2 = "Mule Deer", 
                                         nboot = nboot, dhat = "Dhat4")
  ####  Black Bear - Elk  ####
  bear_elk_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Black Bear"), 
                                        spp2 = filter(dets_smr, Species == "Elk"), 
                                        name1 = "Black Bear", name2 = "Elk", 
                                        nboot = nboot, dhat = "Dhat4")
  bear_elk_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Black Bear"), 
                                          spp2 = filter(dets_fall, Species == "Elk"), 
                                          name1 = "Black Bear", name2 = "Elk", 
                                          nboot = nboot, dhat = "Dhat1")
  # bear_elk_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Black Bear"),
  #                                       spp2 = filter(dets_wtr, Species == "Elk"),
  #                                       name1 = "Black Bear", name2 = "Elk",
  #                                       nboot = nboot, dhat = "Dhat1")
  bear_elk_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Black Bear"), 
                                         spp2 = filter(dets_sprg, Species == "Elk"), 
                                         name1 = "Black Bear", name2 = "Elk", 
                                         nboot = nboot, dhat = "Dhat1")
  ####  Black bear - Moose  ####
  bear_moose_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Black Bear"), 
                                         spp2 = filter(dets_smr, Species == "Moose"), 
                                         name1 = "Black Bear", name2 = "Moose", 
                                         nboot = nboot, dhat = "Dhat4")
  bear_moose_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Black Bear"), 
                                          spp2 = filter(dets_fall, Species == "Moose"), 
                                          name1 = "Black Bear", name2 = "Moose", 
                                          nboot = nboot, dhat = "Dhat1")
  # bear_moose_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Black Bear"),
  #                                       spp2 = filter(dets_wtr, Species == "Moose"),
  #                                       name1 = "Black Bear", name2 = "Moose",
  #                                       nboot = nboot, dhat = "Dhat1")
  bear_moose_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Black Bear"), 
                                          spp2 = filter(dets_sprg, Species == "Moose"), 
                                          name1 = "Black Bear", name2 = "Moose", 
                                          nboot = nboot, dhat = "Dhat1")
  ####  Black bear - White-tailed Deer  ####
  bear_wtd_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Black Bear"), 
                                           spp2 = filter(dets_smr, Species == "White-tailed Deer"), 
                                           name1 = "Black Bear", name2 = "wtd", 
                                           nboot = nboot, dhat = "Dhat4")
  bear_wtd_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Black Bear"), 
                                            spp2 = filter(dets_fall, Species == "White-tailed Deer"), 
                                            name1 = "Black Bear", name2 = "wtd", 
                                            nboot = nboot, dhat = "Dhat1")
  # bear_wtd_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Black Bear"),
  #                                       spp2 = filter(dets_wtr, Species == "White-tailed Deer"),
  #                                       name1 = "Black Bear", name2 = "wtd",
  #                                       nboot = nboot, dhat = "Dhat1")
  #'  Low risk >50 black bears; High risk <50 black bears
  bear_wtd_sprg_over1 <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Black Bear"), 
                                            spp2 = filter(dets_sprg, Species == "White-tailed Deer"), 
                                            name1 = "Black Bear", name2 = "wtd", 
                                            nboot = nboot, dhat = "Dhat1")
  bear_wtd_sprg_over4 <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Black Bear"), 
                                          spp2 = filter(dets_sprg, Species == "White-tailed Deer"), 
                                          name1 = "Black Bear", name2 = "wtd", 
                                          nboot = nboot, dhat = "Dhat4")
  ####  Bobcat - Mule deer  ####
  bob_md_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Bobcat"), 
                                        spp2 = filter(dets_smr, Species == "Mule Deer"), 
                                        name1 = "Bobcat", name2 = "Mule Deer", 
                                        nboot = nboot, dhat = "Dhat4")
  bob_md_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Bobcat"), 
                                          spp2 = filter(dets_fall, Species == "Mule Deer"), 
                                          name1 = "Bobcat", name2 = "Mule Deer", 
                                          nboot = nboot, dhat = "Dhat1")
  # bob_md_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Bobcat"),
  #                                       spp2 = filter(dets_wtr, Species == "Mule Deer"),
  #                                       name1 = "Bobcat", name2 = "Mule Deer",
  #                                       nboot = nboot, dhat = "Dhat1")
  bob_md_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Bobcat"), 
                                         spp2 = filter(dets_sprg, Species == "Mule Deer"), 
                                         name1 = "Bobcat", name2 = "Mule Deer", 
                                         nboot = nboot, dhat = "Dhat1")
  ####  Bobcat - White-tailed Deer  ####
  bob_wtd_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Bobcat"), 
                                       spp2 = filter(dets_smr, Species == "White-tailed Deer"), 
                                       name1 = "Bobcat", name2 = "White-tailed Deer", 
                                       nboot = nboot, dhat = "Dhat4")
  bob_wtd_fall_over1 <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Bobcat"), 
                                        spp2 = filter(dets_fall, Species == "White-tailed Deer"), 
                                        name1 = "Bobcat", name2 = "White-tailed Deer", 
                                        nboot = nboot, dhat = "Dhat1")
  bob_wtd_fall_over4 <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Bobcat"), 
                                         spp2 = filter(dets_fall, Species == "White-tailed Deer"), 
                                         name1 = "Bobcat", name2 = "White-tailed Deer", 
                                         nboot = nboot, dhat = "Dhat4")
  bob_wtd_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Bobcat"),
                                        spp2 = filter(dets_wtr, Species == "White-tailed Deer"),
                                        name1 = "Bobcat", name2 = "White-tailed Deer",
                                        nboot = nboot, dhat = "Dhat1")
  bob_wtd_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Bobcat"), 
                                        spp2 = filter(dets_sprg, Species == "White-tailed Deer"), 
                                        name1 = "Bobcat", name2 = "White-tailed Deer", 
                                        nboot = nboot, dhat = "Dhat1")
  ####  Coyote - Mule deer  ####
  coy_md_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Coyote"), 
                                       spp2 = filter(dets_smr, Species == "Mule Deer"), 
                                       name1 = "Coyote", name2 = "Mule Deer", 
                                       nboot = nboot, dhat = "Dhat4")
  coy_md_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Coyote"), 
                                        spp2 = filter(dets_fall, Species == "Mule Deer"), 
                                        name1 = "Coyote", name2 = "Mule Deer", 
                                        nboot = nboot, dhat = "Dhat4")
  coy_md_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Coyote"),
                                        spp2 = filter(dets_wtr, Species == "Mule Deer"),
                                        name1 = "Coyote", name2 = "Mule Deer",
                                        nboot = nboot, dhat = "Dhat4")
  coy_md_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Coyote"), 
                                        spp2 = filter(dets_sprg, Species == "Mule Deer"), 
                                        name1 = "Coyote", name2 = "Mule Deer", 
                                        nboot = nboot, dhat = "Dhat4")
  ####  Coyote - White-tailed Deer  ####
  coy_wtd_smr_over <- pred_prey_overlap(spp1 = filter(dets_smr, Species == "Coyote"), 
                                        spp2 = filter(dets_smr, Species == "White-tailed Deer"), 
                                        name1 = "Coyote", name2 = "White-tailed Deer", 
                                        nboot = nboot, dhat = "Dhat4")
  coy_wtd_fall_over <- pred_prey_overlap(spp1 = filter(dets_fall, Species == "Coyote"), 
                                          spp2 = filter(dets_fall, Species == "White-tailed Deer"), 
                                          name1 = "Coyote", name2 = "White-tailed Deer", 
                                          nboot = nboot, dhat = "Dhat4")
  coy_wtd_wtr_over <- pred_prey_overlap(spp1 = filter(dets_wtr, Species == "Coyote"),
                                        spp2 = filter(dets_wtr, Species == "White-tailed Deer"),
                                        name1 = "Coyote", name2 = "White-tailed Deer",
                                        nboot = nboot, dhat = "Dhat4")
  coy_wtd_sprg_over <- pred_prey_overlap(spp1 = filter(dets_sprg, Species == "Coyote"), 
                                         spp2 = filter(dets_sprg, Species == "White-tailed Deer"), 
                                         name1 = "Coyote", name2 = "White-tailed Deer", 
                                         nboot = nboot, dhat = "Dhat4")
  
  #'  Save this monster
  coug_prey_overlap <- list(coug_md_smr_over, coug_md_fall_over, coug_md_wtr_over, coug_md_sprg_over1, coug_md_sprg_over4,
                            coug_elk_smr_over, coug_elk_fall_over, coug_elk_wtr_over, coug_elk_sprg_over,
                            coug_moose_smr_over, coug_moose_fall_over, coug_moose_wtr_over, coug_moose_sprg_over,
                            coug_wtd_smr_over, coug_wtd_fall_over1, coug_wtd_fall_over4, coug_wtd_wtr_over1, coug_wtd_wtr_over4, coug_wtd_sprg_over)
  wolf_prey_overlap <- list(wolf_md_smr_over, wolf_md_fall_over, wolf_md_sprg_over,
                            wolf_elk_smr_over, wolf_moose_smr_over, wolf_moose_fall_over, wolf_moose_wtr_over,
                            wolf_wtd_smr_over, wolf_wtd_fall_over, wolf_wtd_wtr_over)
  bear_prey_overlap <- list(bear_md_smr_over, bear_md_fall_over1, bear_md_fall_over4, bear_md_sprg_over,
                            bear_elk_smr_over, bear_elk_fall_over, bear_elk_sprg_over,
                            bear_moose_smr_over, bear_moose_fall_over, bear_moose_sprg_over,
                            bear_wtd_smr_over, bear_wtd_fall_over, bear_wtd_sprg_over1, bear_wtd_sprg_over4)
  bob_prey_overlap <- list(bob_md_smr_over, bob_md_fall_over, bob_md_sprg_over,
                           bob_wtd_smr_over, bob_wtd_fall_over1, bob_wtd_fall_over4, bob_wtd_wtr_over, bear_wtd_sprg_over)
  coy_prey_overlap <- list(coy_md_smr_over, coy_md_fall_over, coy_md_wtr_over, coy_md_sprg_over,
                           coy_wtd_smr_over, coy_wtd_fall_over, coy_wtd_wtr_over, coy_wtd_sprg_over)
  pred_prey_overlap <- list(coug_prey_overlap, wolf_prey_overlap, bear_prey_overlap, bob_prey_overlap, coy_prey_overlap)
  
  save(pred_prey_overlap, file = paste("./Outputs/Temporal Overlap/PredPrey_LowHi_Overlap_", Sys.Date(), ".RData"))
  
  #'  --------------------------------------------
  ####  Single species temporal overlap analysis  ####
  #'  --------------------------------------------
  #'  Function to estimate differences in temporal activity for a species at camera  
  #'  sites where cattle/humans are present (detected) vs absent (not detected).
  #'  Using spp3 to also represent cameras on public vs private land in this function.
  spp_overlap <- function(spp1, spp2, name1, name2, nboot, dhat) { #, i
    #'  Create logical vectors (T/F) indicating which cameras spp1 was detected 
    #'  at with and without spp2
    both.present <- spp1$CameraLocation %in% spp2$CameraLocation
    spp1_dat <- cbind(spp1, both.present) 
    #'  Split out data into camera locations where both are present vs spp2 absent
    spp1_spp2.present <- spp1_dat[spp1_dat$both.present == T,]
    spp1_spp2.absent <- spp1_dat[spp1_dat$both.present == F,]
    
    #'  Review sample size per species- smaller sample will determine which coefficient
    #'  of overlap estimator to use (delta1 for small samples [<50 detection events], 
    #'  delta4 for larger samples [>75 detection events])
    ndet_spp2.present <- nrow(spp1_spp2.present)
    ndet_spp2.absent <- nrow(spp1_spp2.absent)
    print(ndet_spp2.present); print(ndet_spp2.absent)
    
    #' #'  Visualize general temporal activity with density plots
    densityPlot(spp1$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name1, " Daily Activity"))
    densityPlot(spp2$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name2, " Daily Activity"))
    
    #'  Visualize temporal overlap
    saveOverlap <- overlapPlot(spp1_spp2.present$sunTime, spp1_spp2.absent$sunTime, rug = T, 
                               xscale = NA, xcenter = "noon", linet = c(1, 1), linec = c("red", "blue"), 
                               linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " diel activity \nwhen ", name2, " are Present and Absente")) 
    saveDensity <- densityPlot(spp2$sunTime, add = T, xscale = NA, linec = "black", lwd = 2, lty = 2, extend = NULL)
    legend("topleft", c("Anthro activity present", "Anthro activity  absent", "Anthro activity"), lty=c(1, 1, 2), col=c("red", "blue", "black"), bg = "white", bty = "n")
    
    #'  Wrangle density data from wide to long format
    DensityA <- saveOverlap[,1:2] %>%
      mutate(Anthro_Activity = "Present")
    DensityB <- saveOverlap[,c(1,3)] %>%
      mutate(Anthro_Activity = "Absent")
    plotdata <- full_join(DensityA, DensityB, by = "x") %>%
      full_join(saveDensity, by ="x") %>%
      mutate(Species.z = name2)
    
    # plotdata <- full_join(saveOverlap, saveDensity, by = "x")
    # colnames(plotdata) <- c("x", "DensityA", "DensityB", "Anthro_Activity")
    
    #'  Calculate coefficient of overlap
    dhats_spp1.spp2 <- overlapEst(A = spp1_spp2.present$sunTime, 
                                  B = spp1_spp2.absent$sunTime, type = dhat) 
    
    #'  Bootstrap to estimate standard errors
    #'  FYI: smooth = TRUE is default and allows bootstrap to randomly sample from 
    #'  a distribution of times that have a wider range than the original sample
    #'  (see pg. 5 in Overlap package vignette for details) 
    spp1.spp2.boot <- bootstrap(spp1_spp2.present$sunTime, spp1_spp2.absent$sunTime, 
                                nboot, smooth = TRUE, type = dhat)  
    #'  Bootstrap mean will be a little different then detla coefficient due to
    #'  bootstrap bias (BSmean - delta) that needs to be accounted for in 95% CIs
    BSmean <- mean(spp1.spp2.boot)
    
    #'  Bootstrap 95% Confidence Intervals
    #'  norm0 uses the standard deviation of bootstrap results to calculate CI (delta +/- 1.96*SDboot)
    #'  basic0 takes the 2.5% and 97.5% percentiles and adjusts based on BS bias (percentile - BSbias)
    #'  If sampling distribution is normal, norm0 and basic0 should be similar;
    #'  if sampling distribution is skewed (i.e., if delta is close to 0 or 1) then
    #'  basic0 is the better estimator
    #'  Using bootCIlogit instead of bootCI so that bias corrections are done on
    #'  the logit scale, then backtransformed. Without this, 95% CIs can fall
    #'  outside (0, 1) interval. See Overlap vignette for more details.
    CI <- bootCIlogit(dhats_spp1.spp2, spp1.spp2.boot) #dhats_spp1.spp2[i]
    
    #'  Print results
    #'  Effect of spp2 being present
    print("Overlap coefficients when spp2 is present"); print(dhats_spp1.spp2)
    print("Bootstrap mean"); print(BSmean)
    print("Bootstrap 95% CI"); print(CI)
    
    #'  Save as a giant list
    overlap_list <- list(dhats_spp1.spp2, spp1.spp2.boot, CI, ndet_spp2.present, ndet_spp2.absent, plotdata)
    names(overlap_list) <- c("dhats_spp1.spp2", "spp1.spp2.boot", "CI", "ndet_spp2.present", "ndet_spp2.absent", "overlap.plot.data")
    
    return(overlap_list)
  }
  #'  Estimate temporal overlap for a species when cattle are/aren't detected
  #'  Focusing on only OK study area since big difference in number of cameras 
  #'  with cattle in NE vs OK, pooling across study areas could be confounding
  #'  any apparent temporal patterns
  ####  Single-Species Overlap Grazing Season  ####
  coug_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"),
                                 spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                 name1 = "Cougar", name2 = "Cattle", nboot = 10000, dhat = "Dhat1") #i = 1
  
  