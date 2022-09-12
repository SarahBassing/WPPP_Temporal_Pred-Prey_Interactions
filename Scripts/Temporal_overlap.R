  #'  ==================================================================
  #'  Temporal overlap patterns in response to background predation risk
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  June 2022
  #'  ==================================================================
  #'  Script to estimate temporal overlap of prey species in response to varying
  #'  levels of background predation risk in a multi-predator system. 
  #'  ==================================================================
  
  #'  Load libraries
  library(lubridate)
  library(chron)
  library(overlap)
  library(circular)
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
  
  #'  Identity which predators were detected at each camera site per season (e.g., predation risk)
  #'  Split up megadetection data by season acorss years
  smr <- megadata %>%
    filter(Date > "2018-05-31" & Date < "2018-10-01" | Date > "2019-05-31" & Date < "2019-10-01" | Date > "2020-05-31" & Date < "2020-10-01") 
  fll <- megadata %>%
    filter(Date > "2018-09-30" & Date < "2018-12-01" | Date > "2019-09-30" & Date < "2019-12-01" | Date > "2020-09-30" & Date < "2020-12-01") 
  wtr <- megadata %>%
    filter(Date > "2018-11-30" & Date < "2019-04-01" | Date > "2019-11-30" & Date < "2020-04-01" | Date > "2020-11-30" & Date < "2021-04-01") 
  spg <- megadata %>%
    filter(Date > "2019-03-29" & Date < "2019-06-01" | Date > "2020-03-29" & Date < "2020-06-01" | Date > "2021-03-29" & Date < "2021-06-01") 
  
  #'  Function to identify which cameras had a predator detection of interest
  seasonal_dets <- function(seasondata, spp){
    dets <- stations_data$CameraLocation %in% seasondata$CameraLocation[seasondata$Species == spp]
    binary_cov <- as.data.frame(dets) %>% transmute(binary_cov = ifelse(dets == TRUE, 1, 0))
    return(binary_cov)
  }
  #'  Create list of seasonal datasets
  season_list <- list(smr, fll, wtr, spg)
  #'  Generate binary variable indicating whether predator was detected at a camera in each season
  bear_dets <- lapply(season_list, seasonal_dets, spp = "Black Bear")
  bear_det <- as.data.frame(bear_dets)
  colnames(bear_det) <- c("Summer_bear", "Fall_bear", "Winter_bear", "Spring_bear")
  bob_dets <- lapply(season_list, seasonal_dets, spp = "Bobcat")
  bob_det <- as.data.frame(bob_dets)
  colnames(bob_det) <- c("Summer_bobcat", "Fall_bobcat", "Winter_bobcat", "Spring_bobcat")
  coug_dets <- lapply(season_list, seasonal_dets, spp = "Cougar")
  coug_det <- as.data.frame(coug_dets)
  colnames(coug_det) <- c("Summer_cougar", "Fall_cougar", "Winter_cougar", "Spring_cougar")
  coy_dets <- lapply(season_list, seasonal_dets, spp = "Coyote")
  coy_det <- as.data.frame(coy_dets)
  colnames(coy_det) <- c("Summer_coyote", "Fall_coyote", "Winter_coyote", "Spring_coyote")
  lynx_dets <- lapply(season_list, seasonal_dets, spp = "Lynx")
  lynx_det <- as.data.frame(lynx_dets)
  colnames(lynx_det) <- c("Summer_lynx", "Fall_lynx", "Winter_lynx", "Spring_lynx")
  wolf_dets <- lapply(season_list, seasonal_dets, spp = "Wolf")
  wolf_det <- as.data.frame(wolf_dets)
  colnames(wolf_det) <- c("Summer_wolf", "Fall_wolf", "Winter_wolf", "Spring_wolf")
  
  #'  Add covariates representing background predation risk to
  stations_data <- stations_data %>%
    #'  Create habitat complexity index by multiplying TRI by forest habitat covs
    #'  index 1) 250m spatial scale; index 2) 30m/camera site spatial scale
    mutate(Complexity_index1 = TRI_250m * (PercForest*100), # multiply % forest (in decimals) by 100 so on same scale as canopy cover
           Complexity_index2 = TRI * Canopy_Cov)
  #'  Add naive occupancy for each predator per season to statations dataframe
  stations_data <- cbind(stations_data, bear_det, bob_det, coug_det, coy_det, lynx_det, wolf_det) 
  
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
  
  #'  Add categorical variable designating level of background risk based on
  #'  habitat complexity index (values < mean HCI is low risk, values => mean HCI 
  #'  is high risk; mean HCI = )
  stations_data <- mutate(stations_data, backgroundRisk = ifelse(Complexity_index1 < mean(stations_data$Complexity_index1), "Low", "High"))
  #' #'  Save for time-btwn-detection analyses
  # write.csv(stations_data, "./Data/cam_stations_hab_complex_data.csv")
  
  #' #'  Is naive predator occupancy associated with habitat complexity?
  #' bear_index_cor <- glm(bear_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(bear_index_cor)
  #' bob_index_cor <- glm(bob_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(bob_index_cor)
  #' coug_index_cor <- glm(coug_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(coug_index_cor)
  #' coy_index_cor <- glm(coy_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(coy_index_cor)
  #' lynx_index_cor <- glm(lynx_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(lynx_index_cor)
  #' wolf_index_cor <- glm(wolf_det ~ Complexity_index1, family = binomial(link = "logit"), data = stations_data); summary(wolf_index_cor)
  
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
      mutate(PredPrey = "Predator",
             Species = name1,
             BackgroundRisk = "Low")
    DensityB_low <- saveOverlap_lowrisk[,c(1,3)] %>%
      mutate(PredPrey = "Prey",
             Species = name2,
             BackgroundRisk = "Low")
    overlap_lowrisk <- full_join(DensityA_low, DensityB_low, by = c("x", "BackgroundRisk")) 
    
    #'  Overlap when background risk is high
    saveOverlap_highrisk <- overlapPlot(spp1_highrisk$sunTime, spp2_highrisk$sunTime, rug = T, 
                                      xscale = NA, xcenter = "noon", linet = c(1, 1), linec = c("red", "blue"), 
                                      linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " (red) & ", name2, " (blue) \ndiel activity when background risk is high")) 
    
    #'  Wrangle from wide to long format
    DensityA_high <- saveOverlap_highrisk[,1:2] %>%
      mutate(PredPrey = "Predator",
             Species = name1,
             BackgroundRisk = "High")
    DensityB_high <- saveOverlap_highrisk[,c(1,3)] %>%
      mutate(PredPrey = "Prey",
             Species = name2,
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
  nboot <- 10000
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
  #'  Save cougar-prey overlap results
  coug_prey_overlap <- list(coug_md_smr_over, coug_md_fall_over, coug_md_wtr_over, coug_md_sprg_over1, coug_md_sprg_over4,
                            coug_elk_smr_over, coug_elk_fall_over, coug_elk_wtr_over, coug_elk_sprg_over,
                            coug_moose_smr_over, coug_moose_fall_over, coug_moose_wtr_over, coug_moose_sprg_over,
                            coug_wtd_smr_over, coug_wtd_fall_over1, coug_wtd_fall_over4, coug_wtd_wtr_over1, coug_wtd_wtr_over4, coug_wtd_sprg_over)
  
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
  #'  Save wolf-prey overlap results
  wolf_prey_overlap <- list(wolf_md_smr_over, wolf_md_fall_over, wolf_md_sprg_over,
                            wolf_elk_smr_over, wolf_moose_smr_over, wolf_moose_fall_over, 
                            wolf_moose_wtr_over, wolf_wtd_smr_over, wolf_wtd_fall_over, wolf_wtd_wtr_over)
  
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
  #'  Save black bear-prey overlap results
  bear_prey_overlap <- list(bear_md_smr_over, bear_md_fall_over1, bear_md_fall_over4, bear_md_sprg_over,
                            bear_elk_smr_over, bear_elk_fall_over, bear_elk_sprg_over,
                            bear_moose_smr_over, bear_moose_fall_over, bear_moose_sprg_over,
                            bear_wtd_smr_over, bear_wtd_fall_over, bear_wtd_sprg_over1, bear_wtd_sprg_over4)
  
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
  #'  Low risk <50 bobcat; High risk >50 bobcat
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
  #'  Save bobcat-prey overlap results
  bob_prey_overlap <- list(bob_md_smr_over, bob_md_fall_over, bob_md_sprg_over,
                           bob_wtd_smr_over, bob_wtd_fall_over1, bob_wtd_fall_over4,
                           bob_wtd_wtr_over, bob_wtd_sprg_over)
  
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
  #'  Save coyote-prey overlap results
  coy_prey_overlap <- list(coy_md_smr_over, coy_md_fall_over, coy_md_wtr_over, coy_md_sprg_over,
                           coy_wtd_smr_over, coy_wtd_fall_over, coy_wtd_wtr_over, coy_wtd_sprg_over)
  
  #'  Save this monster --> list of lists!
  pred_prey_overlap <- list(coug_prey_overlap, wolf_prey_overlap, bear_prey_overlap, bob_prey_overlap, coy_prey_overlap)
  
  save(pred_prey_overlap, file = paste0("./Outputs/Temporal Overlap/PredPrey_LowHi_Overlap_", Sys.Date(), ".RData"))
  
  
  
  #'  --------------------------------------------
  ####  Single species temporal overlap analysis  ####
  #'  --------------------------------------------
  #'  Function to estimate differences in temporal activity for a species at camera  
  #'  sites where background predation risk is high and low based on level of
  #'  habitat complexity and whether predators were ever detected at a site.
  spp_overlap <- function(spp2_dat, risky_sites, name1, name2, nboot, dhat) { #, i
    #'  Create a logical vector (T/F) indicating whether background risk was high 
    #'  or low at each site where spp2 was detected
    risk <- spp2_dat$CameraLocation %in% risky_sites
    spp2_dat <- cbind(spp2_dat, risk) %>% mutate(risk = ifelse(risk == TRUE, "Low", "High"))
    spp2_lowrisk <- filter(spp2_dat, spp2_dat$risk == "Low")
    spp2_highrisk <- filter(spp2_dat, spp2_dat$risk == "High")
    
    #'  Review sample size per species- smaller sample will determine which coefficient
    #'  of overlap estimator to use (delta1 for small samples [<50 detection events], 
    #'  delta4 for larger samples [>75 detection events])
    ndet_spp2.lowrisk <- nrow(spp2_lowrisk)
    ndet_spp2.highrisk <- nrow(spp2_highrisk)
    print(ndet_spp2.lowrisk); print(ndet_spp2.highrisk)
    
    #'  Visualize general temporal activity with density plots
    densityPlot(spp2_dat$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name2, " Daily Activity"))
    
    #'  Visualize temporal overlap
    saveOverlap <- overlapPlot(spp2_lowrisk$sunTime, spp2_highrisk$sunTime, rug = T, 
                               xscale = NA, xcenter = "noon", linet = c(1, 1), linec = c("red", "blue"), 
                               linew = c(2, 2), main = paste0("Overlap Plot of ", name2, " diel activity when \nbackground risk is low (red) & high (blue)")) 
    legend("topleft", c("Low risk", "High risk"), lty=c(1, 1), col=c("red", "blue"), bg = "white", bty = "n", title = paste0("Risk type: ", name1))
    
    #'  Wrangle density data from wide to long format
    Density_low <- saveOverlap[,1:2] %>%
      mutate(BackgroundRisk = "Low")
    colnames(Density_low) <- c("x", "density", "BackgroundRisk")
    Density_high <- saveOverlap[,c(1,3)] %>%
      mutate(BackgroundRisk = "High")
    colnames(Density_high) <- c("x", "density", "BackgroundRisk")
    #'  Bind into single long data set of density estimates to make custom overlap plots
    plotdata <- rbind(Density_low, Density_high) %>%
      mutate(Species = name2,
             RiskType = name1)
       
    #'  Calculate coefficient of overlap
    dhats_spp2.lowhigh <- overlapEst(A = spp2_lowrisk$sunTime, B = spp2_highrisk$sunTime, type = dhat) 
    
    #'  Bootstrap to estimate standard errors
    #'  FYI: smooth = TRUE is default and allows bootstrap to randomly sample from 
    #'  a distribution of times that have a wider range than the original sample
    #'  (see pg. 5 in Overlap package vignette for details) 
    spp2.lowhigh.boot <- bootstrap(spp2_lowrisk$sunTime, spp2_highrisk$sunTime,
                                   nboot, smooth = TRUE, type = dhat)  
    #'  Bootstrap mean will be a little different then detla coefficient due to
    #'  bootstrap bias (BSmean - delta) that needs to be accounted for in 95% CIs
    BSmean <- mean(spp2.lowhigh.boot)
    
    #'  Bootstrap 95% Confidence Intervals
    #'  norm0 uses the standard deviation of bootstrap results to calculate CI (delta +/- 1.96*SDboot)
    #'  basic0 takes the 2.5% and 97.5% percentiles and adjusts based on BS bias (percentile - BSbias)
    #'  If sampling distribution is normal, norm0 and basic0 should be similar;
    #'  if sampling distribution is skewed (i.e., if delta is close to 0 or 1) then
    #'  basic0 is the better estimator
    #'  Using bootCIlogit instead of bootCI so that bias corrections are done on
    #'  the logit scale, then backtransformed. Without this, 95% CIs can fall
    #'  outside (0, 1) interval. See Overlap vignette for more details.
    CI <- bootCIlogit(dhats_spp2.lowhigh, spp2.lowhigh.boot) #dhats_spp1.spp2[i]
    
    #'  Print results
    #'  Effect of spp2 being present
    print("Overlap coefficients when spp2 is present"); print(dhats_spp2.lowhigh)
    print("Bootstrap mean"); print(BSmean)
    print("Bootstrap 95% CI"); print(CI)
    
    #'  Save as a giant list
    overlap_list <- list(dhats_spp2.lowhigh, spp2.lowhigh.boot, CI, ndet_spp2.lowrisk, ndet_spp2.highrisk, plotdata)
    names(overlap_list) <- c("dhats_spp2.lowhigh", "spp2.lowhigh.boot", "CI", "ndet_spp2.lowrisk", "ndet_spp2.highrisk", "overlap.plot.data")
    
    return(overlap_list)
  }
  #'  Estimate temporal overlap for each prey species based on varying levels of 
  #'  background risk - habitat complexity and predator presence (cougar, wolf,
  #'  black bear, bobcat, and coyote)
  nboot <- 100
  
  ####  Mule deer  ####
  #'  Summer activity
  md_smr_hab_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                 name1 = "Habitat complexity", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4") #i = 1
  md_smr_coug_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$Summer_cougar == 0],
                                 name1 = "Cougar detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4")
  md_smr_wolf_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_wolf == 0],
                                  name1 = "Wolf detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_smr_bear_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_bear == 0],
                                  name1 = "Black bear detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_smr_bob_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_bobcat == 0],
                                  name1 = "Bobcat detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_smr_coy_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_coyote == 0],
                                  name1 = "Coyote detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  #'  Fall activity
  md_fall_hab_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                 name1 = "Habitat complexity", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4") #i = 1
  md_fall_coug_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Fall_cougar == 0],
                                  name1 = "Cougar detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_fall_wolf_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Fall_wolf == 0],
                                  name1 = "Wolf detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_fall_bear_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Fall_bear == 0],
                                  name1 = "Black bear detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_fall_bob_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$Fall_bobcat == 0],
                                 name1 = "Bobcat detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4")
  md_fall_coy_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$Fall_coyote == 0],
                                 name1 = "Coyote detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4")
  #'  Winter activity
  md_wtr_hab_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                 name1 = "Habitat complexity", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4") 
  md_wtr_coug_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Winter_cougar == 0],
                                  name1 = "Cougar detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_wtr_wolf_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Winter_wolf == 0],
                                  name1 = "Wolf detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat1") # 46 md detections at high risk cameras
  # md_wtr_bear_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
  #                                 risky_sites = stations_data$CameraLocation[stations_data$bear_det == 0],
  #                                 name1 = "Black bear detected", name2 = "Mule Deer", 
  #                                 nboot = nboot, dhat = "Dhat4")
  md_wtr_bob_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$Winter_bobcat == 0],
                                 name1 = "Bobcat detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4")
  md_wtr_coy_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$Winter_coyote == 0],
                                 name1 = "Coyote detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4")
  #'  Spring activity
  md_sprg_hab_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                 name1 = "Habitat complexity", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4") #i = 1
  md_sprg_coug_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Spring_cougar == 0],
                                  name1 = "Cougar detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_sprg_wolf_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Spring_wolf == 0],
                                  name1 = "Wolf detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat1") # 50 md detections at high risk cameras
  md_sprg_bear_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Spring_bear == 0],
                                  name1 = "Black bear detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  md_sprg_bob_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$Spring_bobcat == 0],
                                 name1 = "Bobcat detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4")
  md_sprg_coy_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                 risky_sites = stations_data$CameraLocation[stations_data$Spring_coyote == 0],
                                 name1 = "Coyote detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4")
  #'  List all mule deer overlap results together
  md_overlap_list <- list(md_smr_hab_over, md_smr_coug_over, md_smr_wolf_over, md_smr_bear_over, md_smr_bob_over, md_smr_coy_over,
                          md_fall_hab_over, md_fall_coug_over, md_fall_wolf_over, md_fall_bear_over, md_fall_bob_over, md_fall_coy_over,
                          md_wtr_hab_over, md_wtr_coug_over, md_wtr_wolf_over, md_wtr_bob_over, md_wtr_coy_over, #md_wtr_bear_over, 
                          md_sprg_hab_over, md_sprg_coug_over, md_sprg_wolf_over, md_sprg_bear_over, md_sprg_bob_over, md_sprg_coy_over)
  
  ####  Elk  ####
  #'  Summer activity
  elk_smr_hab_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                 risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                 name1 = "Habitat complexity", name2 = "Elk", 
                                 nboot = nboot, dhat = "Dhat4")
  elk_smr_coug_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_cougar == 0],
                                  name1 = "Cougar detected", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat4")
  elk_smr_wolf_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_wolf == 0],
                                  name1 = "Wolf detected", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat1") # 36 elk detections at high risk cameras
  elk_smr_bear_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_bear == 0],
                                  name1 = "Black bear detected", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat4")
  #'  Fall activity
  elk_fall_hab_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
                                  risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                  name1 = "Habitat complexity", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat1") # 32 detections at high HCI cameras
  elk_fall_coug_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Fall_cougar == 0],
                                   name1 = "Cougar detected", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat1") # 36 detections at high risk cameras
  # elk_fall_wolf_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
  #                                  risky_sites = stations_data$CameraLocation[stations_data$Fall_wolf == 0],
  #                                  name1 = "Wolf detected", name2 = "Elk", 
  #                                  nboot = nboot, dhat = "Dhat1") # 4 detections at high risk cameras
  elk_fall_bear_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Fall_bear == 0],
                                   name1 = "Black bear detected", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat1") # 45 detections low risk, 47 detections high risk
  #'  Winter activity
  elk_wtr_hab_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
                                 risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                 name1 = "Habitat complexity", name2 = "Elk", 
                                 nboot = nboot, dhat = "Dhat1") # 37/22 detections low/high risk
  elk_wtr_coug_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Winter_cougar == 0],
                                  name1 = "Cougar detected", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat1") #22/37 detections low/high risk
  # elk_wtr_wolf_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
  #                                 risky_sites = stations_data$CameraLocation[stations_data$Winter_wolf == 0],
  #                                 name1 = "Wolf detected", name2 = "Elk",
  #                                 nboot = nboot, dhat = "Dhat1")
  # elk_wtr_bear_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
  #                                 risky_sites = stations_data$CameraLocation[stations_data$Winter_bear == 0],
  #                                 name1 = "Black bear detected", name2 = "Elk", 
  #                                 nboot = nboot, dhat = "Dhat4")
  #'  Spring activity
  elk_sprg_hab_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
                                  risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                  name1 = "Habitat complexity", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat4") # 57 detections high risk cameras
  elk_sprg_coug_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Spring_cougar == 0],
                                   name1 = "Cougar detected", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat1") # 63 detections low risk cameras
  # elk_sprg_wolf_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
  #                                  risky_sites = stations_data$CameraLocation[stations_data$Spring_wolf == 0],
  #                                  name1 = "Wolf detected", name2 = "Elk", 
  #                                  nboot = nboot, dhat = "Dhat1") # 7 detections high risk cameras
  elk_sprg_bear_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Spring_bear == 0],
                                   name1 = "Black bear detected", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat1") # 39 detections high risk cameras
  
  #'  List all elk overlap results together
  elk_overlap_list <- list(elk_smr_hab_over, elk_smr_coug_over, elk_smr_wolf_over, elk_smr_bear_over,
                          elk_fall_hab_over, elk_fall_coug_over, elk_fall_bear_over, #elk_fall_wolf_over, 
                          elk_wtr_hab_over, elk_wtr_coug_over, #elk_wtr_wolf_over,
                          elk_sprg_hab_over, elk_sprg_coug_over, elk_sprg_bear_over) #elk_sprg_wolf_over, 
  
  ####  Moose  ####
  #'  Summer activity
  moose_smr_hab_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                  risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                  name1 = "Habitat complexity", name2 = "Moose", 
                                  nboot = nboot, dhat = "Dhat4") 
  moose_smr_coug_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Summer_cougar == 0],
                                   name1 = "Cougar detected", name2 = "Moose", 
                                   nboot = nboot, dhat = "Dhat4")
  moose_smr_wolf_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Summer_wolf == 0],
                                   name1 = "Wolf detected", name2 = "Moose", 
                                   nboot = nboot, dhat = "Dhat4")
  moose_smr_bear_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Summer_bear == 0],
                                   name1 = "Black bear detected", name2 = "Moose", 
                                   nboot = nboot, dhat = "Dhat1")
  #'  Fall activity
  moose_fall_hab_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                   risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                   name1 = "Habitat complexity", name2 = "Moose", 
                                   nboot = nboot, dhat = "Dhat4") #i = 1
  moose_fall_coug_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Fall_cougar == 0],
                                    name1 = "Cougar detected", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat4")
  moose_fall_wolf_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Fall_wolf == 0],
                                    name1 = "Wolf detected", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat4") # 64 detections high risk cameras
  moose_fall_bear_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Fall_bear == 0],
                                    name1 = "Black bear detected", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat1")
  #'  Winter activity
  moose_wtr_hab_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
                                  risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                  name1 = "Habitat complexity", name2 = "Moose", 
                                  nboot = nboot, dhat = "Dhat4") #i = 1
  moose_wtr_coug_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Winter_cougar == 0],
                                   name1 = "Cougar detected", name2 = "Moose", 
                                   nboot = nboot, dhat = "Dhat4")
  moose_wtr_wolf_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Winter_wolf == 0],
                                   name1 = "Wolf detected", name2 = "Moose",
                                   nboot = nboot, dhat = "Dhat4")
  # moose_wtr_bear_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
  #                                 risky_sites = stations_data$CameraLocation[stations_data$Winter_bear == 0],
  #                                 name1 = "Black bear detected", name2 = "Moose", 
  #                                 nboot = nboot, dhat = "Dhat4")
  #'  Spring activity
  moose_sprg_hab_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                   risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                   name1 = "Habitat complexity", name2 = "Moose", 
                                   nboot = nboot, dhat = "Dhat4") #i = 1
  moose_sprg_coug_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Spring_cougar == 0],
                                    name1 = "Cougar detected", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat1") # 56 detections high risk cameras
  moose_sprg_wolf_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Spring_wolf == 0],
                                    name1 = "Wolf detected", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat1") # 13 detections high risk cameras
  moose_sprg_bear_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Spring_bear == 0],
                                    name1 = "Black bear detected", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat1") # 60 detections high risk cameras
  #'  List all moose overlap results together
  moose_overlap_list <- list(moose_smr_hab_over, moose_smr_coug_over, moose_smr_wolf_over, moose_smr_bear_over,
                           moose_fall_hab_over, moose_fall_coug_over, moose_fall_wolf_over, moose_fall_bear_over,
                           moose_wtr_hab_over, moose_wtr_coug_over, moose_wtr_wolf_over,
                           moose_sprg_hab_over, moose_sprg_coug_over, moose_sprg_wolf_over, moose_sprg_bear_over)
  
  ####  White-tailed deer  ####
  #'  Summer activity
  wtd_smr_hab_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                  name1 = "Habitat complexity", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4") #i = 1
  wtd_smr_coug_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Summer_cougar == 0],
                                   name1 = "Cougar detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4")
  wtd_smr_wolf_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Summer_wolf == 0],
                                   name1 = "Wolf detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4")
  wtd_smr_bear_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Summer_bear == 0],
                                   name1 = "Black bear detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4")
  wtd_smr_bob_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_bobcat == 0],
                                  name1 = "Bobcat detected", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  wtd_smr_coy_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Summer_coyote == 0],
                                  name1 = "Coyote detected", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  #'  Fall activity
  wtd_fall_hab_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                   name1 = "Habitat complexity", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4") #i = 1
  wtd_fall_coug_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Fall_cougar == 0],
                                    name1 = "Cougar detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4")
  wtd_fall_wolf_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Fall_wolf == 0],
                                    name1 = "Wolf detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4")
  wtd_fall_bear_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Fall_bear == 0],
                                    name1 = "Black bear detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4")
  wtd_fall_bob_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Fall_bobcat == 0],
                                   name1 = "Bobcat detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4")
  wtd_fall_coy_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Fall_coyote == 0],
                                   name1 = "Coyote detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4")
  #'  Winter activity
  wtd_wtr_hab_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                  name1 = "Habitat complexity", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4") #i = 1
  wtd_wtr_coug_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Winter_cougar == 0],
                                   name1 = "Cougar detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4")
  wtd_wtr_wolf_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Winter_wolf == 0],
                                   name1 = "Wolf detected", name2 = "White-tailed Deer",
                                   nboot = nboot, dhat = "Dhat4")
  # wtd_wtr_bear_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
  #                                 risky_sites = stations_data$CameraLocation[stations_data$Winter_bear == 0],
  #                                 name1 = "Black bear detected", name2 = "White-tailed Deer", 
  #                                 nboot = nboot, dhat = "Dhat4")
  wtd_wtr_bob_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Winter_bobcat == 0],
                                  name1 = "Bobcat detected", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  wtd_wtr_coy_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Winter_coyote == 0],
                                  name1 = "Coyote detected", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4")
  #'  Spring activity
  wtd_sprg_hab_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$backgroundRisk == "Low"],
                                   name1 = "Habitat complexity", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4") 
  wtd_sprg_coug_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Spring_cougar == 0],
                                    name1 = "Cougar detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4")
  wtd_sprg_wolf_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Spring_wolf == 0],
                                    name1 = "Wolf detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat1") # 38 detections high risk cameras 
  wtd_sprg_bear_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                    risky_sites = stations_data$CameraLocation[stations_data$Spring_bear == 0],
                                    name1 = "Black bear detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4")
  wtd_sprg_bob_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Spring_bobcat == 0],
                                   name1 = "Bobcat detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4")
  wtd_sprg_coy_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                   risky_sites = stations_data$CameraLocation[stations_data$Spring_coyote == 0],
                                   name1 = "Coyote detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4")
  #'  List all white-tailed deer overlap results together
  wtd_overlap_list <- list(wtd_smr_hab_over, wtd_smr_coug_over, wtd_smr_wolf_over, wtd_smr_bear_over, wtd_smr_bob_over, wtd_smr_coy_over,
                           wtd_fall_hab_over, wtd_fall_coug_over, wtd_fall_wolf_over, wtd_fall_bear_over, wtd_fall_bob_over, wtd_fall_coy_over,
                           wtd_wtr_hab_over, wtd_wtr_coug_over, wtd_wtr_wolf_over, wtd_wtr_bob_over, wtd_wtr_coy_over, #wtd_wtr_bear_over, 
                           wtd_sprg_hab_over, wtd_sprg_coug_over, wtd_sprg_wolf_over, wtd_sprg_bear_over, wtd_sprg_bob_over, wtd_sprg_coy_over)
  
  #'  Save all species-specific overlap plots in one giant list
  prey_overlap <- list(md_overlap_list, elk_overlap_list, moose_overlap_list, wtd_overlap_list)
  save(prey_overlap, file = paste0("./Outputs/Temporal Overlap/PreyOnly_LowHi_Overlap_", Sys.Date(), ".RData"))
  
  
  
  
  