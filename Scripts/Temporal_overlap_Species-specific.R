  #'  ====================================
  #'  Species-specific temporal overlap
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  Sept 2022
  #'  ====================================
  #'  Script to estimate temporal overlap of prey species in response to varying
  #'  levels of background predation risk in a multi-predator system and estimate
  #'  temporal overlap of predator species in response to varying levels of
  #'  habitat complexity. Risk is based on two measures of habitat complexity -  
  #'  an index of terrain ruggedness (TRI) and percentage of forested habitat
  #'  within 250m of each site - and whether each predator species was detected 
  #'  at a camera site (0 = no predator of a given species was detected at a
  #'  camera during given season; 1 = at least 1 predator of a given species was 
  #'  detected at a camera during given season).
  #'  Sites < mean TRI (or % forest) are classified as LOW risk;
  #'  sites >= mean TRI (or % forest) are classified as HIGH risk.
  #'  Sites = 0 predator detection are classified as LOW risk
  #'  sites = 1 predator detection are classified as HIGH risk.
  #'  ====================================
  
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
                    "PercForest", "PercForestMix2", "Canopy_Cov", "Monitoring",
                    "Latitude", "Longitude")) %>%
    mutate(Canopy_Cov = ifelse(is.na(Canopy_Cov), 22, Canopy_Cov)) # fill in missing values based on mean canopy cover (22%)
  
  #'  Make camera coordinates spatial
  cam_locs <- SpatialPoints(cbind(stations_data$Longitude, stations_data$Latitude), 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  #'  Read in mean TRI raster and extract values at each camera site
  mean_tri <- raster("./Shapefiles/mean_tri_250m.tif")
  cam_tri <- raster::extract(mean_tri, cam_locs)
  stations_data$TRI_250m <- cam_tri
  
  #'  Identity which predators were detected at each camera site per season (e.g., predation risk)
  #'  Split up megadetection data by season across years
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
  #'  Add naive occupancy for each predator per season to stations data frame
  stations_data <- cbind(stations_data, bear_det, bob_det, coug_det, coy_det, lynx_det, wolf_det) 
  
  #'  Summary stats on background predation risk
  summary(stations_data)
  hist(stations_data$TRI, breaks = 20)
  abline(v = median(stations_data$TRI), col = "blue")
  abline(v = mean(stations_data$TRI), col = "red")
  hist(stations_data$PercForest, breaks = 20)
  abline(v = median(stations_data$PercForest), col = "blue")
  abline(v = mean(stations_data$PercForest), col = "red")
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
  stations_data <- mutate(stations_data, backgroundRisk_HCI = ifelse(Complexity_index1 < mean(stations_data$Complexity_index1), "Low", "High"))
  stations_data <- mutate(stations_data, backgroundRisk_TRI = ifelse(TRI < mean(stations_data$TRI), "Low", "High"))
  stations_data <- mutate(stations_data, backgroundRisk_For = ifelse(PercForest < mean(stations_data$PercForest), "Low", "High"))
  
  #'  Save for time-btwn-detection analyses
  # write.csv(stations_data, "./Data/cam_stations_hab_complex_data.csv")
  
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
  
  #'  Save for detection data summary and making data available for publication
  smr_dets <- mutate(dets_smr, Season = "Summer")
  fll_dets <- mutate(dets_fall, Season = "Fall")
  wtr_dets <- mutate(dets_wtr, Season = "Winter")
  spg_dets <- mutate(dets_sprg, Season = "Spring")
  seasonal_detections <- as.data.frame(rbind(smr_dets, fll_dets, wtr_dets, spg_dets))
  
  # write.csv(seasonal_detections, paste0("./Outputs/seasonal_detections_", Sys.Date(), ".csv"))
  
  # #'  --------------------------------------------
  ####  Single species temporal overlap analysis  ####
  #'  --------------------------------------------
  #'  Function to estimate differences in temporal activity for a species at camera  
  #'  sites where background predation risk is high and low based on level of
  #'  habitat complexity and whether predators were ever detected at a site.
  spp_overlap <- function(spp2_dat, risky_sites, name1, name2, nboot, dhat, risktype) { 
    #'  Generate variable that is representing high vs low risk
    stations_data$backgroundRisk <- risktype
    risky_sites <- stations_data$CameraLocation[stations_data$backgroundRisk == "Low"]
    
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
  #'  background risk - TRI, % Forest, and predator presence (cougar, wolf,
  #'  black bear, bobcat, and coyote)
  #'  Change predator background risk from 0/1 to low/high
  stations_data <- stations_data %>%
    mutate(across(c(16:39), ~ ifelse(. == 0, "Low", "High")))
  
  #'  Assign number of bootstraps
  nboot <- 10000
  
  ####  Mule deer  ####
  #'  Summer activity
  md_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                 name1 = "TRI", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4", 
                                 risktype = stations_data$backgroundRisk_TRI)
  md_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                 name1 = "PercForest", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4", 
                                 risktype = stations_data$backgroundRisk_For)
  md_smr_coug_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                  name1 = "Cougar detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$Summer_cougar)
  md_smr_wolf_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                  name1 = "Wolf detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$Summer_wolf)
  md_smr_bear_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                  name1 = "Black bear detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$Summer_bear)
  md_smr_bob_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                 name1 = "Bobcat detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4",
                                 risktype = stations_data$Summer_bobcat)
  md_smr_coy_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Mule Deer"),
                                 name1 = "Coyote detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4",
                                 risktype = stations_data$Summer_coyote)
  #'  Fall activity
  md_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                  name1 = "TRI", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$backgroundRisk_TRI) 
  md_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                  name1 = "PercForest", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$backgroundRisk_For) 
  md_fall_coug_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                   name1 = "Cougar detected", name2 = "Mule Deer", 
                                   nboot = nboot, dhat = "Dhat4",
                                   risktype = stations_data$Fall_cougar)
  md_fall_wolf_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                   name1 = "Wolf detected", name2 = "Mule Deer", 
                                   nboot = nboot, dhat = "Dhat4",
                                   risktype = stations_data$Fall_wolf)
  md_fall_bear_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                   name1 = "Black bear detected", name2 = "Mule Deer", 
                                   nboot = nboot, dhat = "Dhat4",
                                   risktype = stations_data$Fall_bear)
  md_fall_bob_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                  name1 = "Bobcat detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$Fall_bobcat)
  md_fall_coy_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Mule Deer"),
                                  name1 = "Coyote detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$Fall_coyote)
  #'  Winter activity
  md_wtr_tri_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                 name1 = "TRI", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4", 
                                 risktype = stations_data$backgroundRisk_TRI)
  md_wtr_for_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                 name1 = "PercForest", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4", 
                                 risktype = stations_data$backgroundRisk_For) #' 57 md detections in high forest habitat
  md_wtr_coug_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                  name1 = "Cougar detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$Winter_cougar)
  md_wtr_wolf_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                  name1 = "Wolf detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat1",
                                  risktype = stations_data$Winter_wolf) # 46 md detections at high risk cameras
  # md_wtr_bear_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
  #                                 name1 = "Black bear detected", name2 = "Mule Deer", 
  #                                 nboot = nboot, dhat = "Dhat4",
  #                                 risktype = stations_data$Winter_bear)
  md_wtr_bob_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                 name1 = "Bobcat detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4",
                                 risktype = stations_data$Winter_bobcat)
  md_wtr_coy_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Mule Deer"),
                                 name1 = "Coyote detected", name2 = "Mule Deer", 
                                 nboot = nboot, dhat = "Dhat4",
                                 risktype = stations_data$Winter_coyote)
  #'  Spring activity
  md_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                  name1 = "TRI", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_TRI)
  md_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                  name1 = "PercForest", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_For)
  md_sprg_coug_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                   name1 = "Cougar detected", name2 = "Mule Deer", 
                                   nboot = nboot, dhat = "Dhat4",
                                   risktype = stations_data$Spring_cougar)
  md_sprg_wolf_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                   name1 = "Wolf detected", name2 = "Mule Deer", 
                                   nboot = nboot, dhat = "Dhat4",
                                   risktype = stations_data$Spring_wolf) # 50 md detections at high risk cameras
  md_sprg_bear_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                   name1 = "Black bear detected", name2 = "Mule Deer", 
                                   nboot = nboot, dhat = "Dhat4",
                                   risktype = stations_data$Spring_bear)
  md_sprg_bob_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                  risky_sites = stations_data$CameraLocation[stations_data$Spring_bobcat == 0],
                                  name1 = "Bobcat detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$Spring_bobcat)
  md_sprg_coy_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Mule Deer"),
                                  name1 = "Coyote detected", name2 = "Mule Deer", 
                                  nboot = nboot, dhat = "Dhat4",
                                  risktype = stations_data$Spring_coyote)
  #'  List all mule deer overlap results together
  md_overlap_list <- list(md_smr_tri_over, md_smr_for_over, md_smr_coug_over, md_smr_wolf_over, md_smr_bear_over, md_smr_bob_over, md_smr_coy_over,
                          md_fall_tri_over, md_fall_for_over, md_fall_coug_over, md_fall_wolf_over, md_fall_bear_over, md_fall_bob_over, md_fall_coy_over,
                          md_wtr_tri_over, md_wtr_for_over, md_wtr_coug_over, md_wtr_wolf_over, md_wtr_bob_over, md_wtr_coy_over, #md_wtr_bear_over, 
                          md_sprg_tri_over, md_sprg_for_over, md_sprg_coug_over, md_sprg_wolf_over, md_sprg_bear_over, md_sprg_bob_over, md_sprg_coy_over)
  
  ####  Elk  ####
  #'  Summer activity
  elk_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                  name1 = "TRI", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_TRI)
  elk_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                  name1 = "PercForest", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_For)
  elk_smr_coug_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                   name1 = "Cougar detected", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Summer_cougar)
  elk_smr_wolf_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                   name1 = "Wolf detected", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$Summer_wolf) # 36 elk detections at high risk cameras
  elk_smr_bear_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Elk"),
                                   name1 = "Black bear detected", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Summer_bear)
  #'  Fall activity
  elk_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
                                   name1 = "TRI", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$backgroundRisk_TRI) # 26 detections at high TRI cameras
  elk_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
                                   name1 = "PercForest", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$backgroundRisk_For) # 44 detections at high forest cameras
  elk_fall_coug_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
                                    name1 = "Cougar detected", name2 = "Elk", 
                                    nboot = nboot, dhat = "Dhat1", 
                                    risktype = stations_data$Fall_cougar) # 40 detections at high risk cameras
  # elk_fall_wolf_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
  #                                  name1 = "Wolf detected", name2 = "Elk",
  #                                  nboot = nboot, dhat = "Dhat1",
  #                                  risktype = stations_data$Fall_wolf) # 6 detections at high risk cameras
  elk_fall_bear_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Elk"),
                                    name1 = "Black bear detected", name2 = "Elk", 
                                    nboot = nboot, dhat = "Dhat1", 
                                    risktype = stations_data$Fall_bear) # 10 detections high risk
  #'  Winter activity
  elk_wtr_tri_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
                                  name1 = "TRI", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat1", 
                                  risktype = stations_data$backgroundRisk_TRI) # 19 detections high TRI
  elk_wtr_for_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
                                  name1 = "PercForest", name2 = "Elk", 
                                  nboot = nboot, dhat = "Dhat1", 
                                  risktype = stations_data$backgroundRisk_For) # 25 detections high Forest
  elk_wtr_coug_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
                                   name1 = "Cougar detected", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$Winter_cougar) #29 detections low risk
  # elk_wtr_wolf_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
  #                                 name1 = "Wolf detected", name2 = "Elk",
  #                                 nboot = nboot, dhat = "Dhat1", 
  #                                 risktype = stations_data$Winter_wolf)
  # elk_wtr_bear_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Elk"),
  #                                 name1 = "Black bear detected", name2 = "Elk",
  #                                 nboot = nboot, dhat = "Dhat4", 
  #                                 risktype = stations_data$Winter_bear)
  #'  Spring activity
  elk_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
                                   name1 = "TRI", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI) # 55 detections high TRI cameras
  elk_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
                                   name1 = "PercForest", name2 = "Elk", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For) # 60 cameras low Forest cameras
  elk_sprg_coug_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
                                    name1 = "Cougar detected", name2 = "Elk", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$Spring_cougar) # 61 detections low risk cameras
  # elk_sprg_wolf_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
  #                                  name1 = "Wolf detected", name2 = "Elk",
  #                                  nboot = nboot, dhat = "Dhat1", 
  #                                  risktype = stations_data$Spring_wolf) # 7 detections high risk cameras
  elk_sprg_bear_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Elk"),
                                    name1 = "Black bear detected", name2 = "Elk", 
                                    nboot = nboot, dhat = "Dhat1", 
                                    risktype = stations_data$Spring_bear) # 35 detections high risk cameras
  
  #'  List all elk overlap results together
  elk_overlap_list <- list(elk_smr_tri_over, elk_smr_for_over, elk_smr_coug_over, elk_smr_wolf_over, elk_smr_bear_over,
                           elk_fall_tri_over, elk_fall_for_over, elk_fall_coug_over, elk_fall_bear_over, #elk_fall_wolf_over, 
                           elk_wtr_tri_over, elk_wtr_for_over, elk_wtr_coug_over, #elk_wtr_wolf_over,
                           elk_sprg_tri_over, elk_sprg_for_over, elk_sprg_coug_over, elk_sprg_bear_over) #elk_sprg_wolf_over, 
  
  ####  Moose  ####
  #'  Summer activity
  moose_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                    name1 = "TRI", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_TRI)
  moose_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                    name1 = "PercForest", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_For)
  moose_smr_coug_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                     name1 = "Cougar detected", name2 = "Moose", 
                                     nboot = nboot, dhat = "Dhat4", 
                                     risktype = stations_data$Summer_cougar)
  moose_smr_wolf_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                     name1 = "Wolf detected", name2 = "Moose", 
                                     nboot = nboot, dhat = "Dhat4", 
                                     risktype = stations_data$Summer_wolf)
  moose_smr_bear_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Moose"),
                                     name1 = "Black bear detected", name2 = "Moose", 
                                     nboot = nboot, dhat = "Dhat1", 
                                     risktype = stations_data$Summer_bear)
  #'  Fall activity
  moose_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                     name1 = "TRI", name2 = "Moose", 
                                     nboot = nboot, dhat = "Dhat4", 
                                     risktype = stations_data$backgroundRisk_TRI)  # 66 high TRI cameras
  moose_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                     name1 = "PercForest", name2 = "Moose", 
                                     nboot = nboot, dhat = "Dhat4", 
                                     risktype = stations_data$backgroundRisk_For) 
  moose_fall_coug_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                      name1 = "Cougar detected", name2 = "Moose", 
                                      nboot = nboot, dhat = "Dhat4", 
                                      risktype = stations_data$Fall_cougar)
  moose_fall_wolf_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                      name1 = "Wolf detected", name2 = "Moose", 
                                      nboot = nboot, dhat = "Dhat4", 
                                      risktype = stations_data$Fall_wolf) # 64 detections high risk cameras
  moose_fall_bear_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Moose"),
                                      name1 = "Black bear detected", name2 = "Moose",
                                      nboot = nboot, dhat = "Dhat1",
                                      risktype = stations_data$Fall_bear)
  #'  Winter activity
  moose_wtr_tri_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
                                    name1 = "TRI", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_TRI)  # 62 high TRI cameras
  moose_wtr_for_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
                                    name1 = "PercForest", name2 = "Moose", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_For) #72 low TRI cameras
  moose_wtr_coug_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
                                     name1 = "Cougar detected", name2 = "Moose", 
                                     nboot = nboot, dhat = "Dhat4", 
                                     risktype = stations_data$Winter_cougar) # 56 high risk cameras
  moose_wtr_wolf_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
                                     name1 = "Wolf detected", name2 = "Moose",
                                     nboot = nboot, dhat = "Dhat4", 
                                     risktype = stations_data$Winter_wolf) # 48 high risk cameras
  # moose_wtr_bear_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Moose"),
  #                                 name1 = "Black bear detected", name2 = "Moose", 
  #                                 nboot = nboot, dhat = "Dhat4", 
  #                                 risktype = stations_data$Winter_bear)
  #'  Spring activity
  moose_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                     name1 = "TRI", name2 = "Moose", 
                                     nboot = nboot, dhat = "Dhat4", 
                                     risktype = stations_data$backgroundRisk_TRI)  # 51 high TRI cameras
  moose_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                     name1 = "PercForest", name2 = "Moose", 
                                     nboot = nboot, dhat = "Dhat4", 
                                     risktype = stations_data$backgroundRisk_For) # 60 low forest cameras
  moose_sprg_coug_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                      name1 = "Cougar detected", name2 = "Moose", 
                                      nboot = nboot, dhat = "Dhat4", 
                                      risktype = stations_data$Spring_cougar) # 56 detections high risk cameras
  moose_sprg_wolf_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                      name1 = "Wolf detected", name2 = "Moose", 
                                      nboot = nboot, dhat = "Dhat1", 
                                      risktype = stations_data$Spring_wolf) # 13 detections high risk cameras
  moose_sprg_bear_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Moose"),
                                      name1 = "Black bear detected", name2 = "Moose", 
                                      nboot = nboot, dhat = "Dhat4", 
                                      risktype = stations_data$Spring_bear) # 60 detections high risk cameras
  #'  List all moose overlap results together
  moose_overlap_list <- list(moose_smr_tri_over, moose_smr_for_over, moose_smr_coug_over, moose_smr_wolf_over, moose_smr_bear_over,
                             moose_fall_tri_over, moose_fall_for_over, moose_fall_coug_over, moose_fall_wolf_over, moose_fall_bear_over,
                             moose_wtr_tri_over, moose_wtr_for_over, moose_wtr_coug_over, moose_wtr_wolf_over,
                             moose_sprg_tri_over, moose_sprg_for_over, moose_sprg_coug_over, moose_sprg_wolf_over, moose_sprg_bear_over)
  
  ####  White-tailed deer  ####
  #'  Summer activity
  wtd_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                  name1 = "TRI", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_TRI)
  wtd_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                  name1 = "PercFores", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_For)
  wtd_smr_coug_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                   name1 = "Cougar detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Summer_cougar)
  wtd_smr_wolf_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                   name1 = "Wolf detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Summer_wolf)
  wtd_smr_bear_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                   name1 = "Black bear detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Summer_bear)
  wtd_smr_bob_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                  name1 = "Bobcat detected", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$Summer_bobcat)
  wtd_smr_coy_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "White-tailed Deer"),
                                  name1 = "Coyote detected", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$Summer_coyote)
  #'  Fall activity
  wtd_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                   name1 = "TRI", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  wtd_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                   name1 = "PercForest", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  wtd_fall_coug_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                    name1 = "Cougar detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$Fall_cougar)
  wtd_fall_wolf_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                    name1 = "Wolf detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$Fall_wolf)
  wtd_fall_bear_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                    name1 = "Black bear detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$Fall_bear)
  wtd_fall_bob_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                   name1 = "Bobcat detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Fall_bobcat)
  wtd_fall_coy_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "White-tailed Deer"),
                                   name1 = "Coyote detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Fall_coyote)
  #'  Winter activity
  wtd_wtr_tri_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                  name1 = "TRI", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_TRI) 
  wtd_wtr_for_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                  name1 = "PercForest", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_For)
  wtd_wtr_coug_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                   name1 = "Cougar detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Winter_cougar)
  wtd_wtr_wolf_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                   name1 = "Wolf detected", name2 = "White-tailed Deer",
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Winter_wolf)
  # wtd_wtr_bear_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
  #                                 name1 = "Black bear detected", name2 = "White-tailed Deer", 
  #                                 nboot = nboot, dhat = "Dhat4", 
  #                                 risktype = stations_data$Winter_bear)
  wtd_wtr_bob_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                  name1 = "Bobcat detected", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$Winter_bobcat)
  wtd_wtr_coy_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "White-tailed Deer"),
                                  name1 = "Coyote detected", name2 = "White-tailed Deer", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$Winter_coyote)
  #'  Spring activity
  wtd_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                   name1 = "TRI", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  wtd_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                   name1 = "PercForest", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  wtd_sprg_coug_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                    name1 = "Cougar detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$Spring_cougar)
  wtd_sprg_wolf_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                    name1 = "Wolf detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat1", 
                                    risktype = stations_data$Spring_wolf) # 38 detections high risk cameras 
  wtd_sprg_bear_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                    name1 = "Black bear detected", name2 = "White-tailed Deer", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$Spring_bear)
  wtd_sprg_bob_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                   name1 = "Bobcat detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Spring_bobcat)
  wtd_sprg_coy_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "White-tailed Deer"),
                                   name1 = "Coyote detected", name2 = "White-tailed Deer", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$Spring_coyote)
  #'  List all white-tailed deer overlap results together
  wtd_overlap_list <- list(wtd_smr_tri_over, wtd_smr_for_over, wtd_smr_coug_over, wtd_smr_wolf_over, wtd_smr_bear_over, wtd_smr_bob_over, wtd_smr_coy_over,
                           wtd_fall_tri_over, wtd_fall_for_over, wtd_fall_coug_over, wtd_fall_wolf_over, wtd_fall_bear_over, wtd_fall_bob_over, wtd_fall_coy_over,
                           wtd_wtr_tri_over, wtd_wtr_for_over, wtd_wtr_coug_over, wtd_wtr_wolf_over, wtd_wtr_bob_over, wtd_wtr_coy_over, #wtd_wtr_bear_over, 
                           wtd_sprg_tri_over, wtd_sprg_for_over, wtd_sprg_coug_over, wtd_sprg_wolf_over, wtd_sprg_bear_over, wtd_sprg_bob_over, wtd_sprg_coy_over)
  
  #'  Save all species-specific overlap plots in one giant list
  prey_overlap <- list(md_overlap_list, elk_overlap_list, moose_overlap_list, wtd_overlap_list)
  save(prey_overlap, file = paste0("./Outputs/Temporal Overlap/PreyOnly_TRI_Forest_Pred_Overlap_", Sys.Date(), ".RData"))
  
  
  ####  Predator Activity Overlap  ####
  ####  Black bear  ####
  #'  Summer activity
  bear_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Black Bear"),
                                 name1 = "TRI", name2 = "Black Bear", 
                                 nboot = nboot, dhat = "Dhat4", 
                                 risktype = stations_data$backgroundRisk_TRI)
  bear_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Black Bear"),
                                 name1 = "PercForest", name2 = "Black Bear", 
                                 nboot = nboot, dhat = "Dhat4", 
                                 risktype = stations_data$backgroundRisk_For)
  #'  Fall activity
  bear_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Black Bear"),
                                   name1 = "TRI", name2 = "Black Bear", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  bear_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Black Bear"),
                                   name1 = "PercForest", name2 = "Black Bear", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  #'  Spring activity
  bear_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Black Bear"),
                                   name1 = "TRI", name2 = "Black Bear", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  bear_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Black Bear"),
                                   name1 = "PercForest", name2 = "Black Bear", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  
  #'  List all black bear overlap results together
  bear_overlap_list <- list(bear_smr_tri_over, bear_smr_for_over, 
                            bear_fall_tri_over, bear_fall_for_over, 
                            bear_sprg_tri_over, bear_sprg_for_over)
  
  ####  Bobcat  ####
  #'  Summer activity
  bob_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Bobcat"),
                                   name1 = "TRI", name2 = "Bobcat", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  bob_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Bobcat"),
                                   name1 = "PercForest", name2 = "Bobcat", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  #'  Fall activity
  bob_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Bobcat"),
                                    name1 = "TRI", name2 = "Bobcat", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_TRI)
  bob_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Bobcat"),
                                    name1 = "PercForest", name2 = "Bobcat", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_For)
  #'  Winter activity
  bob_wtr_tri_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Bobcat"),
                                    name1 = "TRI", name2 = "Bobcat", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_TRI)
  bob_wtr_for_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Bobcat"),
                                    name1 = "PercForest", name2 = "Bobcat", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_For)
  #'  Spring activity
  bob_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Bobcat"),
                                   name1 = "TRI", name2 = "Bobcat", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  bob_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Bobcat"),
                                   name1 = "PercForest", name2 = "Bobcat", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  
  #'  List all bobcat overlap results together
  bob_overlap_list <- list(bob_smr_tri_over, bob_smr_for_over, 
                           bob_fall_tri_over, bob_fall_for_over, 
                           bob_wtr_tri_over, bob_wtr_for_over, 
                           bob_sprg_tri_over, bob_sprg_for_over)
  
  ####  Cougar  ####
  #'  Summer activity
  coug_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Cougar"),
                                  name1 = "TRI", name2 = "Cougar", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_TRI)
  coug_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Cougar"),
                                  name1 = "PercForest", name2 = "Cougar", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_For)
  #'  Fall activity
  coug_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Cougar"),
                                   name1 = "TRI", name2 = "Cougar", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$backgroundRisk_TRI)
  coug_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Cougar"),
                                   name1 = "PercForest", name2 = "Cougar", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  #'  Winter activity
  coug_wtr_tri_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Cougar"),
                                  name1 = "TRI", name2 = "Cougar", 
                                  nboot = nboot, dhat = "Dhat1", 
                                  risktype = stations_data$backgroundRisk_TRI)
  coug_wtr_for_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Cougar"),
                                  name1 = "PercForest", name2 = "Cougar", 
                                  nboot = nboot, dhat = "Dhat4", 
                                  risktype = stations_data$backgroundRisk_For)
  #'  Spring activity
  coug_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Cougar"),
                                   name1 = "TRI", name2 = "Cougar", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  coug_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Cougar"),
                                   name1 = "PercForest", name2 = "Cougar", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  
  #'  List all cougar overlap results together
  coug_overlap_list <- list(coug_smr_tri_over, coug_smr_for_over, 
                            coug_fall_tri_over, coug_fall_for_over, 
                            coug_wtr_tri_over, coug_wtr_for_over, 
                            coug_sprg_tri_over, coug_sprg_for_over)
  
  ####  Coyote  ####
  #'  Summer activity
  coy_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Coyote"),
                                   name1 = "TRI", name2 = "Coyote", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  coy_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Coyote"),
                                   name1 = "PercForest", name2 = "Coyote", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  #'  Fall activity
  coy_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Coyote"),
                                    name1 = "TRI", name2 = "Coyote", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_TRI)
  coy_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Coyote"),
                                    name1 = "PercForest", name2 = "Coyote", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_For)
  #'  Winter activity
  coy_wtr_tri_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Coyote"),
                                   name1 = "TRI", name2 = "Coyote", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_TRI)
  coy_wtr_for_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Coyote"),
                                   name1 = "PercForest", name2 = "Coyote", 
                                   nboot = nboot, dhat = "Dhat4", 
                                   risktype = stations_data$backgroundRisk_For)
  #'  Spring activity
  coy_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Coyote"),
                                    name1 = "TRI", name2 = "Coyote", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_TRI)
  coy_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Coyote"),
                                    name1 = "PercForest", name2 = "Coyote", 
                                    nboot = nboot, dhat = "Dhat4", 
                                    risktype = stations_data$backgroundRisk_For)
  
  #'  List all coyote overlap results together
  coy_overlap_list <- list(coy_smr_tri_over, coy_smr_for_over, 
                           coy_fall_tri_over, coy_fall_for_over, 
                           coy_wtr_tri_over, coy_wtr_for_over, 
                           coy_sprg_tri_over, coy_sprg_for_over)
  
  ####  Wolf  ####
  #'  Summer activity
  wolf_smr_tri_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Wolf"),
                                   name1 = "TRI", name2 = "Wolf", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$backgroundRisk_TRI)
  wolf_smr_for_over <- spp_overlap(spp2_dat = filter(dets_smr, Species == "Wolf"),
                                   name1 = "PercForest", name2 = "Wolf", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$backgroundRisk_For)
  #'  Fall activity
  wolf_fall_tri_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Wolf"),
                                    name1 = "TRI", name2 = "Wolf", 
                                    nboot = nboot, dhat = "Dhat1", 
                                    risktype = stations_data$backgroundRisk_TRI)
  wolf_fall_for_over <- spp_overlap(spp2_dat = filter(dets_fall, Species == "Wolf"),
                                    name1 = "PercForest", name2 = "Wolf", 
                                    nboot = nboot, dhat = "Dhat1", 
                                    risktype = stations_data$backgroundRisk_For)
  #'  Winter activity
  wolf_wtr_tri_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Wolf"),
                                   name1 = "TRI", name2 = "Wolf", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$backgroundRisk_TRI)
  wolf_wtr_for_over <- spp_overlap(spp2_dat = filter(dets_wtr, Species == "Wolf"),
                                   name1 = "PercForest", name2 = "Wolf", 
                                   nboot = nboot, dhat = "Dhat1", 
                                   risktype = stations_data$backgroundRisk_For)
  #'  Spring activity
  wolf_sprg_tri_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Wolf"),
                                    name1 = "TRI", name2 = "Wolf", 
                                    nboot = nboot, dhat = "Dhat1", 
                                    risktype = stations_data$backgroundRisk_TRI)
  wolf_sprg_for_over <- spp_overlap(spp2_dat = filter(dets_sprg, Species == "Wolf"),
                                    name1 = "PercForest", name2 = "Wolf", 
                                    nboot = nboot, dhat = "Dhat1", 
                                    risktype = stations_data$backgroundRisk_For)
  
  
  #'  List all wolf overlap results together
  wolf_overlap_list <- list(wolf_smr_tri_over, wolf_smr_for_over, 
                           wolf_fall_tri_over, wolf_fall_for_over, 
                           wolf_wtr_tri_over, wolf_wtr_for_over, 
                           wolf_sprg_tri_over, wolf_sprg_for_over)
  
  #'  Save all species-specific overlap plots in one giant list
  pred_overlap <- list(bear_overlap_list, bob_overlap_list, coug_overlap_list, coy_overlap_list, wolf_overlap_list)
  save(pred_overlap, file = paste0("./Outputs/Temporal Overlap/PredOnly_TRI_Forest_Overlap_", Sys.Date(), ".RData"))
  
  
  
  #'  Up next: Figures&Tables_OverlapEstiamtes.R script
  