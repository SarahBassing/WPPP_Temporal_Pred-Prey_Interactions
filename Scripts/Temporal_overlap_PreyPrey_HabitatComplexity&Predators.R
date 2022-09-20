  #'  ==================================================================
  #'  Prey-Prey temporal overlap - habitat complexity & predator detections
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  Sept 2022
  #'  ==================================================================
  #'  Script to estimate temporal overlap of prey species in response to varying
  #'  levels of background predation risk in a multi-predator system. Risk is
  #'  based on two measures of habitat complexity: an index of terrain ruggedness 
  #'  (TRI) and percentage of forested habitat within 250m of each site and
  #'  whether each predator species was detected at a camera site (0 = no predator
  #'  of a given species was detected at camera during given season; 1 = at least 
  #'  1 predator of a given species was detected at camera during given season).
  #'  Sites < mean TRI (or % forest) are classified as LOW risk;
  #'  sites >= mean TRI (or % forest) are classified as HIGH risk.
  #'  Sites = 0 predator detection are classified as LOW risk
  #'  sites = 1 predator detection are classified as HIGH risk.
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
  #'  Add naive occupancy for each predator per season to statations dataframe
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
  
  #' #'  Save for time-btwn-detection analyses
  write.csv(stations_data, "./Data/cam_stations_hab_complex_data.csv")
  
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