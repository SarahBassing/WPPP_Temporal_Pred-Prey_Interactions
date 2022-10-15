  #'  ======================================
  #'  Map of study area and camera locations
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  October 2022
  #'  ======================================
  
  #'  Load libraries
  library(ggplot2)
  library(ggspatial)
  library(cowplot)
  library(sf)
  library(raster)
  library(tidyverse)
  library(viridis)
  
  
  #'  Get some basics pulled together to be used across most figures
  #'  -----------------------------
  ####  Spatial Data for Mapping  ####
  #'  -----------------------------
  #'  Define projections
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Read in study area data and reproject
  WA <- st_read("./Shapefiles/Washington_State_Boundary/WA_State_Geospatial_Open_Data", layer = "WA_State_Boundary") %>%
    st_transform(crs = wgs84)
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = wgs84)
  OK_SA$NAME <- "Okanogan"
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = wgs84)
  NE_SA$NAME <- "Northeast"
 
  projection(WA)
  projection(OK_SA)
  extent(OK_SA)
  extent(NE_SA)
  
  #'  Centroid of each polygon
  st_centroid(OK_SA)
  st_centroid(NE_SA)
  
  #'  Reduced DEM raster resolution
  tri_low <- raster("./Shapefiles/WA DEM rasters/WPPP_TRI_low.tif")
  tri_p_low <- rasterToPoints(tri_low)
  tri_p_df <- as.data.frame(tri_p_low)
  colnames(tri_p_df) <- c("x", "y", "value")
  projection(tri_low)
  
  ####  Map study area and camera locations  ####
  #'  ========================================
  #'  Read in camera locations
  station_covs <- read.csv("./Data/cam_stations_data.csv") %>%
    dplyr::select(-X)
  
  #'  Make camera location data spatial and reproject to study area projection
  cams <- st_as_sf(station_covs, coords = c("Longitude", "Latitude"), crs = wgs84) #%>%
    #st_transform(cams, crs = sa_proj)
  cams <- mutate(cams, Year2 = ifelse(Year == "Year1", "2018-2019", "2019-2020"),
                       Year2 = ifelse(Year == "Year3", "2020-2021", Year2 )) 
  
  #'  City locations
  city <- c("Chewelah, \n   WA", "Winthrop, \n   WA")
  x <- c(-117.7193, -120.1096)
  y <- c(48.28302, 48.42966)
  city_df <- as.data.frame(cbind(city, x, y))
  city_sf <- st_as_sf(city_df, coords = c("x", "y"), crs = wgs84) #%>%
    #st_transform(crs = sa_proj)
  
  
  #'  Plot state of WA with study areas
  #'  https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
  WA_SA_map <- ggplot() + 
    geom_sf(data = WA, fill = "gray95") +
    #' #'  Label map of WA with "Washington State"
    #' geom_sf_text(data = WA, aes(label = JURISDIC_3, hjust = 0.8, vjust = 3), size = 12) +
    geom_sf(data = OK_SA, fill = "grey25", color = "grey20") +
    geom_sf_text(data = OK_SA, aes(label = NAME, hjust = 1.15, vjust = 1.6), size = 4.5) +
    geom_sf(data = NE_SA, fill = "grey25", color = "grey20") +
    geom_sf_text(data = NE_SA, aes(label = NAME, hjust = 0.55, vjust = 3), size = 4.5) +
    #'  Get rid of lines and axis labels
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
          #'  No margins around figure
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  plot(WA_SA_map)
  
  #'  Plot study areas against DEM with camera locations
  cam_SA_map <- ggplot() +
    geom_raster(data = tri_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) + 
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey80", high = "gray10") + #tan4
    #'  Add study area outlines and label with their names
    geom_sf(data = OK_SA, fill = NA, color="black", size = 0.80) +
    #'  Note the hjust & vjust need to change based on font size and coord_sf
    #'  DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
    # geom_sf_text(data = OK_SA, aes(label = NAME, hjust = 1.3, vjust = 7), size = 7) + #vjust = -6.90
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    # geom_sf_text(data = NE_SA) + #, aes(label = NAME, hjust = 1.3, vjust = -5.5), size = 7) + #vjust = -6.5
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = cams, aes(color = Year2), shape = 19, size = 2.5) +
    #'  Change camera data aesthetics (make sure it's colorblind friendly)
    scale_color_viridis_d(alpha = 0.6) +
    # geom_sf(data = city_sf, col = "black", shape = 1, size = 2.5) +
    # geom_sf_text(data = city_sf, aes(label = city, hjust = -0.05, vjust = 1), size = 3) +
    labs(colour = "Camera\nlocations") +
    #'  Constrain plot to two study areas plus some room on the side & bottom
    coord_sf(xlim = c(-116.8, -121.1), ylim = c(47.3, 49.0), expand = FALSE) +
    #'  Get rid of lines and axis names
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16)) +
    labs(x = "Longitude", y = "Latitude") +
    #'  Add north arrow
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering) +
    #'  Add scale bar (be sure to double check the scale)
    annotation_scale(location = "bl", width_hint = 0.5)
  plot(cam_SA_map)
  
  #'  Build plot with map of study areas and inset map of WA
  #'  https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  #'  Requires "cowplot" package
  #'  Don't use png or other calls to save while plotting- formatting gets messed up
  #'  Use export option in Plot window and formatting holds
  # tiff(file = "./Outputs/Figures/Maps/StudyAreas_Cameras1820_05.31.22.tiff",
  #     width = 1000, height = 691) 
  tiff(file = "./Outputs/Figures/StudyAreaMap_Cameras_2018-2020.tiff",
       units="in", width=11, height=6.5, res=800, compression = 'lzw') 
  StudyArea_Map <- ggdraw(cam_SA_map) +
    draw_plot(
      {
        WA_SA_map #+
        #'  Label map of WA with "Washington State"
        #'  hjust & vjust will depend on inset map's width/height specified below
        # geom_sf_text(data = WA, aes(label = JURISDIC_3, hjust = 0.5, vjust = 1), size = 5)  
      },
      #'  Distance along a (0,1) x-axis to draw the left edge of the plot
      x = 0.6,
      #'  Distance along a (0,1) y-axis to draw the bottom edge of the plot
      y = 0.20,
      #'  Width & height of the plot expressed as proportion of the entire ggdraw object
      #'  THIS DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
      width = 0.25,
      height = 0.25)
  plot(StudyArea_Map)
  dev.off()
