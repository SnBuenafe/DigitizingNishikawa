# Loading packages and objects
pacman::p_load("tidyverse", "sf", "raster", "terra", "rnaturalearth", "RColorBrewer", "cmocean", "magrittr", "patchwork")
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true")) # to be used for plotting
source("fSpatPlan_Convert2PacificRobinson.R")

# Load species data
spp <- read.csv("Output/Species_Data.csv")
tow <- read.csv("Output/Tow_Data.csv")

#### Defining preliminaries for plotting and converting to an sf object ####

# Projections
longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
cCRS <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
rob_pacific <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

land <- ne_countries(scale = 'large', returnclass = 'sf') %>% 
  fSpatPlan_Convert2PacificRobinson()

# Plotting necessities (pulled from github.com/valentinitnelav and edited for use)
NE_graticules_rob <- spTransform(NE_graticules, CRSobj = rob_pacific)
NE_box_rob        <- spTransform(NE_box, CRSobj = cCRS)

prj.coord <- rgdal::project(cbind(lbl.Y$lon, lbl.Y$lat), proj = cCRS)
lbl.Y.prj <- cbind(prj.coord, lbl.Y)
names(lbl.Y.prj)[1:2] <- c("X.prj","Y.prj")

prj.coord <- rgdal::project(cbind(lbl.X$lon, lbl.X$lat), proj = rob_pacific)
lbl.X.prj <- cbind(prj.coord, lbl.X)
names(lbl.X.prj)[1:2] <- c("X.prj","Y.prj")


#### Projecting data: Robinson Projection && Make Data Gridded 1x1 ####
# Show all data as gridded data.
make_GriddedData <- function(df, species_name, season_name) {
  longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  rob_pacific <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  df <- df %>% 
    filter(species == species_name, season == season_name)
  
  if (dim(df)[1] != 0) {
    df_crs <- st_crs(longlat)
    
    df_sf <- st_as_sf(df, coords = c("longitude", "latitude"), crs = df_crs)
    
    # Make grid around the limits of the area
    df_poly <- df_sf %>% 
      st_make_grid(cellsize = c(1,1), offset = st_bbox(df_sf)[c("xmin", "ymin")] - 0.5) %>% 
      st_as_sf()
    
    # "Intersect" grids and points (from csv); TRUE if points are contained within the grid cells.
    idx <- st_contains(df_poly, df_sf, sparse = FALSE) %>%
      rowSums() %>% 
      as.logical()
    
    # Project to Pacific-centered Robinson
    df_poly2 <- df_poly[idx,] %>% 
      fSpatPlan_Convert2PacificRobinson() %>% 
      as_tibble()
    
    df_tmp <- df %>% as_tibble() %>% 
      cbind(., df_poly2) %>% 
      st_as_sf(sf_column_name = "x") %>% 
      dplyr::rename(geometry = x)
  }
  else {
    df_tmp <- NULL
  }

  return(df_tmp)
}
make_SpeciesSeasonPlot <- function(df) {
  plot <- ggplot() + geom_sf(data = df, aes(fill = as.factor(abundance)), color = "grey64", size = 0.01) +
    geom_sf(data = land, fill = "grey20", color = NA, size = 0.01) +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
    
    # add graticules projected to Robinson
    geom_path(data=NE_graticules_rob, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25) +
    
    # add graticule labels – latitude and longitude
    geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
    geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
    # the default, ratio = 1 in coord_fixed ensures that one unit on the x-axis is the same length as one unit on the y-axis
    
    scale_fill_manual(values = c("#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"),
                        aesthetics = "fill",
                        name = "Abundance") + theme_classic() + theme(axis.title = element_blank())
  
  return(plot)
}

# Show 0s as points, and blend CPUE category 3 and 4 together. So, instead of having 5 categories, we have 3, just for plotting purposes.
blend_SpeciesSeasonPlot <- function(df, species_name, season_name) {
  
  df_filtered <- df %>% 
    dplyr::filter(species == species_name, season == season_name) %>% 
    dplyr::filter(abundance == 0)
  
  df_crs <- st_crs(longlat)
  
  df_sf <- st_as_sf(df_filtered, coords = c("longitude", "latitude"), crs = df_crs) %>% 
    fSpatPlan_Convert2PacificRobinson()
  
  df_tmp <- df %>% 
    dplyr::filter(abundance > 0) %>% 
    make_GriddedData(., species_name, season_name)
  
  plot <- ggplot()
  if(!is_null(df_tmp)) {
    df_tmp %<>% 
      dplyr::mutate(category = case_when(abundance >= 3 ~ as.integer(3),
                                        TRUE ~ abundance))
    plot <- plot +
      geom_sf(data = df_tmp, aes(fill = as.factor(category)), color = "grey64", size = 0.01)
  }
  
  plot <- plot +
    geom_sf(data = df_sf, size = 0.01, color = "grey64") +
    geom_sf(data = land, fill = "grey20", color = NA, size = 0.01) +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
    
    # add graticules projected to Robinson
    geom_path(data=NE_graticules_rob, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25) +
    
    # add graticule labels – latitude and longitude
    geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
    geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
    # the default, ratio = 1 in coord_fixed ensures that one unit on the x-axis is the same length as one unit on the y-axis
    
    scale_fill_manual(values = c("#c7e9b4", "#41b6c4", "#253494"),
                      aesthetics = "fill",
                      name = "Abundance") + theme_classic() + theme(axis.title = element_blank())
  
  return(plot)
}

# Adding + 0.5 degree to both latitude and longitude
# Because points in the .csv file represent the lower, left of each 1x1 grid cell
# And we want the points to be in the center

## Species Data
spp_tmp <- spp %>% 
  dplyr::mutate(longitude = longitude + 0.5, latitude = latitude + 0.5)

## Towing Data
tow_tmp <- tow %>% 
  dplyr::mutate(longitude = longitude + 0.5, latitude = latitude + 0.5)

#### Plotting Gridded Data ####
# Saves the plots to Figures/ (make sure you have the folder/create the folder!)
season_list <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")
save_plots <- function(species_name, blend = FALSE) { # Default: blend = FALSE; use blend = TRUE if you want to blend the plots (see above for explanation)
  for (i in 1:length(season_list)) {
    if (isTRUE(blend)) {
      blend_SpeciesSeasonPlot(df = spp_tmp, species_name, season_list[i]) %>% 
        ggsave(filename = paste0(species_name, "_", season_list[i], "_Blended.png"), path = "Figures/", width = 21, height = 21, dpi = 300)
    }
    else {
      spp_tmp %>% 
        make_GriddedData(., species_name, season_list[i]) %>% 
        make_SpeciesSeasonPlot(.) %>% 
        ggsave(filename = paste0(species_name, "_", season_list[i], ".png"), path = "Figures/", width = 21, height = 21, dpi = 300)
    }

  }
}
# Skipjack Tuna
save_plots("skipjack-tuna", blend = TRUE)
# Blue Marlin
save_plots("blue-marlin", blend = TRUE)
# Yellowfin Tuna
save_plots("yellowfin-tuna", blend = TRUE)
# Albacore
save_plots("albacore", blend = TRUE)
# Shortbill Spearfish
save_plots("shortbill-spearfish", blend = TRUE)
# Frigate Tuna
save_plots("frigate-tuna", blend = TRUE)
# Bigeye Tuna
save_plots("bigeye-tuna", blend = TRUE)
# Swordfish
save_plots("swordfish", blend = TRUE)
# Striped Marlin
save_plots("striped-marlin", blend = TRUE)
# Sauries
save_plots("sauries", blend = TRUE)
# Sailfish
save_plots("sailfish", blend = TRUE)
# Longfin escolar
save_plots("longfin-escolar", blend = TRUE)
# Bluefin tuna
save_plots("bluefin-tuna", blend = TRUE)
# Little tuna
save_plots("little-tuna", blend = TRUE)
# Southern Bluefin Tuna
save_plots("southern-bluefin-tuna", blend = TRUE)
# Slender tuna
save_plots("slender-tuna", blend = TRUE)
# Bonitos
save_plots("bonitos", blend = TRUE)
# Black marlin
save_plots("black-marlin", blend = TRUE)

#### Map towing effort ####
make_GriddedEffort <- function(df, effort_name, season_name) {
  longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  rob_pacific <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  df <- df %>% 
    filter(category == effort_name, season == season_name)
  
  df_crs <- st_crs(longlat)
  
  df_sf <- st_as_sf(df, coords = c("longitude", "latitude"), crs = df_crs)
  
  df_poly <- df_sf %>% 
    st_make_grid(cellsize = c(1,1), offset = st_bbox(df_sf)[c("xmin", "ymin")] - 0.5) %>% 
    st_as_sf()
  
  idx <- st_contains(df_poly, df_sf, sparse = FALSE) %>%
    rowSums() %>% 
    as.logical()
  
  df_poly2 <- df_poly[idx,] %>% 
    fSpatPlan_Convert2PacificRobinson() %>% 
    as_tibble()
  
  df_tmp <- df %>% as_tibble() %>% 
    cbind(., df_poly2) %>% 
    st_as_sf(sf_column_name = "x") %>% 
    dplyr::rename(geometry = x)
  
  return(df_tmp)
}
make_EffortSeasonPlot <- function(df) {
  plot <- ggplot() + geom_sf(data = df, aes(fill = as.factor(abundance)), color = "grey64", size = 0.01) +
    geom_sf(data = land, fill = "grey20", color = NA, size = 0.01) +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
    # add graticules projected to Robinson
    geom_path(data=NE_graticules_rob, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25) +
    # add graticule labels – latitude and longitude
    geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
    geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
    # the default, ratio = 1 in coord_fixed ensures that one unit on the x-axis is the same length as one unit on the y-axis
    scale_fill_manual(values = c("#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"),
                        aesthetics = "fill",
                        name = "Abundance") + theme_classic() + theme(axis.title = element_blank())
  
  return(plot)
}
#### Plot for tows ####
season_list <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")
save_Effortplots <- function(effort_name) {
  for (i in 1:length(season_list)) {
    tow_tmp %>% 
      make_GriddedEffort(., effort_name, season_list[i]) %>% 
      make_EffortSeasonPlot(.) %>% 
      ggsave(filename = paste0(effort_name, "_", season_list[i], ".png"), path = "Figures/", width = 21, height = 21, dpi = 300)
  }
}

save_Effortplots("tows")
save_Effortplots("volume")

#### Check the order of abundance ####
spp %>% group_by(species) %>% summarize(total_abundance = sum(abundance)) %>% arrange(total_abundance)

#### Make bar plots: proportion of positive samples per season ####
plot_PositiveProportion <- function(df, species_name, season_name) {
  df1 <- df %>% filter(species == species_name, season == season_name) %>% 
    dplyr::select(latitude, abundance) %>% 
    count(latitude, abundance, name = "cases") %>% 
    group_by(latitude) %>% 
    summarize(PositiveProportion = sum(cases[abundance > 0])/sum(cases))
  
  g <- ggplot(df1, aes(x = latitude, y = PositiveProportion))
  
  if(season_name == "jan-mar") {
    g <- g + geom_col(fill = "#D8E8AD", color = "grey20", size = 0.1)
  } else if(season_name == "apr-jun") {
    g <- g + geom_col(fill = "#49A793", color = "grey20", size = 0.1)
  } else if(season_name == "jul-sept") {
    g <- g + geom_col(fill = "#346588", color = "grey20", size = 0.1)
  } else if(season_name == "oct-dec") {
    g <- g + geom_col(fill = "#251B32", color = "grey20", size = 0.1)
  }
  
  plot <- g + coord_flip(ylim = c(0,0.6), xlim = c(-50, 50)) +
   theme_classic() + theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
                             plot.background = element_rect(fill = "transparent", color = NA), # all rectangles
                             panel.grid.major = element_line(color = "grey64", size = 0.1)
                            )
  return(plot)
}

# Skipjack
season_list <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")
file_name <- c("_Spring", "_Summer", "_Autumn", "_Winter")
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "skipjack-tuna", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/SkipjackTuna", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Blue Marlin
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "blue-marlin", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/BlueMarlin", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Yellowfin Tuna
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "yellowfin-tuna", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/YellowfinTuna", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Albacore
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "albacore", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/Albacore", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Shortbill Spearfish
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "shortbill-spearfish", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/ShortbilLSpearish", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Frigate tuna
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "frigate-tuna", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/FrigateTuna", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Bigeye tuna
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "bigeye-tuna", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/BigeyeTuna", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Swordfish
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "swordfish", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/Swordfish", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Striped Marlin
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "striped-marlin", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/StripedMarlin", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Sauries
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "sauries", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/Sauries", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Sailfish
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "sailfish", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/Sailfish", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Longfin Escolar
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "longfin-escolar", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/LongfinEscolar", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Bluefin Tuna
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "bluefin-tuna", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/BluefinTuna", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Little Tuna
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "little-tuna", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/LittleTuna", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Southern Bluefin Tuna
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "southern-bluefin-tuna", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/SouthernBluefinTuna", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Slender Tuna
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "slender-tuna", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/SlenderTuna", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Bonitos
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "bonitos", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/Bonitos", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}

# Black Marlin
for (i in 1:length(season_list)) {
  plot_PositiveProportion(spp, "black-marlin", season_list[i]) %>% ggsave(filename = paste0("Figures/barplots/BlackMarlin", file_name[i], ".png"), bg = "transparent", width = 3, height = 7)
}
