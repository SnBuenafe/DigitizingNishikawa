# Loading packages and objects
pacman::p_load("tidyverse", "sf", "raster", "terra", "rnaturalearth", "RColorBrewer", "cmocean", "magrittr")
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
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

# Plotting necessities
NE_graticules_rob <- spTransform(NE_graticules, CRSobj = rob_pacific)
NE_box_rob        <- spTransform(NE_box, CRSobj = cCRS)

prj.coord <- rgdal::project(cbind(lbl.Y$lon, lbl.Y$lat), proj = cCRS)
lbl.Y.prj <- cbind(prj.coord, lbl.Y)
names(lbl.Y.prj)[1:2] <- c("X.prj","Y.prj")

prj.coord <- rgdal::project(cbind(lbl.X$lon, lbl.X$lat), proj = rob_pacific)
lbl.X.prj <- cbind(prj.coord, lbl.X)
names(lbl.X.prj)[1:2] <- c("X.prj","Y.prj")



#### Creating barplots (total abundance per season) ####
df <- spp %>% group_by(species, season) %>% 
  summarize(total_abundance = sum(abundance))

df$season <- factor(df$season, levels = c("jan-mar", "apr-jun", "jul-sept", "oct-dec"))

StackedBar <- ggplot(df, aes(fill=species, y=total_abundance, x=season)) + 
  geom_bar(position="stack", stat="identity") +
  cmocean::scale_fill_cmocean(name = "deep",
                              discrete = TRUE) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"))
ggsave("03_StackedBarplot.png",
       plot = StackedBar, width = 21, height = 21, dpi = 300,
       path = "Figures/")

# To get top 8 most abundant species
df %>% group_by(species) %>% summarize(total = sum(total_abundance)) %>% arrange(total) %>% print(n = Inf)
#### Projecting data: Robinson Projection && Make Data Gridded 1x1 ####
make_GriddedData <- function(df, species_name, season_name) {
  longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  rob_pacific <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  df <- df %>% 
    filter(species == species_name, season == season_name)
  
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
season_list <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")
save_plots <- function(species_name) {
  for (i in 1:length(season_list)) {
    spp_tmp %>% 
      make_GriddedData(., species_name, season_list[i]) %>% 
      make_SpeciesSeasonPlot(.) %>% 
      ggsave(filename = paste0(species_name, "_", season_list[i], ".png"), path = "Figures/", width = 21, height = 21, dpi = 300)
  }
}

# Skipjack Tuna
save_plots("skipjack-tuna")
# Blue Marlin
save_plots("black-marlin")
# Yellowfin Tuna
save_plots("yellowfin-tuna")
# Albacore
save_plots("albacore")
# Shortbill Spearfish
save_plots("shortbill-spearfish")
# Frigate Tuna
save_plots("frigate-tuna")
# Bigeye Tuna
save_plots("bigeye-tuna")
# Swordfish
save_plots("swordfish")
# Striped Marlin
save_plots("striped-marlin")
# Sauries
save_plots("sauries")
# Sailfish
save_plots("sailfish")
# Longfin escolar
save_plots("longfin-escolar")
# Bluefin tuna
save_plots("bluefin-tuna")
# Little tuna
save_plots("little-tuna")
# Southern Bluefin Tuna
save_plots("southern-bluefin-tuna")
# Slender tuna
save_plots("slender-tuna")
# Bonitos
save_plots("bonitos")
# Black marlin
save_plots("black-marlin")

#### Map towing effort ####
make_GriddedEffort <- function(df, effort_name, season_name) {
  longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  rob_pacific <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  df <- df %>% 
    filter(effort == effort_name, season == season_name)
  
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
    scale_colour_manual(values = c("#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"),
                        aesthetics = "color",
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

#### Check the order of abundance
spp %>% group_by(species) %>% summarize(total_abundance = sum(abundance)) %>% arrange(total_abundance)
