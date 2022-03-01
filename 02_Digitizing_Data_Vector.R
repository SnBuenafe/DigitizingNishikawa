library(tidyverse)
library(magrittr)
library(sf)
source("fSpatPlan_Convert2PacificRobinson.R")

spp <- read_csv("Output/Species_Data.csv")

#### Projecting data: Robinson Projection && Make Data Gridded 1x1 ####

# Adding + 0.5 degree to both latitude and longitude
# Because points in the .csv file represent the lower, left of each 1x1 grid cell
# And we want the points to be in the center

## Species Data
spp_tmp <- spp %>% 
  dplyr::mutate(longitude = longitude + 0.5, latitude = latitude + 0.5)

# Show all data as gridded data.
season_list <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")
make_GriddedData <- function(df, species_name, season_name) {
  longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  rob_pacific <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  df <- df %>% 
    filter(species == species_name, season == season_name)

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
    dplyr::rename(geometry  = x)

  return(df_tmp)
}
save_RObjects <- function(species_name, season_name) {
  x <- make_GriddedData(df = spp_tmp, species_name, season_name)
  
  saveRDS(x, file = paste0("Output/Vector/VectorFile_", species_name, "_", season_name, ".rds"))
  
  return(x)
}

# Save files. If you want to save them as objects in the environment do the following, for example:
# obj <- save_RObjects("skipjack-tuna", "jan-mar")
# Skipjack tuna
for(i in 1:length(season_list)) {
  save_RObjects("skipjack-tuna", season_list[i])
}

# Blue marlin
for(i in 1:length(season_list)) {
  save_RObjects("blue-marlin", season_list[i])
}

# Yellowfin tuna
for(i in 1:length(season_list)) {
  save_RObjects("yellowfin-tuna", season_list[i])
}

# Albacore
for(i in 1:length(season_list)) {
  save_RObjects("albacore", season_list[i])
}

# Shortbill spearfish
for(i in 1:length(season_list)) {
  save_RObjects("shortbill-spearfish", season_list[i])
}

# Frigate tuna
for(i in 1:length(season_list)) {
  save_RObjects("frigate-tuna", season_list[i])
}

# Bigeye tuna
for(i in 1:length(season_list)) {
  save_RObjects("bigeye-tuna", season_list[i])
}

# Swordfish
for(i in 1:length(season_list)) {
  save_RObjects("swordfish", season_list[i])
}

# Striped marlin
for(i in 1:length(season_list)) {
  save_RObjects("striped-marlin", season_list[i])
}

# Sauries
for(i in 1:length(season_list)) {
  save_RObjects("sauries", season_list[i])
}

# Sailfish
for(i in 1:length(season_list)) {
  save_RObjects("sailfish", season_list[i])
}

# Longfin escolar
for(i in 1:length(season_list)) {
  save_RObjects("longfin-escolar", season_list[i])
}

# Bluefin tuna
for(i in 1:length(season_list)) {
  save_RObjects("bluefin-tuna", season_list[i])
}

# Little tuna
for(i in 1:length(season_list)) {
  save_RObjects("little-tuna", season_list[i])
}

# Southern bluefin tuna
for(i in 1:length(season_list)) {
  save_RObjects("southern-bluefin-tuna", season_list[i])
}

# Slender tuna
for(i in 1:length(season_list)) {
  save_RObjects("slender-tuna", season_list[i])
}

# Bonitos
for(i in 1:length(season_list)) {
  save_RObjects("bonitos", season_list[i])
}

# Black marlin
for(i in 1:length(season_list)) {
  save_RObjects("black-marlin", season_list[i])
}
