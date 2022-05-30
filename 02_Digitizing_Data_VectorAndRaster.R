library(tidyverse)
library(magrittr)
library(sf)
library(terra)
library(stars)
source("fSpatPlan_Convert2PacificRobinson.R")
sf::sf_use_s2(FALSE)

spp <- read_csv("Output/Species_Data.csv")
tow <- read_csv("Output/Tow_Data.csv")

CWP <- st_read("Data/CWP_GRID/cwp-grid-map-1deg_x_1deg.shp") %>% 
  fSpatPlan_Convert2PacificRobinson() %>% 
  dplyr::select(CWP_CODE) %>% 
  mutate(CWP_Code = row_number()) 
  
FAO <- st_read("Data/FAO_AREAS/FAO_AREAS_CWP.shp") %>% 
  fSpatPlan_Convert2PacificRobinson() %>% 
  dplyr::select(F_AREA) %>% 
  mutate(FAO_Code = row_number())

#### Saving as vector files ####
# Vector files are in 1x1 resolution sf objects (polygons)
# And are projected in the Robinson Projection

# Adding + 0.5 degree to both latitude and longitude
# Because points in the .csv file represent the lower, left of each 1x1 grid cell
# And we want the points to be in the center

##### Species Data #####
spp_tmp <- spp %>% 
  dplyr::mutate(latitude = latitude + 0.5) %>% 
  dplyr::mutate(longitude = ifelse(longitude == 180, yes = -179.5,
                       no = longitude + 0.5))

# Show all data as gridded data.
season_list <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")
make_GriddedData <- function(df, species_name, season_name, projected = TRUE) { # Default is projected using Robinson Projection
  longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  rob_pacific <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  df <- df %>% 
    filter(species == species_name, season == season_name)

  df_crs <- st_crs(longlat)
  
  df_sf <- st_as_sf(df, coords = c("longitude", "latitude"), crs = df_crs)
  
  X <- st_coordinates(df_sf)[,1]
  Y <- st_coordinates(df_sf)[,2]
  
  df_sf %<>% bind_cols(., X, Y) %>% 
    rename(longitude = ...5, latitude = ...6) %>% 
    relocate(geometry, .after = latitude)
  
  # Make grid around the limits of the area
  df_poly <- df_sf %>% 
    st_make_grid(cellsize = c(1,1), offset = st_bbox(df_sf)[c("xmin", "ymin")] - 0.5) %>% 
    st_as_sf()
  
  # "Intersect" grids and points (from csv); TRUE if points are contained within the grid cells.
  idx <- st_contains(df_poly, df_sf, sparse = FALSE) %>%
    rowSums() %>% 
    as.logical()
  
  # Project to Pacific-centered Robinson
  df_poly2 <- st_join(x = df_poly[idx,], y = df_sf)
  
  if (isTRUE(projected)){
    df_poly2 %<>% 
      fSpatPlan_Convert2PacificRobinson() %>% 
      as_tibble()
  } else {
    df_poly2 %<>% 
      as_tibble()
  }

  df_final <- df_poly2 %>% 
    st_as_sf(sf_column_name = "geometry")
  
  return(df_final)
}

intersect_FAO <- function(grid) {
  
  FAO_tibble <- FAO %>% 
    st_drop_geometry() %>% 
    as_tibble()
  
  CWP_tibble <- CWP %>% 
    st_drop_geometry() %>% 
    as_tibble()
  
  # Getting nearest feature
  x <- st_nearest_feature(x = grid, y = FAO)
  y <- st_nearest_feature(x = grid, y = CWP)
  
  grid %<>% 
    mutate(FAO_Code = x, CWP_Code = y) %>% 
    as_tibble() %>% 
    left_join(., FAO_tibble, by = "FAO_Code") %>% 
    left_join(., CWP_tibble, by = "CWP_Code") %>% 
    st_as_sf(sf_column_name = "geometry") %>% 
    relocate(geometry, .after = CWP_CODE) %>%  # change the order of columns
    select(-FAO_Code, -CWP_Code) %>% 
    rename(FAO_Major_Fishing_Areas = F_AREA, FAO_CWP_Code = CWP_CODE)
  
  return(grid)
}

save_RObjects <- function(species_name, season_name) {
  x <- make_GriddedData(df = spp_tmp, species_name, season_name) %>% 
    intersect_FAO()
  
  saveRDS(x, file = paste0("Output/Vector/VectorFile_", species_name, "_", season_name, ".rds"))
  
  csv <- x %>% 
    st_drop_geometry()
  
  write_csv(csv, file = paste0("Output/CSV/CSVFile_", species_name, "_", season_name, ".csv"))
}

##### Save files. If you want to save them as objects in the environment do the following, for example: #####
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


#### Saving as raster files ####
save_RasterObj <- function(species_name, season_name, projected = TRUE) {

  if(projected == TRUE) {
    raster <- make_GriddedData(spp_tmp, species_name, season_name, projected) %>% 
      intersect_FAO() %>% 
      mutate(FAO_Major_Fishing_Areas = as.numeric(FAO_Major_Fishing_Areas),
             FAO_CWP_Code = as.numeric(FAO_CWP_Code)) %>% 
      st_rasterize(.) %>% 
      as(., "Raster") %>% 
      terra::rast()
  } else {
    raster <- make_GriddedData(spp_tmp, species_name, season_name, projected) %>% 
      st_rasterize(.) %>% 
      as(., "Raster") %>% 
      terra::rast()
  }

  
  terra::writeRaster(raster, paste0("Output/Raster/RasterFile_", species_name, "_", season_name,".tif"), filetype="GTiff", overwrite=TRUE)
  
  return(raster)
}

# Save files. If you want to save them as objects in the environment do the following, for example:
# obj <- save_RasterObj("skipjack-tuna", "jan-mar")
# Skipjack tuna
season_list <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")
for(i in 1:length(season_list)) {
  save_RasterObj("skipjack-tuna", season_list[i], projected = TRUE)
}

# Blue marlin
for(i in 1:length(season_list)) {
  save_RasterObj("blue-marlin", season_list[i], projected = TRUE)
}

# Yellowfin tuna
for(i in 1:length(season_list)) {
  save_RasterObj("yellowfin-tuna", season_list[i], projected = TRUE)
}

# Albacore
for(i in 1:length(season_list)) {
  save_RasterObj("albacore", season_list[i], projected = TRUE)
}

# Shortbill spearfish
for(i in 1:length(season_list)) {
  save_RasterObj("shortbill-spearfish", season_list[i], projected = TRUE)
}

# Frigate tuna
for(i in 1:length(season_list)) {
  save_RasterObj("frigate-tuna", season_list[i], projected = TRUE)
}

# Bigeye tuna
for(i in 1:length(season_list)) {
  save_RasterObj("bigeye-tuna", season_list[i], projected = TRUE)
}

# Swordfish
for(i in 1:length(season_list)) {
  save_RasterObj("swordfish", season_list[i], projected = TRUE)
}

# Striped marlin
for(i in 1:length(season_list)) {
  save_RasterObj("striped-marlin", season_list[i], projected = TRUE)
}

# Sauries
for(i in 1:length(season_list)) {
  save_RasterObj("sauries", season_list[i], projected = TRUE)
}

# Sailfish
for(i in 1:length(season_list)) {
  save_RasterObj("sailfish", season_list[i], projected = TRUE)
}

# Longfin escolar
for(i in 1:length(season_list)) {
  save_RasterObj("longfin-escolar", season_list[i], projected = TRUE)
}

# Bluefin tuna
for(i in 1:length(season_list)) {
  save_RasterObj("bluefin-tuna", season_list[i], projected = TRUE)
}

# Little tuna
for(i in 1:length(season_list)) {
  save_RasterObj("little-tuna", season_list[i], projected = TRUE)
}

# Southern bluefin tuna
for(i in 1:length(season_list)) {
  save_RasterObj("southern-bluefin-tuna", season_list[i], projected = TRUE)
}

# Slender tuna
for(i in 1:length(season_list)) {
  save_RasterObj("slender-tuna", season_list[i], projected = TRUE)
}

# Bonitos
for(i in 1:length(season_list)) {
  save_RasterObj("bonitos", season_list[i], projected = TRUE)
}

# Black marlin
for(i in 1:length(season_list)) {
  save_RasterObj("black-marlin", season_list[i], projected = TRUE)
}


##### Tow Data #####
tow_tmp <- tow %>% 
  dplyr::mutate(latitude = latitude + 0.5) %>% 
  dplyr::mutate(longitude = ifelse(longitude == 180, yes = -179.5,
                                   no = longitude + 0.5))
# Vector
make_GriddedData <- function(df, category_name, season_name, projected = TRUE) { # Default is projected using Robinson Projection
  longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  rob_pacific <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  df <- df %>% 
    filter(category == category_name, season == season_name)
  
  df_crs <- st_crs(longlat)
  
  df_sf <- st_as_sf(df, coords = c("longitude", "latitude"), crs = df_crs)
  
  X <- st_coordinates(df_sf)[,1]
  Y <- st_coordinates(df_sf)[,2]
  
  df_sf %<>% bind_cols(., X, Y) %>% 
    rename(longitude = ...5, latitude = ...6) %>% 
    relocate(geometry, .after = latitude)
  
  # Make grid around the limits of the area
  df_poly <- df_sf %>% 
    st_make_grid(cellsize = c(1,1), offset = st_bbox(df_sf)[c("xmin", "ymin")] - 0.5) %>% 
    st_as_sf()
  
  # "Intersect" grids and points (from csv); TRUE if points are contained within the grid cells.
  idx <- st_contains(df_poly, df_sf, sparse = FALSE) %>%
    rowSums() %>% 
    as.logical()
  
  # Project to Pacific-centered Robinson
  df_poly2 <- st_join(x = df_poly[idx,], y = df_sf)
  
  if (isTRUE(projected)){
    df_poly2 %<>% 
      fSpatPlan_Convert2PacificRobinson() %>% 
      as_tibble()
  } else {
    df_poly2 %<>% 
      as_tibble()
  }
  
  df_final <- df_poly2 %>% 
    st_as_sf(sf_column_name = "geometry")
  return(df_final)
}

save_RObjects <- function(category_name, season_name) {
  x <- make_GriddedData(df = tow_tmp, category_name, season_name) %>% 
    intersect_FAO()
  
  saveRDS(x, file = paste0("Output/Vector/VectorFile_", category_name, "_", season_name, ".rds"))
  
  csv <- x %>% 
    st_drop_geometry()
  
  write_csv(csv, file = paste0("Output/CSV/CSVFile_", category_name, "_", season_name, ".csv"))
}
# Tows
for(i in 1:length(season_list)) {
  save_RObjects("tows", season_list[i])
}
# Volume
for(i in 1:length(season_list)) {
  save_RObjects("volume", season_list[i])
}
# Raster
save_RasterObj <- function(category_name, season_name, projected = TRUE) {
  if(projected == TRUE) {
    raster <- make_GriddedData(tow_tmp, category_name, season_name, projected) %>% 
      intersect_FAO() %>% 
      mutate(FAO_Major_Fishing_Areas = as.numeric(FAO_Major_Fishing_Areas),
             FAO_CWP_Code = as.numeric(FAO_CWP_Code)) %>% 
      st_rasterize(.) %>% 
      as(., "Raster") %>% 
      terra::rast()
  } else{
    raster <- make_GriddedData(tow_tmp, category_name, season_name, projected) %>% 
      st_rasterize(.) %>% 
      as(., "Raster") %>% 
      terra::rast()
  }
  
  terra::writeRaster(raster, paste0("Output/Raster/RasterFile_", category_name, "_", season_name,".tif"), filetype="GTiff", overwrite=TRUE)
  
  return(raster)
}
# Tows
for(i in 1:length(season_list)) {
  save_RasterObj("tows", season_list[i], projected = TRUE)
}
# Volume
for(i in 1:length(season_list)) {
  save_RasterObj("volume", season_list[i], projected = TRUE)
}
