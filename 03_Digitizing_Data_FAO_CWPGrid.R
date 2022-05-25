#devtools::install_github("https://github.com/openfigis/RFigisGeo")
# Loading packages and objects
pacman::p_load("rgdal", "RFigisGeo", "sf", "ggplot2", "tidyverse", "magrittr")
source("fSpatPlan_Convert2PacificRobinson.R")

# Intersecting vector files FAO's CWP Grid and Fishing Areas

# Download FAO major fishing areas here: 
FAO <- st_read("FAO_AREAS_CWP/FAO_AREAS_CWP.shp") %>% 
  fSpatPlan_Convert2PacificRobinson() 

FAO_tibble <- FAO %>% 
  as_tibble() %>% dplyr::select(-geometry) %>% # change to dataframe
  rename(FAO_SURFACE = SURFACE)

# See: https://www.fao.org/cwp-on-fishery-statistics/archivedhandbook/general-concepts/major-fishing-areas-general/en/
CWPGrid <- RFigisGeo::createCWPGrid(res = '1deg_x_1deg', xmin = -179, xmax = 180, ymin = -90, ymax = 90) %>% 
  st_as_sf(coords = c("X_COORD", "Y_COORD")) %>% 
  fSpatPlan_Convert2PacificRobinson() 

CWP_tibble <- CWPGrid %>% 
  as_tibble() %>% dplyr::select(-geometry) %>%  # change to dataframe
  mutate(CWP_Code = row_number()) %>% 
  rename(CWP_SURFACE = SURFACE)

FAOGrid <- function(season) {
  grid <- readRDS(paste0("Output/Vector/VectorFile_albacore_", season, ".rds")) %>% 
    mutate(rowNumber = row_number())
  
  # Getting nearest feature
  x <- st_nearest_feature(x = grid, y = FAO)
  y <- st_nearest_feature(x = grid, y = CWPGrid)
  
  grid %<>% 
    mutate(ID = x, CWP_Code = y) %>% 
    as_tibble() %>% 
    left_join(., FAO_tibble, by = "ID") %>% 
    left_join(., CWP_tibble, by = "CWP_Code") %>% 
    dplyr::select(-species, -season, -abundance, -ID, -CWP_Code, -geometry)
  
  return(grid)
}

janMar <- FAOGrid("jan-mar")
aprJun <- FAOGrid("apr-jun")
julSept <- FAOGrid("jul-sept")
octDec <- FAOGrid("oct-dec")

# Load file
season_list <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")
gridList <- list(janMar, aprJun, julSept, octDec)
species_list <- unique(read_csv("Output/Species_Data.csv")$species)

for(j in 1:length(species_list)) {
  for(i in 1:length(season_list)) {
    file <- readRDS(paste0("Output/Vector/VectorFile_", species_list[j] ,"_", season_list[i], ".rds")) %>% 
      as_tibble() %>% dplyr::mutate(rowNumber = row_number()) %>% 
      left_join(., gridList[[i]]) %>% 
      st_as_sf(sf_column_name = "geometry")
    saveRDS(file, paste0("Output/Vector/VectorFile_", species_list[j],"_", season_list[i], ".rds"))
  }
}

# Do this for the tow data as well
effort_list <- c("tows", "volume")
for(j in 1:length(effort_list)) {
  for(i in 1:length(season_list)) {
    file <- readRDS(paste0("Output/Vector/VectorFile_", effort_list[j] ,"_", season_list[i], ".rds")) %>% 
      as_tibble() %>% dplyr::mutate(rowNumber = row_number()) %>% 
      left_join(., gridList[[i]]) %>% 
      st_as_sf(sf_column_name = "geometry")
    saveRDS(file, paste0("Output/Vector/VectorFile_", effort_list[j],"_", season_list[i], ".rds"))
  }
}
