library(tidyverse)
library(magrittr)
library(doParallel)

#### Digitize Spp Data and Tow Data ####

spp <- c("albacore", "bigeye-tuna", "black-marlin", "blue-marlin", "bluefin-tuna", "bonitos", "frigate-tuna", "little-tuna", "longfin-escolar", "sailfish", "sauries", "shortbill-spearfish", "skipjack-tuna", "slender-tuna", "southern-bluefin-tuna", "striped-marlin", "swordfish", "yellowfin-tuna")
season <- c("jan-mar", "apr-jun", "jul-sept", "oct-dec")

dir <- "RawData/"
file_list <- list.files(dir) 

# make sequences for longitude and latitude:
l1 <- seq(from = 100, to = -180, by = -1)
l2 <- seq(from = 179, to = 101, by = -1)
longitude <- c("coordinates", l1, l2)
latitude <- seq(from = 50, to = -50, by = -1)

ncores <- detectCores() - 1 
cl <- makeCluster(ncores)
registerDoParallel(cl)

csv <- vector("list", length = (length(spp) * length(season)))

per_species <- foreach(i = 1:length(file_list), .packages = c('tidyverse')) %dopar% {
  
  # Create empty tibble
  csv_blank <- tibble(species = character(),
                      season = character(),
                      latitude = numeric(),
                      longitude = numeric(),
                      abundance = numeric())
  
  file <- paste0(dir, file_list[i])
  file_name <- str_replace(file_list[i], ".csv", "")
  species_name <- (str_split(file_name, pattern = "_") %>% unlist())[1]
  season_name <- (str_split(file_name, pattern = "_") %>% unlist())[2]
  
  df <- read.csv(file) %>% 
    as_tibble() %>% 
    select(-X) # removing first column
    
  # make sure that dimensions are 361 (longitude) x 101 (latitude)
  df_new <- df[1:length(latitude), 1:length(longitude)] %>% 
    cbind(latitude, .) %>% 
    as_tibble()
  colnames(df_new) <- longitude
  
  # Populate .csv tibble
  for(j in 2:ncol(df_new)){
    for(k in 1:nrow(df_new)){
      if(!is.na(df_new[k, j])){
        csv_blank <- csv_blank %>% 
          add_row(species = species_name,
                  season = season_name,
                  latitude = df_new$coordinates[k],
                  longitude = as.numeric(colnames(df_new)[j]),
                  abundance = df_new[k, j] %>% pull())
      }
    }
  }
  
  csv[[i]] <- csv_blank
}
stopCluster(cl)

temp <- do.call(rbind, per_species)

# Do small change on Bonitos
temp <- temp %>%  
  add_row(species = "bonitos", season = "jul-sept", longitude = -106, latitude = -8, abundance = 0)

spp_data <- temp %>% filter(species %in% spp)
write_csv(spp_data, "Output/Species_Data.csv")

tow_data <- temp %>% filter(!species %in% spp) %>% 
  rename(category = species, effort = abundance)
write_csv(tow_data, "Output/Tow_Data.csv")
