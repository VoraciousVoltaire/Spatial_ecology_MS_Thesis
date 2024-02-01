# Script to calculate distances between colonies 

# Loading relevant packages----
library(dplyr)
library(sf)
library(sp)
library(tidyverse)
library(geosphere) # for calculating 
library(spData) # for world 

# change

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
mylocs <- readRDS(paste0(datadir, "/SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"/summaryTable.rds"))

# Merging and curating for NBS (Non-Breeding Season)----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"

# Function to calculate distance using Vincenty formula----

vincenty_distance <- function(lat1, lon1, lat2, lon2) {
  p1 <- cbind(lon1, lat1)
  p2 <- cbind(lon2, lat2)
  
}

# For loop for calulating distances between each colony

col_names <- c("Colony_1", "Colony_2", "Distance_in_km")
output_df <- data.frame(matrix(ncol = length(col_names)))
colnames(output_df) <- col_names

for(i in unique(nbs_mylocs$colony)){ # first for loop starts
  sub <- nbs_mylocs %>% filter(colony == i)
  lon1 <- unique(sub$col_lon)
  lat1 <- unique(sub$col_lat)
  p1 <- cbind(lon1, lat1)
  for(j in unique(nbs_mylocs[nbs_mylocs$colony != i,]$colony)){ # second for loop starts
    sub <- nbs_mylocs %>% filter(colony == j)
    lon2 <- unique(sub$col_lon)
    lat2 <- unique(sub$col_lat)
    p2 <- cbind(lon2, lat2)
    
  distance <- (distVincentySphere(p1, p2))/1000 # converting meters to kilometers
  loop_df <- data.frame("Colony_1" = i, "Colony_2" = j, "Distance_in_km" = distance)
  output_df <- rbind(output_df, loop_df)
  
  } # second for loop ends
} # first for loop ends

View(output_df[-1,])

# Plotting colony locations----

# Defining land----
land <- as(world, "Spatial")
land_sf <- st_as_sf(land)
ggplot() + 
  geom_sf(data = land_sf, color = "lightgrey") +
  coord_sf(xlim = c(min(nbs_mylocs$col_lon) - 5, max(nbs_mylocs$col_lon) + 5),
           ylim = c(min(nbs_mylocs$col_lat) - 5, max(nbs_mylocs$col_lat) + 5)) + 
  geom_point(data = nbs_mylocs, aes(x = col_lon, y = col_lat, col = colony), size = 3) + 
  theme_minimal() + theme(panel.grid = element_blank()) +
  ggtitle("Locations of northern fulmar colonies") +
  labs(x = "Latitude", y = "Longitude")

plot(x = unique(nbs_mylocs$col_lon), y = unique(nbs_mylocs$col_lat))
plot(land, col = "lightgrey", alpha = 0.5, add = T)

