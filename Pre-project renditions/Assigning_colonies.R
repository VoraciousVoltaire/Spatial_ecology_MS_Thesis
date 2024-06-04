# Assigning colonies to specific OSPAR regions 

# Clearing environment
remove(list = ls())

# Checking if the median positions of all individuals lies close to their colony locations

# Loading essential packages----
library(dplyr)
library(rnaturalearth)
library(sf)

# Loading essential packages----
library(dplyr)
library(tidyverse)
library(geosphere)

# Loading data-----
datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
mylocs <- readRDS(paste0(datadir, "/SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"/summaryTable.rds"))

# Merging and curating for NBS (Non-Breeding Season)----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") %>%
  dplyr::mutate(year = year(timestamp))
names(indiv_merged_df)[names(indiv_merged_df) == "ring"] <- "individ_id"

# Creating a new combined dataset which merges tracks from Skalfandi and Langanes (distance between colonies: 130.429 km)
colonies_to_combine <- c("Skjalfandi", "Langanes")
indiv_merged_df$colony[indiv_merged_df$colony %in% colonies_to_combine] <- "Combined"

indiv_merged_df$col_lon[indiv_merged_df$colony == "Combined"] <- mean(unique(indiv_merged_df$col_lon[indiv_merged_df$colony == "Combined"]))
indiv_merged_df$col_lat[indiv_merged_df$colony == "Combined"] <- mean(unique(indiv_merged_df$col_lat[indiv_merged_df$colony == "Combined"]))

col_names <- c("Colony", "Distance_in_km")
output_df <- data.frame(matrix(ncol = length(col_names)))
colnames(output_df) <- col_names

for(i in unique(indiv_merged_df$colony)){
  
  sub <- indiv_merged_df %>% filter(indiv_merged_df$colony == i,) 
  
  p1 <- cbind(unique(sub$col_lon), unique(sub$col_lat))
  p2 <- cbind(median(sub$lon), median(sub$lat))
  
  distance <- (distVincentySphere(p1, p2))/1000 # converting meters to kilometers
  loop_df <- data.frame("Colony" = i, "Distance_in_km" = distance)
  output_df <- rbind(output_df, loop_df)
}
 
View(output_df[-1,]) # Thus, except for Alkefjellet, colony locations are good enough markers;
# Side note: focus on Alkefjellet 

# Now bounding respective regions and checking proximity to each colony

# Creating a new dataset with Alkefjellet's median longlat replacing its col longlat
df <- indiv_merged_df |> dplyr::select(c("colony","col_lon","col_lat"))
df <- df |> group_by(colony) |> summarize(lat = mean(col_lat), lon = mean(col_lon))
sub <- indiv_merged_df |> filter(colony == "Alkefjellet",)
df[df$colony == "Alkefjellet",]$lat <- median(sub$lat)
df[df$colony == "Alkefjellet",]$lon <- median(sub$lon)

# First displaying colonies in Robinson projection
# Defining land----
land <- as(world, "Spatial")
land_sf <- st_as_sf(land)

shp <- ne_countries(type = "countries", 
                    scale = "medium",
                    returnclass = "sf") |>
                    st_transform("ESRI:54030") 

df_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
df_rob <- st_transform(df_sf, crs = crs(shp))
coords_df <- as.data.frame(sf::st_coordinates(df_rob))
coords_df <- cbind(df_rob$colony, coords_df)
names(coords_df)[names(coords_df) == "df_rob$colony"] <- "Colony"

shp |> ggplot() + geom_sf() + 
  coord_sf(expand = F, xlim = c(min(coords_df$X) - 250000, max(coords_df$X) + 250000), 
           ylim = c(min(coords_df$Y) - 250000, max(coords_df$Y) + 250000)) + 
  geom_point(data = coords_df, aes(x = X, y = Y, col = Colony)) 




# For help----
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

# Creating a new dataframe for OSPAR regions and their coordinates

OSPAR_df <- data.frame("Region" = c("North Sea", "Faroe Islands", "Iceland", "Arctic Cananda", "Svalbard"),
                       "lon" = c(3, -6.97, -19.02, -107.07, 16), "lat" = c(56, 61.89, 64.93, 65.81, 78))

col_names <- c("Colony", "Region", "Distance_in_km")
output_df <- data.frame(matrix(ncol = length(col_names)))
colnames(output_df) <- col_names

for(i in unique(df$colony)){ # first for loop starts
  
  sub <- df %>% filter(colony == i)
  lon1 <- sub$lon 
  lat1 <- sub$lat
  p1 <- cbind(lon1, lat1)
  
  for(j in OSPAR_df$Region){ # second for loop starts
    sub <- OSPAR_df %>% filter(Region == j)
    lon2 <- sub$lon
    lat2 <- sub$lat
    p2 <- cbind(lon2, lat2)
    
    distance <- (distVincentySphere(p1, p2))/1000 # converting meters to kilometers
    loop_df <- data.frame("Colony" = i, "Region" = j, "Distance_in_km" = distance)
    output_df <- rbind(output_df, loop_df)
    
  } # second for loop ends
} # first for loop ends

View(final_df)
final_df <- output_df[-1,] %>% group_by(Colony) %>% filter(Distance_in_km == min(Distance_in_km)) 

# Creating a new data frame with the obtained results

ecoqo <- data.frame("Region" = c("English Channel", "Skagerrak", "Faroe Islands", "Iceland", "Svalbard"),
                    "EcoQO" = c(68, 49, 40.5, 27.6, 22.5))

# Creating a new plastic exposure risk score dataset

input <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/exposure_scores_by_population_calculated_on_ind_basis.csv")
to_combine <- c("Skjalfandi", "Langanes")
input$population[input$population %in% to_combine] <- "Combined"

r_1 <- input %>% filter(population == "Inishkea"| population == "Little.Saltee") %>% transmute(pers = mean(population_exposure))
r_2 <- input %>% filter(population == "Jarsteinen") %>% transmute(pers = population_exposure)
r_3 <- input %>% filter(population == "Eynhallow"| population == "Faroe.Islands" | population == "Isle.of.Canna") %>% transmute(pers = mean(population_exposure))
r_4 <- input %>% filter(population == "Combined"| population == "Jan.Mayen") %>% transmute(pers = mean(population_exposure))
r_5 <- input %>% filter(population == "Alkefjellet"| population == "Bjørnøya") %>% transmute(pers = mean(population_exposure))

pers <- data.frame("Region" = c("English Channel", "Skagerrak", "Faroe Islands", "Iceland", "Svalbard"),
                   "pers" = c(unique(r_1$pers), unique(r_2$pers), unique(r_3$pers), unique(r_4$pers), unique(r_5$pers)))

merged_df <- merge(ecoqo, pers, by = "Region")
View(merged_df)

# Correlation test----
correlation <- cor.test(merged_df$EcoQO, merged_df$pers, method = "kendall")
print(correlation)
