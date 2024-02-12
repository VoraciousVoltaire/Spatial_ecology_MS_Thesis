# Hypothesis 2: assorting distribution points according to the region in which they fall on, then calculating separate home ranges for 'em, then 
# multiplying 'em with the plastic raster and getting a separate pers for each region and then finally perform the correlation test 

# Clearing environment----
remove(list = ls())

# Loading essential packages----
library(dplyr)
library(rnaturalearth)
library(sf)

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

# Loading shapefiles to overlay on the land raster----

PAME_shapefile <- st_read(paste0(datadir,"/PAME/modified_LME.shp"))
plot(PAME_shapefile$geometry)
OSPAR_shapefile <- st_read(paste0(datadir, "/OSPAR/OSPAR_subregions_20160418_3857.shp"))
st_crs(OSPAR_shapefile)

ind_merge_sf <- indiv_merged_df %>% sf::st_as_sf(coords = c('lon', 'lat'), crs = 4326)

# Plotting pooled distribution----
worldmap <- ggplot2::map_data('world')

p1 <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group)) +
  geom_sf(data = ind_merge_sf, aes(colour = colony), size = .1) +
  coord_sf(xlim = c(-90,90), ylim = c(10, 90)) +
  theme_light() +
  ggtitle("Northern fulmar pooled distribution")
p1

# Defining land----
land <- as(world, "Spatial")
land_sf <- st_as_sf(land)

shp <- ne_countries(type = "countries", 
                    scale = "medium",
                    returnclass = "sf") |>
  st_transform("ESRI:54030") 

# Overlaying shapefiles atop distribution 
p2 <- shp |> ggplot() + geom_sf() 
p3 <- p2 + geom_sf(data = ind_merge_sf, aes(colour = colony), size = .1) 

t_PAME <- PAME_shapefile |> st_transform("ESRI:54030")
p4 <- t_PAME |> ggplot() + geom_sf(aes(color = "red"))
p4

# Moment of truth
p5 <- p2 + geom_sf(data = t_PAME, fill = "lightblue") + geom_sf(data = ind_merge_sf, aes(col = colony), size = .1) 
p5

multipolygon <- PAME_shapefile
pame <- multipolygon[1:15,c(3,10)]

# Defining a for loop results data frame----

results_df <- data.frame(matrix(nrow = length(indiv_merged_df$individ_id), ncol = 4))
colnames(results_df) <- c("individ_id", "lon", "lat", "Location")

# for loop for all individuals----
for(i in 1:length(indiv_merged_df$individ_id)){ # First for loop starts
  
# Creating an sf object
sub_lon = indiv_merged_df$lon[i]
sub_lat =  indiv_merged_df$lat[i]
  
for(j in pame$LME_NAME){ # Second for loop starts
  sub_2 <- pame %>% filter(LME_NAME == j) 
  point <- st_sfc(st_point(c(sub_lon, sub_lat)), crs = st_crs(sub_2))
  
  # # Filter out non-polygon geometries
  # polygon_geoms <- sub_2[st_geometry_type(sub_2$geometry) %in% c("POLYGON", "MULTIPOLYGON"), ]
  # class(polygon_geoms$geometry)
  
  # Repair invalid geometries
  sub_2 <- st_make_valid(sub_2)
  
  # Using st_within
  is_inside <- st_within(point, sub_2, sparse = F)
  
  # pol.x = attributes(sub_2$geometry)$bbox[c(1,3)] # bbox x coordinates
  # pol.y = attributes(sub_2$geometry)$bbox[c(2,4)] # bbox y coordinates
  
  if(is_inside == T){ # Second if loop starts
    
    to_bind <- data.frame("individ_id" = indiv_merged_df$individ_id[i], "lon" = sub_lon, "lat" = sub_lat, "Location" = j)
    results_df <- rbind(results_df, to_bind)
    
  } # Second if loop ends
  
  
} # Second for loop ends
} # First for loop ends

View(results_df)
final_df <- merge(indiv_merged_df, results_df, by = c("individ_id", "lon", "lat"))
final_df_clean <- na.omit(final_df)

# What I want is to represent each ind track position to be listed inside one of the PAME_shapefile$LME_NAME; i.e. mutate a separate column 
# with an ifelse condition inside the brackets


# final_df_clean is to be used for calculating kernels based on a monthly basis (just copy its script and replace colony with Location)

