# Images for report----

# Loading essential packages----
library(spData)
library(tidyverse)
library(ggplot2)
library(leaflet)
library(rnaturalearth)
library(sf)
library(sp)
library(raster)
library(mapproj)
library(rayshader)
library(maps)
library(RColorBrewer)
library(viridis)
library(cowplot)

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
mylocs <- readRDS(paste0(datadir, "/SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"/summaryTable.rds"))
indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"
all_data <- nbs_mylocs
df <- all_data[!is.na(all_data$timestamp), ]
df <- df %>% dplyr::mutate(month = month(timestamp))
months <- sort(unique(df$month))
df_mod <- df %>% mutate(tracking_year = ifelse(month < 7, year - 1, year))
df_mod_sf <- df_mod %>% 
  sf::st_as_sf(coords = c("col_lon", "col_lat"), crs = 4326) 

# Changing to different projection system
df_mod_robinson <- st_transform(df_mod_sf, crs = "+proj=robin")

# Check the new projection
st_crs(df_mod_robinson)


# Defining land
land <- as(world, "Spatial")
world_robinson <- st_transform(world, crs = "+proj=robin")

# Check the new projection
st_crs(world_robinson)


# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# 1. Distribution of colonies overlaid on top of five specific PAME regions----

# Plotting

# 1.1 Leaflet plot

# Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Convert world map data to WGS84 projection
world <- st_transform(world, "+proj=longlat +datum=WGS84")

# Load water bodies data
water <- ne_download(category = "physical", type = "ocean", returnclass = "sf")

# Convert water bodies data to WGS84 projection
water <- st_transform(water, "+proj=longlat +datum=WGS84")

# Data
colonies <- data.frame(
  colony = c("Alkefjellet", "Bjørnøya", "Eynhallow", "Faroe Islands", "Inishkea", "Isle of Canna", 
             "Jan Mayen", "Jarsteinen", "Langanes", "Little Saltee", "Skjalfandi"),
  lat = c(79.6, 74.5, 59.1, 62.0, 54.1, 57.1, 70.9, 59.2, 66.4, 52.1, 66.0),
  lon = c(18.5, 19.0, -3.12, -6.80, -10.2, -6.53, -8.72, 5.17, -14.6, -6.62, -17.4)
)
library(sf)

# Convert colony locations to sf object
colonies_sf <- st_as_sf(colonies, coords = c("lon", "lat"), crs = 4326)

# Transform to Robinson projection
colonies_robinson <- st_transform(colonies_sf, crs = "+proj=robin")

# Create leaflet map
m <- leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  setView(lng = -10, lat = 60, zoom = 4) %>%  # Centered over the northeast Atlantic
  addPolygons(data = world, fillOpacity = 0, color = "black", weight = 1) %>%  # Add country boundaries
  addPolygons(data = water, fillColor = "lightblue", fillOpacity = 0.5, color = "transparent") %>%  # Add water bodies in lightblue
  addMarkers(data = colonies, ~lon, ~lat, label = ~colony, options = markerOptions(iconSize = c(10, 10)))  # Add markers with colony names and adjust marker size

# Display the map
m

# Mimicking Clark et al. (2023) script----

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data/") # directory where plastics raster is stored

land <- as(world, "Spatial")

# Distributions of all birds 
bird_dist <- raster::raster("08_all_species_distribution.tif")
b <- bird_dist
b[is.na(b)] <- 0 
b_sum1 <- b/sum(raster::getValues(b))
b_sum1[is.na(bird_dist)] <- NA
b_sum1[b_sum1 == 0] <- NA

plastics <- raster::raster("00_PlasticsRaster.tif")

## rescale to 1
plastics2 <- plastics
plastics2[is.na(plastics2)] <- 0 
p_sum1    <- plastics2/sum(raster::getValues(plastics2))
p_sum1[is.na(plastics)] <- NA

yelblus <- c(brewer.pal(n = 9, name = "YlGnBu"),"#00172e")
col_birds <- c(colorRampPalette(yelblus)(1000))

# define Robinson projection
proj <- "+proj=robin"

# project shapefile and raster to Robinson projection
land_sf <- sf::st_as_sf(land)
world_prj <- land_sf %>% sf::st_transform(proj)
plot(world_prj$geometry)
b_dens_proj <- raster::projectRaster(b_sum1, crs = proj)

# convert raster to dataframe for ggplot
bat_df <- raster::rasterToPoints(b_dens_proj) %>% as_tibble()

# if you can simply mask it then create a polygon that is a square larger than the map area
# then cookie cut out the map region

# define a box to clip shapefile to - make slightly bigger than required to avoid edge clip effects
CP <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90)) %>%
  sf::st_as_sfc() %>%
  sf::st_segmentize(1) %>% st_set_crs(4326)
CP_prj <- CP %>% sf::st_transform(crs = proj)

# get the bounding box in transformed coordinates and expand by 10%
xlim <- sf::st_bbox(CP_prj)[c("xmin", "xmax")]*1.2
ylim <- sf::st_bbox(CP_prj)[c("ymin", "ymax")]*1.2

# turn into enclosing rectangle
encl_rect <- 
  list(
    cbind(
      c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]), 
      c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
    )
  ) %>%
  sf::st_polygon() %>%
  sf::st_sfc(crs = proj)

# calculate the area outside the earth outline as the difference
# between the enclosing rectangle and the earth outline

cookie <- sf::st_difference(encl_rect, CP_prj)

p1 <- ggplot() +
  theme_minimal() +
  theme(legend.position = "none")+
  geom_sf(aes(), colour = NA, fill = "grey75", 
          data = world_prj) +
  geom_sf(aes(),
             data = df_mod_robinson, colour="red", pch=18, size=1)+
  geom_sf(aes(), fill = "white", color = "black", 
          data = cookie, size = .5) +
  xlab("") + ylab("");p1

p1 <- ggplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_sf(aes(), colour = NA, fill = "grey75", 
          data = world_prj) +
  geom_sf(aes(), data = df_mod_robinson, colour = "red", pch = 18, size = 1) +
  geom_sf(aes(), fill = "white", color = "black", 
          data = cookie, size = 0.5) +
  xlab("") + ylab("") +
  coord_sf(xlim = c(-20, 20), ylim = c(40, 80))  # Zoom to Northeast Atlantic

# plastic capped to 10%
p_dens_proj <- raster::projectRaster(plastics, crs = proj)

# convert raster to dataframe for ggplot
p_df <- raster::rasterToPoints(p_dens_proj) %>% as_tibble()

sp3 <- ggplot() +
  theme_minimal() +
  scale_fill_viridis(option = "mako", direction = -1, 
                     name = "log (plastic count / km²)",
                     labels = scales::comma) + 
  theme(legend.position = c(0.5, 0),  # Shift the legend vertically
        legend.justification = c(0.5, -0.45),  # Justify the legend to the bottom and center
        legend.box.margin = margin(10, 0, 0, 0),  # Adjust margin around the legend
        plot.title = element_text(hjust = 0.5, size = 17),  # Increase title font size and shift downwards
        legend.key.height = unit(0.2, "in"),
        legend.key.width = unit(2, "in"),  # Make the legend narrower
        legend.direction = "horizontal",
        legend.title = element_text(vjust = 0.8)) +  # Shift the legend title vertically
  geom_raster(aes(x = x, y = y, fill = X00_PlasticsRaster),  
              data = p_df) +
  geom_sf(aes(), colour = NA, fill = "grey75", 
          data = world_prj) +
  geom_sf(aes(), fill = "white", color = "black", 
          data = cookie, size = 0.5) +
  geom_sf(data = colonies_robinson, color = "red", size = 2) +  
  xlab("") + ylab("") +
  ggtitle("Global distribution of marine floating plastic debris")

print(sp3)

# 2. Overlaying selected PAME shapefiles on top of bird distributions

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
PAME_shapefile <- st_read(paste0(datadir,"/PAME/modified_LME.shp"))
plot(PAME_shapefile$geometry)
robinson_PAME <- st_transform(PAME_shapefile, crs = proj)
plot(robinson_PAME$geometry)






















