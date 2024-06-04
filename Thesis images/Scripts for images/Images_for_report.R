# Images for report

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
library(adehabitatHR)
library(dplyr)

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
ind_merge_sf <- df_mod_sf

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
m

# Mimicking Clark et al. (2023) script----

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data/") # directory where plastics raster is stored

land <- as(world, "Spatial")

# # Distributions of all birds 
# bird_dist <- raster::raster("08_all_species_distribution.tif")
# b <- bird_dist
# b[is.na(b)] <- 0 
# b_sum1 <- b/sum(raster::getValues(b))
# b_sum1[is.na(bird_dist)] <- NA
# b_sum1[b_sum1 == 0] <- NA

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
  theme(legend.position = c(0.5, 0.02),  # Shift the legend vertically
        legend.justification = c(0.5, -0.45),  # Justify the legend to the bottom and center
        legend.box.margin = margin(10, 0, 0, 0),  # Adjust margin around the legend
        plot.title = element_text(hjust = 0.5, size = 17, element_text(face = "bold")),  # Increase title font size and shift downwards
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

sp3 <- ggplot() +
  theme_minimal() +
  scale_fill_viridis(option = "inferno", direction = -1, 
                     name = "log (plastic count / km²)",
                     labels = scales::comma) + 
  theme(legend.position = c(0.5, 0.02),  # Shift the legend vertically
        legend.justification = c(0.5, -0.45),  # Justify the legend to the bottom and center
        legend.box.margin = margin(10, 0, 0, 0),  # Adjust margin around the legend
        plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),  # Make the ggtitle bold
        legend.key.height = unit(0.2, "in"),
        legend.key.width = unit(2, "in"),  # Make the legend narrower
        legend.direction = "horizontal",
        legend.title = element_text(vjust = 0.8, face = "bold")) +  # Shift the legend title vertically
  geom_raster(aes(x = x, y = y, fill = X00_PlasticsRaster),  
              data = p_df) +
  geom_sf(aes(), colour = NA, fill = "grey75", 
          data = world_prj) +
  geom_sf(aes(), fill = "white", color = "black", 
          data = cookie, size = 0.5) +  
  xlab("") + ylab("") +
  ggtitle("Global distribution of marine floating plastic debris")


print(sp3)
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/")
ggsave("Inferno_Global_distribution_of_marine_floating_plastic_debris.png", height = 8, width = 14, unit = "in", dpi = 900, plot = sp3, bg = "white")

# 2. Overlaying selected PAME shapefiles on top of bird distributions

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
PAME_shapefile <- st_read(paste0(datadir,"/PAME/modified_LME.shp"))
plot(PAME_shapefile$geometry)
robinson_PAME <- st_transform(PAME_shapefile, crs = proj)
plot(robinson_PAME$geometry)
PAME_curated <- PAME_shapefile[c(2,8,9,12,14),]
PAME_curated

# 3. Pooled northern fulmar distribution and 

worldmap <- ggplot2::map_data('world')
p1 <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group, fill = "lightgrey")) +
  geom_sf(data = ind_merge_sf, aes(), size = 0.1) +
  coord_sf(xlim = c(-90, 90), ylim = c(30, 90)) +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold")  # Set bold face for axis titles
  ) +
  ggtitle("Northern fulmar pooled distribution") +
  xlab("Longitude") +
  ylab("Latitude"); p1

# Load water polygons
water_polygons <- rnaturalearth::ne_download(
  category = "physical",
  type = "ocean",
  scale = 110,
  returnclass = "sf"
)

# Create the plot

p1 <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  geom_sf(data = ind_merge_sf, aes(), size = 0.1) +
  coord_sf(xlim = c(-90, 90), ylim = c(30, 90)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold")  # Set bold face for axis titles
  ) +
  ggtitle("Northern fulmar pooled distribution") +
  xlab("Longitude") +
  ylab("Latitude")

# Display the plot
print(p1) # this is OK; now just have to add contours of 95% UD home ranges in different colours- and the colour of the colony in the legend has to be the same as the contour

# Plotting colonies in different colours on top of this

p1_with_colonies <- p1 + geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat, color = colony))
p1_with_colonies

# Density gradient----

library(ggplot2)
library(MASS)  # For kde2d function

# Extract X and Y coordinates from the geometry column
ind_merge_coords <- as.data.frame(st_coordinates(ind_merge_sf))

# Compute density values using kde2d function
dens <- kde2d(ind_merge_coords$X, ind_merge_coords$Y)

# Convert density values to a data frame format
dens_df <- expand.grid(x = dens$x, y = dens$y)
dens_df$z <- c(dens$z)

# Density plot with adjusted transparency for cells with zero density
grad_p <- ggplot(dens_df, aes(x = x, y = y, fill = z)) +
  geom_raster(aes(alpha = ifelse(z == 0, 0, 1))) +  # Set transparency for zero density cells
  scale_fill_viridis(option = "viridis") +  # Viridis color scale
  theme_void()  # Remove background and axis

# Overlaying density plot on top of p1
overlay_plot <- p1 + 
  annotation_custom(ggplotGrob(grad_p), xmin = -90, xmax = 90, ymin = 30, ymax = 90)

# Displaying the overlay plot
print(overlay_plot)

##################
worldmap <- ggplot2::map_data('world')
PAME_plot <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold")  # Set bold face for axis titles
  ) +
  ggtitle("Northern fulmar pooled distribution") +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_sf(data = PAME_shapefile, aes(), fill = "NA") +
  geom_sf(data = PAME_curated, aes(), fill = "white") +
  geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat, color = colony)) +
  coord_sf(xlim = c(-95, 95), ylim = c(30, 90))
PAME_plot

# Your existing code...

PAME_plot <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold"),  # Set bold face for axis titles
    legend.title = element_text(face = "bold"),  # Set bold face for legend title
    legend.text = element_text(face = "bold")  # Set bold face for legend text
  ) +
  ggtitle("Relevant Protection of the Arctic Marine Environment (PAME) zones") +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_sf(data = PAME_shapefile, aes(), fill = NA, color = "black", size = 7) +
  geom_sf(data = PAME_curated, aes(), fill = "white", color = "black", size = 7) +
  geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat, color = colony)) +
  geom_sf_label(data = PAME_curated, aes(label = LME_NAME), 
                size = 3.5, color = "black", fill = "white", label.padding = unit(0.1, "lines"),
                label.size = 0.5) +  # Adjust size and aesthetics of the label
  coord_sf(xlim = c(-95, 95), ylim = c(30, 90)) +  # Set limits of the plot
  labs(color = "Colony")  # Set legend title
print(PAME_plot)
ggsave("Relevant_PAME_zones.png", height = 8, width= 12, unit = "in", dpi = 900, plot = PAME_plot, bg = "white")

# Now overlaying this with all year round 95% KDEs

output_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/vectors/"
colonies_to_combine <- c("Skjalfandi", "Langanes")
df_mod$colony[df_mod$colony %in% colonies_to_combine] <- "Iceland"

df_mod <- df_mod %>% filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp))
for(i in unique(df_mod$colony)){ # First for loop begins
  sub <- df_mod %>% filter(colony == i)
  tracks_wgs <- subset(sub, lat < 90)  
  
  
  # Kernel density estimation----
  
  if(nrow(tracks_wgs) > 50){
    
    if(min(tracks_wgs$lon) <= -179 ){ lon_min <- -180
    } else {lon_min <- floor(min(tracks_wgs$lon))-1 }
    
    if(max(tracks_wgs$lon) >= 179){ lon_max <- 180
    } else { lon_max <- ceiling(max(tracks_wgs$lon))+1 }
    
    if(min(tracks_wgs$lat) <= -89 ){ lat_min <- -90 
    } else { lat_min <- floor(min(tracks_wgs$lat))-1 }
    
    if(max(tracks_wgs$lat) >= 89){ lat_max <- 90
    } else { lat_max <- ceiling(max(tracks_wgs$lat))+1 }
    
    so.grid <- expand.grid( LON = seq(lon_min, lon_max, by=1), 
                            LAT = seq(lat_min, lat_max, by=1))
    
    sp::coordinates(so.grid) <- ~LON+LAT
    crs(so.grid) <- proj_wgs84
    
    # Setting a colony-centric crs
    mean_loc <- geosphere::geomean(cbind(tracks_wgs$lon,tracks_wgs$lat))
    DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                             " +lat_0=",mean_loc[2])) 
    
    so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj)
    coords <- so.grid.proj@coords
    
    c <- min(coords[,1])-1000000   ## to check my min lon
    d <- max(coords[,1])+1000000   ## to check my max lon
    
    e <- min(coords[,2])-1000000   ## to check my min lat
    f <- max(coords[,2])+1000000   ## to check my max lat
    
    a <- seq(c, d, by=10000)
    b <- seq(e, f, by=10000)
    null.grid <- expand.grid(x = a,y = b)
    sp::coordinates(null.grid) <- ~x+y
    sp::gridded(null.grid) <- TRUE
    
    # Converting tracks_wgs into a spatial points data frame
    sp::coordinates(tracks_wgs) <- ~lon+lat
    crs(tracks_wgs) <- proj_wgs84
    
    tracks <- sp::spTransform(tracks_wgs, CRS = DgProj)
    
    kudl <- adehabitatHR::kernelUD(tracks[,"colony"], 
                                   grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
    hr_vertices <- adehabitatHR::getverticeshr(kudl, percent = 50) # changed from 95 to 50
    
    # Convert SpatialPolygonsDataFrame to sf object
    hr_sf <- sf::st_as_sf(hr_vertices)
    hr_sf <- st_transform(hr_sf, crs = 4326)
    
    # Save the shapefile using sf package
    sf::st_write(hr_sf, paste0(output_dir, "/colony_", i, "_home_range.shp"), append = F)
    
  } # first if loop ends
} # first for loop ends

# List all shapefiles in the output directory
shapefiles <- list.files(output_dir, pattern = ".shp", full.names = TRUE)

# Create an empty data frame to store colony names and corresponding shapefiles
colony_data <- data.frame()

for (shapefile in shapefiles) {
  # Read the shapefile
  data <- st_read(shapefile)
  
  # Extract boundaries
  boundaries <- st_boundary(data)
  
  # Define the output filename
  output_filename <- paste0(output_dir, gsub(".shp", "_boundary.shp", basename(shapefile)))
  
  # Write the boundaries to a new shapefile
  st_write(boundaries, output_filename)
}
boundary_files <- list.files(output_dir, pattern = "_home_range_boundary\\.shp$", full.names = TRUE)
boundary_files

# Read shapefiles into sf objects
boundary_files <- lapply(boundary_files, st_read)

# Combine all sf objects into one
all_boundary_data_sf <- do.call(rbind, boundary_files)

 # Nous sommes ici
PAME_final_plot <- PAME_plot + geom_sf(data = all_boundary_data_sf, aes(color = id, geometry = geometry), fill = NA, size = 0.2) + coord_sf(xlim = c(-95, 95), ylim = c(30, 90)) +labs(color = "Colony") + 
  geom_sf_label(data = PAME_curated, aes(label = LME_NAME), 
                size = 3.5, color = "black", fill = "white", label.padding = unit(0.1, "lines"),
                label.size = 0.5)   # Adjust size and aesthetics of the label
PAME_final_plot
ggsave("Relevant_PAME_zones_with_contours.png", height = 8, width= 12, unit = "in", dpi = 900, plot = PAME_final_plot, bg = "white")

# library(ggplot2)
# library(rnaturalearth)
# 
# # Load land polygons from rnaturalearth
# land_polygons <- ne_countries(scale = "medium", returnclass = "sf")
# 
# # Create the plot
# p1 <- ggplot() +
#   geom_sf(data = land_polygons, fill = "lightgrey") +
#   geom_sf(data = water_polygons, fill = "lightblue", color = "black") +
#   geom_sf(data = ind_merge_sf, aes(), size = 0.1) + 
#   coord_sf(xlim = c(-90, 90), ylim = c(30, 90)) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
#     axis.title = element_text(face = "bold")  # Set bold face for axis titles
#   ) +
#   ggtitle("Northern fulmar pooled distribution") +
#   xlab("Longitude") +
#   ylab("Latitude")
# 
# # Display the plot
# print(p1)

# Copy-pasting 95% UD script----

# Defining output directory
output_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/vectors/"

# Defining land

land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

colonies_to_combine <- c("Skjalfandi", "Langanes")
df_mod$colony[df_mod$colony %in% colonies_to_combine] <- "Iceland"
df_mod_nbs <- df_mod %>%
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) 

for(i in unique(df_mod_nbs$colony)){ # First for loop begins
  sub <- df_mod_nbs %>% filter(colony == i)
  tracks_wgs <- subset(sub, lat < 90)  
  
    
    # Kernel density estimation----
    
    if(nrow(tracks_wgs) > 50){
      
      if(min(tracks_wgs$lon) <= -179 ){ lon_min <- -180
      } else {lon_min <- floor(min(tracks_wgs$lon))-1 }
      
      if(max(tracks_wgs$lon) >= 179){ lon_max <- 180
      } else { lon_max <- ceiling(max(tracks_wgs$lon))+1 }
      
      if(min(tracks_wgs$lat) <= -89 ){ lat_min <- -90 
      } else { lat_min <- floor(min(tracks_wgs$lat))-1 }
      
      if(max(tracks_wgs$lat) >= 89){ lat_max <- 90
      } else { lat_max <- ceiling(max(tracks_wgs$lat))+1 }
      
      so.grid <- expand.grid( LON = seq(lon_min, lon_max, by=1), 
                              LAT = seq(lat_min, lat_max, by=1))
      
      sp::coordinates(so.grid) <- ~LON+LAT
      crs(so.grid) <- proj_wgs84
      
      # Setting a colony-centric crs
      mean_loc <- geosphere::geomean(cbind(tracks_wgs$lon,tracks_wgs$lat))
      DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                               " +lat_0=",mean_loc[2])) 
      
      so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj)
      coords <- so.grid.proj@coords
      
      c <- min(coords[,1])-1000000   ## to check my min lon
      d <- max(coords[,1])+1000000   ## to check my max lon
      
      e <- min(coords[,2])-1000000   ## to check my min lat
      f <- max(coords[,2])+1000000   ## to check my max lat
      
      a <- seq(c, d, by=10000)
      b <- seq(e, f, by=10000)
      null.grid <- expand.grid(x = a,y = b)
      sp::coordinates(null.grid) <- ~x+y
      sp::gridded(null.grid) <- TRUE
      
      # Converting tracks_wgs into a spatial points data frame
      sp::coordinates(tracks_wgs) <- ~lon+lat
      crs(tracks_wgs) <- proj_wgs84
      
      tracks <- sp::spTransform(tracks_wgs, CRS = DgProj)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"colony"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      hr_vertices <- adehabitatHR::getverticeshr(kudl, percent = 95)
      
      # Convert SpatialPolygonsDataFrame to sf object
      hr_sf <- sf::st_as_sf(hr_vertices)
      hr_sf <- st_transform(hr_sf, crs = 4326)
      
      # Save the shapefile using sf package
      sf::st_write(hr_sf, paste0(output_dir, "/colony_", i, "_home_range.shp"), append = F)
      
} # first if loop ends
} # first for loop ends

# List all shapefiles in the output directory
shapefiles <- list.files(output_dir, pattern = ".shp", full.names = TRUE)

# Create an empty data frame to store colony names and corresponding shapefiles
colony_data <- data.frame()

# Loop through each shapefile
for (file in shapefiles) {
  # Read the shapefile
  shapefile_data <- sf::st_read(file)
  
  # Extract colony name from the file name
  colony_name <- gsub("_home_range.shp", "", gsub(".*colony_", "", file))
  
  # Add colony name and shapefile path to the data frame
  colony_data <- rbind(colony_data, data.frame(colony = colony_name, shapefile = file))
}

# Now, loop through each colony and add it to the plot
for (i in 1:nrow(colony_data)) {
  # Read the shapefile
  colony_shapefile <- sf::st_read(colony_data[i, "shapefile"])
  
  # Extract colony name
  colony_name <- colony_data[i, "colony"]
  
  # Add contour shape to the plot
  p1_with_colonies <- p1_with_colonies + geom_sf(data = colony_shapefile, fill = NA, color = "black", size = 1) +
    geom_sf(data = colony_shapefile, aes(fill = as.factor(colony_name)), show.legend = "polygon")
}
p1_with_colonies

# Update the plot aesthetics and legend
p2 <- p1_with_colonies +
  scale_fill_manual(values = rainbow(nrow(colony_data))) +  # Set colors for colonies
  labs(fill = "Colony") +  # Set legend label
  theme(legend.position = "right")  # Adjust legend position

# Display the plot
print(p2)

# Now overlaying PAME shapefiles on top of this----



# try
# Required packages
library(dplyr)
library(adehabitatHR)
library(sp)

# Original ggplot code
p1 <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  geom_sf(data = ind_merge_sf, aes(), size = 0.1) +
  coord_sf(xlim = c(-90, 90), ylim = c(30, 90)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold")  # Set bold face for axis titles
  ) +
  ggtitle("Northern fulmar pooled distribution") +
  xlab("Longitude") +
  ylab("Latitude")

# Define the output directory
output_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/vectors/"

# Get a list of shapefiles in the output directory
shapefiles <- list.files(output_dir, pattern = "\\.shp$", full.names = TRUE)
for (shapefile in shapefiles) {
  # Read the shapefile
  data <- st_read(shapefile)
  
  # Extract boundaries
  boundaries <- st_boundary(data)
  
  # Define the output filename
  output_filename <- paste0(output_dir, gsub(".shp", "_boundary.shp", basename(shapefile)))
  
  # Write the boundaries to a new shapefile
  st_write(boundaries, output_filename)
}
boundary_files <- list.files(output_dir, pattern = "_home_range_boundary\\.shp$", full.names = TRUE)

# Read shapefiles into sf objects
boundary_files <- lapply(shapefiles, st_read)

# Combine all sf objects into one
all_boundary_data_sf <- do.call(rbind, boundary_files)


# the last straw



# Plotting

final_plot <- p1_with_colonies + geom_sf(data = all_boundary_data_sf, aes(color = id, geometry = geometry), fill = NA, size = 0.2) + coord_sf(xlim = c(-90, 90), ylim = c(30, 90)) + labs(color = "Colony") +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))
ggsave("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/northern_fulmar_pooled_distribution_v8.png", height = 8, width = 12, unit = "in", dpi = 900, plot = final_plot, bg = "white")

final_plot_2 <- final_plot + 
theme(
  legend.background = element_rect(fill = "white", color = "black"),
  plot.title = element_text(hjust = 0.5,face = "bold", color = "black", background = "white"),
  legend.title = element_text(color = "black", background = "white")
) +
  ggtitle("Northern fulmar distribution") +
  theme(plot.title = element_text(vjust = 0.5, size = 14))



plot(all_boundary_data_sf$geometry, add = T)
plot(ind_merge_sf$geometry)

p3 <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  geom_sf(data = ind_merge_sf, aes(), size = 0.1) +
  coord_sf(xlim = c(-90, 90), ylim = c(30, 90)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold")  # Set bold face for axis titles
  ) +
  ggtitle("Northern fulmar pooled distribution") +
  xlab("Longitude") +
  ylab("Latitude") +
  # Plot boundary files with unique border colors
  geom_sf(data = all_boundary_data_sf, aes(color = id, geometry = geometry), size = 0.2) +
  scale_color_manual(values = colony_colors, 
                     name = "Colony",
                     labels = colony_names)

# Show the plot
print(p2)

correct_exposure_scores_by_month <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/correct_exposure_scores_by_month.csv")
var_df <- correct_exposure_scores_by_month %>% group_by(population) %>% summarise(var_PERS = var(exposure_score))
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/")
write.csv(var_df, "month-wise_variance_PERS.csv")

second_hyp_df <- read.csv("final_df_for_hyp_2.csv")
total_n_tracks_PAME <- second_hyp_df %>% group_by(LME_NAME) %>% summarise(unique_inds = length(unique(individ_id)),n = n())
colnames(total_n_tracks_PAME) <- c("PAME zone","Total number of unique individuals", "Total number of relocations")
View(total_n_tracks_PAME)
write.csv(total_n_tracks_PAME, "total_n_tracks_PAME.csv")
stats_second_hyp_df <- second_hyp_df %>% group_by(LME_NAME, colony) %>% summarise(n = n(), unique = length(unique(individ_id)))
View(stats_second_hyp_df)
total <- second_hyp_df %>% group_by(LME_NAME) %>% summarise(total = n())
proportion_stats <- merge(stats_second_hyp_df, total, by = "LME_NAME") 
View(proportion_stats)
percentage_colony <- proportion_stats %>%
  mutate(percentage = sprintf("%.2f", (n / total) * 100))
View(percentage_colony)
write.csv(percentage_colony, "percentage_colony.csv")

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Supplementary_info/csv/")
total_n_tracks <- df_mod %>% group_by(colony) %>% summarise(unique_inds = length(unique(individ_id)), total_n_tracks = n())
colnames(total_n_tracks) <- c("Colony", "Total number of unique individuals","Total number of relocations")
View(total_n_tracks)
write.csv(total_n_tracks, "total_n_relocations.csv")

# Image 4.0----
# Jitter plot----

# Loading data: 

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/latest_right_attempt/outputs/csv")
df <- read.csv("correct_exposure_scores_by_individual_according_to_tracking_yr.csv")
corrected_df <- df[df$exposure_score > 10,]
analysis_df <- corrected_df[!(corrected_df$population == "Jan.Mayen" & 
                                corrected_df$individual == "NOS.4181597" & 
                                corrected_df$tracking_year == 2018), ]
analysis_df <- analysis_df[,-c(4,5)]
colnames(analysis_df) <- c("Colony", "Individ_id", "Tracking_year", "pers")
analysis_df$Colony[analysis_df$Colony == "Combined"] <- "Iceland"
analysis_df$Colony <- gsub("\\.", " ", analysis_df$Colony)
analysis_df$Colony <- as.factor(analysis_df$Colony)
order_to_match <- c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Iceland", "Faroe Islands", "Jarsteinen", "Eynhallow", "Isle of Canna", "Inishkea", "Little Saltee")
reversed_order_to_match <- rev(order_to_match)
analysis_df$Colony <- factor(analysis_df$Colony, levels = reversed_order_to_match)

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/")
color_vector <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF",
                      "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")

jitter_plot <- ggplot(analysis_df, aes(x = pers, y = Colony, col = Colony)) +
  geom_jitter(width = 0, height = 0.2) +  # Add jitter for better visualization
  scale_color_manual(values = color_vector) +  # Specify colors manually
  labs(x = "Plastic Exposure Risk Score (PERS)", y = "Colony") +
  theme_minimal() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  ggtitle("Tracking-year wise individual-specific PERS across Northern fulmar colonies") +
  guides(col = guide_legend(reverse = TRUE))  # Reverse order of legend

jitter_plot

ggsave("Jitter_plot_for_PERS.png", height = 8, width = 14, unit = "in", dpi = 900, plot = jitter_plot, bg = "white")
           