# New script for organized image 1: pooled distribution of NF points with 50% UD contours overlaid

# Extracting 50% UD contours----

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




# Plots----

p1 <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  geom_sf(data = ind_merge_sf, aes(), size = 0.0000000001) +
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

p1_with_colonies <- p1 
p1_with_colonies

# Adding 50% contours on top of this

final_plot <- p1_with_colonies + geom_sf(data = all_boundary_data_sf, aes(color = id), fill = NA, size = 3, lwd = 1) + 
  geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat, color = Colony), color = "white", shape = 17, size = 5) +
  geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat, color = Colony), shape = 17, size = 4) +
  scale_colour_manual(values = rev(colony_colors$Color)) +
  coord_sf(xlim = c(-90, 90), ylim = c(30, 90)) + labs(color = "Colony")
final_plot
ggsave("Northern_fulmar_pooled_distribution.png", height = 8, width= 12, unit = "in", dpi = 900, plot = final_plot, bg = "white")
