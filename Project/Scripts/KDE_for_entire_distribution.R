# Image 1: Kernel density estimation (KDE) for the entire distribution----

# Loading essential packages----
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(ggplot2)
library(spData)
library(rnaturalearth)
library(viridisLite)
library(viridis)

graphics.off()  # Close all graphics devices
dev.new()       # Open a new graphics device

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
to_combine <- c("Skjalfandi", "Langanes")
df_mod$colony[df_mod$colony %in% to_combine] <- "Combined"

# Adding a column for species----
df_pooled <- df_mod %>% mutate(species = "Northern fulmar")

# Defining land----
land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# Kernel density estimation according to species----

sub <- df_pooled
  
  for(k in 1:length(unique(sub$tracking_year))){ # First for loop begins
    tracks_wgs <- sub[sub$tracking_year == unique(sub$tracking_year)[k],]

    if(nrow(tracks_wgs) > 50){ # First if loop begins
      
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
      
      kudl <- adehabitatHR::kernelUD(tracks[,"species"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      
      image(kudl)
  
      # Changing from estUDm class to spatial pixels data frame
      kern100 <- adehabitatHR::estUDm2spixdf(kudl)
      plot(kern100)
      kern100_stk <- raster::stack(kern100)
      plot(kern100_stk)
      rast <- kern100_stk
      rast[is.na(rast)] <- 0
      rast <- rast/sum(raster::getValues(rast))
      plot(rast)
      
      wgs84_rast <- raster::projectRaster(rast, crs = proj_wgs84, over = T)
      
      plot(wgs84_rast)
      plot(land, add = T)
        
      # Save raster----
      setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/pooled/")
      raster::writeRaster(wgs84_rast, filename = paste0("r1_",unique(sub$tracking_year)[k],".tif"),
                            format = "GTiff", overwrite = T)
        
      # Save image----
      png(filename = paste0("p1_",(unique(sub$tracking_year))[k],".png"))
      plot(wgs84_rast)
      plot(land, add = T, col = "#66000000")
      dev.off()
        
      
    } # First if loop ends
    
  } #  First for loop ends

# List all raster files starting with "r1_" and ending with ".tif" in the folder
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/pooled/")
raster_files <- list.files(pattern = "^r1_.*\\.tif$")

# Read all rasters into a list
raster_list <- lapply(raster_files, raster)
raster_list

# Defining worldmap----
worldmap <- ggplot2::map_data('world')

# Defining water polygons
water_polygons <- rnaturalearth::ne_download(
  category = "physical",
  type = "ocean",
  scale = 110,
  returnclass = "sf"
)

# Now doing the final plot----
# Define extent_of_desired_extent
extent_of_desired_extent <- extent(c(-90, 90, 30, 90))
# Crop land polygon to desired extent
cropped_land <- crop(land, extent_of_desired_extent)
plot(cropped_land)

# Loop through each raster
for(i in 1:length(raster_list)){
  
  # Set the extent of the raster to match the cropped land
  extent(raster_list[[i]]) <- extent(cropped_land)
  
  # Plot raster
  plot(raster_list[[i]])
  
  # Plot cropped land
  plot(cropped_land, add = TRUE)
  
}



for(i in 1:length(raster_list)){
  
  plot(raster_list[[i]])
  plot(cropped_land, add = T)
  
}


# Define subplot dimensions (adjust width and height as needed)
par(mfrow = c(3, 5), xpd = TRUE, mar = c(0.5, 0.5, 0.5, 0.5), oma = c(0.5, 0.5, 0.5, 0.5), mai = c(0.1, 0.1, 0.3, 0.4))  # Decrease size of subplots

# Define the years for the titles
years <- 2006:2020

# Loop through each raster
for (i in 1:length(raster_list)) {
  
  # Set the extent of the raster to match the cropped land
  extent(raster_list[[i]]) <- extent(cropped_land)
  
  # Plot raster (consider alternative legend formatting)
  plot(raster_list[[i]], main = years[i], axes = FALSE, col.regions = function(x) colorRampPalette(c("blue", "yellow", "red"))(raster_list[[i]]@data$values),  # Color ramp example
       # Option 2 (using format function):
       scales = list(x = ~ format(., scientific = TRUE)))  # Apply scientific notation formatting for axis labels
  
  # Option 3 (external legend with scales):
  # plot the raster without the legend
  # ...
  # Create a legend using scales::format_scientific() for formatting
  # ...
  

# Plot cropped land
plot(cropped_land, add = TRUE, lwd = 0.5)
}


# Display the final plot
print(final_plot)
ggsave("Tracking_year_wise_desnity_distribution", bg = "white", height = 8, width = 14, unit = "in", dpi = 900, plot = p1)