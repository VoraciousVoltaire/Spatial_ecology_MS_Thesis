# figuring out the deal with Jan Mayen 2020 and Alkefjellet 

# Loading essential packages
library(adehabitatHR)
library(tidyverse)
library(geosphere)
library(spData)
library(sf)
library(sp)
library(raster)
library(stringr)
library(terra)

# defining output directory----
datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/"

# Defining land----
land <- as(world, "Spatial")

# Defining wgs84 projection----
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# Defining df_mod----
# Creating an sf object----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"
nbs_mylocs_sf <- nbs_mylocs %>% 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) 

# Assorting by months----

all_data <- nbs_mylocs
df <- all_data[!is.na(all_data$timestamp), ]
df <- df %>% dplyr::mutate(month = month(timestamp))
months <- sort(unique(df$month))

# Adding a column for nbs tracking year----

df_mod <- df %>% mutate(tracking_year = ifelse(month < 7, year - 1, year))

# Condensed datasets----
df_mod_5 <- df_mod[df_mod$colony == "Alkefjellet" | (df_mod$colony == "Jan Mayen" & df_mod$tracking_year == 2020),]
df_mod_6 <- df_mod[(df_mod$colony == "Alkefjellet" & (df_mod$tracking_year == 2019 | df_mod$tracking_year == 2020)) | (df_mod$colony == "Jan Mayen" & df_mod$tracking_year == 2020),]
df_mod_7 <- df_mod[df_mod$colony == "Jan Mayen" & df_mod$tracking_year == 2020,]

for(i in unique(df_mod_7$colony)){ # First for loop begins
  sub <- df_mod_7 %>% filter(colony == i)
  
  
  for(j in 1:length(unique(sub$tracking_year))){ # Second for loop begins
    tracks_nbs <- sub[sub$tracking_year == unique(sub$tracking_year)[j],]
    
    # Setting a null grid for kernel estimation 
    if(nrow(tracks_nbs) > 4){ # First if loop starts
      
      # Setting a colony-centered crs 
      median_loc <- cbind(median(tracks_nbs$lon), median(tracks_nbs$lat))
      DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",median_loc[1]," +lat_0=",median_loc[2]," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m ", sep="")) 
      
      if(min(tracks_nbs$lon) <= -179 ){ lon_min <- -180
      } else {lon_min <- floor(min(tracks_nbs$lon))-1 }
      
      if(max(tracks_nbs$lon) >= 179){ lon_max <- 180
      } else { lon_max <- ceiling(max(tracks_nbs$lon))+1 }
      
      if(min(tracks_nbs$lat) <= -89 ){ lat_min <- -90 
      } else { lat_min <- floor(min(tracks_nbs$lat))-1 }
      
      if(max(tracks_nbs$lat) >= 89){ lat_max <- 90
      } else { lat_max <- ceiling(max(tracks_nbs$lat))+1 }
      
      so.grid <- expand.grid(LON = seq(lon_min, lon_max, by=1), 
                             LAT = seq(lat_min, lat_max, by=1))
      
      sp::coordinates(so.grid) <- ~LON+LAT
      sp::proj4string(so.grid) <- sp::proj4string(land)
      
      so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj) # transforming null grid to new crs
      coords <- so.grid.proj@coords
      
      c <- min(coords[,1])-1000000   ## to check my min lon
      d <- max(coords[,1])+1000000   ## to check my max lon
      
      e <- min(coords[,2])-1000000   ## to check my min lat
      f <- max(coords[,2])+1000000   ## to check my max lat
      
      a <- seq(c, d, by=10000)
      b <- seq(e, f, by=10000)
      null.grid <- expand.grid(x = a,y = b)
      sp::coordinates(null.grid) <- ~ x + y
      sp::gridded(null.grid) <- TRUE
      
      # Creating a colony-specific spatial points data frame
      sp::coordinates(tracks_nbs) <- ~ lon + lat
      sp::proj4string(tracks_nbs) <- sp::proj4string(land)
      
      tracks <- sp::spTransform(tracks_nbs, CRS = DgProj)
      tracks$individ_id <- factor(tracks$individ_id)
      
      # Kernel estimation
      kudl <- kernelUD(tracks[,"tracking_year"], grid = null.grid, h = 200000) ## smoothing factor equals 200 km for GLS data
      
      
      # Now, calculating pers according to tracking year----
      
      
      vud <- adehabitatHR::getvolumeUD(kudl)
      
      fud <- vud[[1]]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- data.frame(hr95)
      coordinates(hr95) <- coordinates(fud)
      sp::gridded(hr95) <- TRUE
      
      kde_spixdf <- adehabitatHR::estUDm2spixdf(kudl)
      kern95 <- kde_spixdf
      
      stk_100 <- raster::stack(kern95)
      stk_95 <- raster::stack(hr95)
      
      sum_all_100 <- stk_100[[1]]
      sum_all_95 <- stk_95[[1]]
      
      sum_all_raw <- sum_all_100 * sum_all_95
      
      rast <- sum_all_raw/sum(raster::getValues(sum_all_raw))
      rast[rast == 0] <- NA
      
      x.matrix <- is.na(as.matrix(rast))
      colNotNA <- which(colSums(x.matrix) != nrow(rast))
      rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
      
      croppedExtent <- raster::extent(rast, 
                                      r1 = rowNotNA[1]-2, 
                                      r2 = rowNotNA[length(rowNotNA)]+2,
                                      c1 = colNotNA[1]-2, 
                                      c2 = colNotNA[length(colNotNA)]+2)
      
      cropped <- raster::crop(rast, croppedExtent)
      cropped[is.na(cropped)] <- 0
      
      # Changing land's projection for masking
      mask_proj <- sp::spTransform(land, DgProj)
      mask_proj_pol <- as(mask_proj, "SpatialPolygons") 
      
      rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = TRUE)
      rast_mask <- rast_mask_na
      rast_mask[is.na(rast_mask)] <- 0
      rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
      
      rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE)
      rast_mask_final2 <- rast_mask_final 
      
      mask_wgs84 <- projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
      
  
      
    
      
      
      
      KDE_ref <- paste0(i, "_",unique(sub$tracking_year)[j])
      raster::writeRaster(rast_mask_final, filename = paste0(datadir,"input_rasters/tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      raster::writeRaster(mask_wgs84, filename = paste0(datadir,"output/unique_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      png(filename = paste0(datadir,"output/unique_distributions/",KDE_ref,".png"))
      plot(mask_wgs84, main = paste0("",i))
      plot(land, add = T, col = "#66000000")
      dev.off()
      
     
      
      
    } # First if loop ends
  } # Second for loop ends
} # First for loop ends

# No clue where the error lies now