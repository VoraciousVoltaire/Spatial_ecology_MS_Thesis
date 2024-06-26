# Removing irma and trying to recreate the result. Just removing rows that have IRMA as their entry in loc_type and recalculating kernels on 
# a month-wise basis, multiplying them with the plastic raster allowed by calculating a median of pers for those 6 nbs months and arriving at 
# a single nbs pers for each colony. 

# This is a (relatively) clean script which calculates kernels on a month-wise basis, then sums months
# into one non-breeding season raster for each colony.

# Clearing environment----
remove(list = ls())

# Loading relevant packages----

library(dplyr)
library(lubridate)
library(here)
library(ggplot2)
library(utils)
library(sf)
library(sp)
library(raster)
library(terra)
library(dplyr)
library(tidyverse)
library(adehabitatHR) 
library(spData) # for world vector
library(cowplot)
library(stringr)
library(RColorBrewer)
library(viridis)
library(viridisLite)


# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/"
mylocs <- readRDS(paste0(datadir, "SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"summaryTable.rds"))

# Exploring data----
mylocs %>%
  dplyr::left_join(summary_info, join_by(ring)) %>%
  dplyr::mutate(year = year(timestamp)) %>%
  dplyr::group_by(colony,year) %>%
  dplyr::summarise(n_birds_tracked = n_distinct(ring))

# Creating an sf object----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 

# Just needed to add this line----
indiv_merged_df <- indiv_merged_df[indiv_merged_df$loc_type != "IRMA",]
unique(indiv_merged_df$loc_type)

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

# Setting a directory for month-wise output rasters
dir_kernels <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/GLS_outputs"

# Defining land----

land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# In case, one wants to exactly replicate what the petrels' paper did
# setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/")
# land <- sf::read_sf(dsn = "baselayer", layer = "ne_10m_land")
# plot(land)
# proj_wgs84 <- crs(land)

# For loop starts here----

for(i in unique(df$colony)){ # First for loop begins
  sub <- df %>% filter(colony == i)
  
  for(month_number in 1:length(months)){ # Second for loop begins
    tracks_wgs <- sub[sub$month == months[month_number],] 
    
    sub <- subset(sub, lat < 90)  
    png(filename = paste0(dir_kernels, "/distributions/", i, ".png"))  
    plot(lat~lon, data = sub, type= "n", asp = 1, 
         xlim=c(min(all_data$lon) + 0.2,max(all_data$lon) - 0.2), 
         ylim=c(min(all_data$lat) + 0.5, max(all_data$lat) - 0.5), 
         main="", frame = T, xlab="", ylab="")
    plot(land, col='lightgrey', add= T)
    points(lat~lon, data = sub, pch = 16, cex = 0.5, col="blue")
    points(lat~lon, data = tracks_wgs, pch = 16, cex = 0.5, col="red")
    dev.off()
    
    # Kernel density estimation----
    
    if(nrow(tracks_wgs) > 4){
      
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
      tracks$month <- factor(tracks@data$month)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"month"], 
                                     grid = null.grid, h = 200000)  # smoothing factor equals 200 km for GLS data
      
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
      
      KDE_ref <- paste0(i, "_", months[month_number])
      raster::writeRaster(rast_mask_final, filename = paste0(datadir,"GLS_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      raster::writeRaster(mask_wgs84, filename = paste0(dir_kernels,"/unique_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      png(filename = paste0(dir_kernels,"/unique_distributions/",KDE_ref,".png"))
      plot(mask_wgs84, main = paste0("",i))
      plot(land, add = T, col = "#66000000")
      dev.off()
      
      View(df[df$colony == "Little Saltee",])
      
    } # First if loop ends
  } # Second for loop ends
} # First for loop ends

# Multiplication script starts here

# Loading data----
dir_demClasses <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/GLS_tifs/" # input directory
dir_1by1 <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/GLS_outputs/multiplication_rasters" # output directory
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data")
plastics <- raster::raster("00_PlasticsRaster.tif")

r <- raster()
r

yelblus <- c(brewer.pal(n = 9, name = "YlGnBu"),"#00172e")
cols <- colorRampPalette(yelblus)(255)
colsviri <- cols[20:255]

colsinf <- rev(inferno(200))

plastics2 <- plastics
plastics2[is.na(plastics2)] <- 0 
p_sum1    <- plastics2/sum(raster::getValues(plastics2))
p_sum1[is.na(plastics)] <- NA

RES <- res(plastics) 
R <- 6371007.2 
lat <- raster::yFromRow(plastics, 1:nrow(plastics)) # latitude of the centroid of each cell (in degrees, need to be converted in radians)
area <- (sin(pi/180*(lat + RES[2]/2)) - sin(pi/180*(lat - RES[2]/2))) * (RES[1] * pi/180) * R^2
r_area <- raster::setValues(plastics, rep(area, each = ncol(plastics))) # gives the area of each grid cell in meters 
plot(r_area, col = colsviri)

files <- list.files(dir_demClasses, full.names = TRUE, pattern=".*\\.tif$"); files

nas <- as.data.frame(files)
nas$name <- NA
nas$vals <- NA
nas$nas <- NA
nas$files <- NULL

dat <- data.frame()
result <- c()

for (i in 1:length(files)){ # First for loop begins
  
  a <- raster(files[i])
  name <- a@data@names[1]
  nas$name[i] <- name
  a[is.na(a)] <- 0 
  
  b <- sum(raster::getValues(a)) 
  
  a_proj <- raster::projectRaster(a, plastics, method = "bilinear")
  print(a_proj)
  
  a_proj2 <- a_proj * r_area / 100000000 # rescaling the values in each cell
  a_proj2[is.na(a_proj2)] <- 0 
  
  c <- sum(values(a_proj2))
  
  na_a_proj2 <- a_proj2
  na_a_proj2[na_a_proj2 > 0] <- 1
  na_a_proj2[is.na(na_a_proj2)] <- 0
  nas$vals[i] <- sum(na_a_proj2@data@values)
  
  na_p_sum1 <- p_sum1
  na_p_sum1[is.na(na_p_sum1)] <- 1
  na_p_sum1[na_p_sum1 != 1] <- 0
  na_over <- na_a_proj2 * na_p_sum1
  nas$nas[i] <- sum(na_over@data@values) 
  
  a_proj2[is.na(plastics)] <- NA
  
  raster_name_1 <- gsub(dir_demClasses, "", files[i])
  print(raster_name_1) 
  
  raster_name_2 <- paste0(dir_1by1, raster_name_1)
  
  raster::writeRaster(a_proj2, filename=raster_name_2, format="GTiff", overwrite=TRUE)
  
  over <- a_proj2 * p_sum1
  
  over_score <- over
  summary(over_score@data@values)
  over_score[is.na(over_score)] <- 0
  exposure_score <- round(sum(raster::getValues(over_score))*1000000,4)
  
  png(paste0(dir_1by1,"/maps/",name,".png"), width=1399,height=455)
  par(mfrow=c(1,2))
  plot(a_proj2,main=paste0(name," distribution"),col=colsviri,legend=F)
  plot(over,main=paste0("Exposure score = ", exposure_score),
       col=colsinf,legend=F)
  
  dev.off()
  
  if(sum(na_over@data@values)> 0){ # First if loop starts
    png(paste0(dir_1by1,"/na_maps/",name,".png"), width=1399,height=455)
    par(mfrow=c(1,2))
    plot(a_proj2,main=name,col=colsviri,legend=F)
    plot(na_over,main=paste0("Exposure Nas = ",sum(na_over@data@values)),
         col=colsinf,legend=F)
    dev.off()
  } # First if loop ends
  
  name_split <- strsplit(name,"_")[[1]]
  month <- name_split[length(name_split)]
  population <- paste(name_split[1])
  
  check <- cbind(b,c)
  result <- cbind(population, month, check, exposure_score)
  dat <- rbind(dat, as.data.frame(result))
  print(i)
} # First for loop ends

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/GLS_outputs/csv/")
write.csv(dat, "GLS_exposure_scores_by_month.csv",
          row.names = F)  
nas$percent_na <- nas$nas/nas$vals*100
head(nas)
write.csv(nas, "GLS_NAs_table.csv",
          row.names = F) 

#2 don't exist ----------this isn't required in our case (except for Alkefjellet which has a very high variance)
nas_no_na <- subset(nas,vals > 0)
nas_no_na
mean(nas_no_na$percent_na)

exposure_score_csv <- read.csv("GLS_exposure_scores_by_month.csv")
pop_exposure <- exposure_score_csv %>%
  group_by(population) %>%  
  summarise(population_exposure = round(mean(exposure_score), 4)) %>%
  data.frame() 
write.csv(pop_exposure, "GLS_exposure_scores_by_population.csv",
          row.names = F) 
Species_exposure_score <- mean(pop_exposure$population_exposure)
Species_exposure_score 

