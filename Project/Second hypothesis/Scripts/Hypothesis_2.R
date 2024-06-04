# Hypothesis 2: assorting distribution points according to the region in which they fall on, then calculating separate home ranges for 'em, then 
# multiplying 'em with the plastic raster and getting a separate pers for each region and then finally perform the correlation test 

# Clearing environment----
remove(list = ls())

# Loading essential packages----
library(dplyr)
library(rnaturalearth)
library(sf)
library(dplyr)
library(lubridate)
library(here)
library(ggplot2)
library(utils)
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
library(ppcor)

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
indiv_merged_df$colony[indiv_merged_df$colony %in% colonies_to_combine] <- "Iceland"
ind_merge_sf <- indiv_merged_df %>% sf::st_as_sf(coords = c('lon', 'lat'), crs = 4326)


# Plotting pooled distribution----
worldmap <- ggplot2::map_data('world')

# Loading shapefiles to overlay on the land raster----

PAME_shapefile <- st_read(paste0(datadir,"/PAME/modified_LME.shp"))
plot(PAME_shapefile$geometry)
OSPAR_shapefile <- st_read(paste0(datadir, "/OSPAR/OSPAR_subregions_20160418_3857.shp"))
st_crs(OSPAR_shapefile)

p1 <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group)) +
  geom_sf(data = ind_merge_sf, aes(), size = .1) +
  coord_sf(xlim = c(-90,90), ylim = c(10, 90)) +
  theme_light() +
  ggtitle("Northern fulmar pooled distribution") +
  theme(plot.title = element_text(hjust = 0.5, size = 17))

p1 <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "grey") +
  geom_sf(data = ind_merge_sf, aes(), size = 0.1) +
  coord_sf(xlim = c(-90, 90), ylim = c(30, 90)) +
  theme_minimal() +  # Or any other suitable theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")  # Set bold face for axis titles
  ) +
  ggtitle("Northern fulmar pooled distribution") +
  xlab("Longitude") +
  ylab("Latitude")


# Plot the map
print(p1)


ind_merge_sf_2 <- st_transform(ind_merge_sf, crs = "ESRI:54030")

# Defining land----
land <- as(world, "Spatial")
land_sf <- st_as_sf(land)

shp <- ne_countries(type = "countries", 
                    scale = "medium",
                    returnclass = "sf") |>
  st_transform("ESRI:54030") 

# Overlaying shapefiles atop distribution 
p2 <- shp |> ggplot() + geom_sf() 
p2
p3 <- p2 + geom_sf(data = ind_merge_sf_2, aes(colour = colony), size = .1) 
p3

getwd()
ggsave("test.jpeg", plot = p1, width = 18, height = 12, units = "in", dpi = 300)
ggsave("test_2.jpeg", plot = p1, width = 6, height = 4, units = "in", dpi = 300)


t_PAME <- PAME_shapefile |> st_transform("ESRI:54030")
p4 <- t_PAME |> ggplot() + geom_sf(aes(color = "red"))
p4

# Moment of truth
t_PAME$LME_NAME
curated_list <- c("North Sea", "Faroe Plateau", "Iceland Shelf and Sea", "Barents Sea", "Canadian Eastern Arctic - West Greenland")
curated_PAME <- t_PAME %>% 
  filter(t_PAME$LME_NAME %in% curated_list,)
p5 <- p2 + geom_sf(data = curated_PAME, fill = "lightblue") 
# + geom_sf(data = ind_merge_sf, aes(col = colony), size = .1) 
p5 # with curated PAME regions displayed

multipolygon <- PAME_shapefile
pame <- multipolygon[1:15,c(3,10)]


# Defining a for loop results data frame----

long_results_df <- data.frame(matrix(ncol = 4))
colnames(long_results_df) <- c("individ_id", "lon", "lat", "Location")

# for loop for all individuals----
for(i in 1:length(indiv_merged_df$individ_id)){ # First for loop starts
  
for(j in pame$LME_NAME){ # Second for loop starts
  sub_2 <- pame %>% filter(LME_NAME == j) 
  point <- st_sfc(st_point(c(indiv_merged_df$lon[i], indiv_merged_df$lat[i])), crs = st_crs(sub_2))
  
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
    
    to_bind <- data.frame("individ_id" = indiv_merged_df$individ_id[i], "lon" = indiv_merged_df$lon[i], "lat" = indiv_merged_df$lat[i], "Location" = j)
    long_results_df <- rbind(long_results_df, to_bind)
    
  } # Second if loop ends
  
  else{ # First else loop starts
    long_results_df <- long_results_df 
  } # First else loop ends
  
  
} # Second for loop ends
} # First for loop ends

View(long_results_df)

long_results_df[46894,]
nrow(indiv_merged_df) - 57746
indiv_merged_df[indiv_merged_df$lon == "-9.5195" & indiv_merged_df$lat == "60.9932",]
# So start next loop with the row 58,676
final_df <- merge(indiv_merged_df, long_results_df[-1,], by = c("individ_id", "lon", "lat"))
final_df_clean <- na.omit(final_df)

write.csv(long_results_df, "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/second_hyp_tillnow_df.csv", row.names = F)

# What I want is to represent each ind track position to be listed inside one of the PAME_shapefile$LME_NAME; i.e. mutate a separate column 
# with an ifelse condition inside the brackets


# final_df_clean is to be used for calculating kernels based on a monthly basis (just copy its script and replace colony with Location)

# New for loop starting from row 58,676 of indiv_merged_df----

long_results_4 <- data.frame(matrix(ncol = 4))
colnames(long_results_4) <- c("individ_id", "lon", "lat", "Location")

for(i in 159618:length(indiv_merged_df$individ_id)){ # First for loop starts
  
  for(j in pame$LME_NAME){ # Second for loop starts
    sub_2 <- pame %>% filter(LME_NAME == j) 
    point <- st_sfc(st_point(c(indiv_merged_df$lon[i], indiv_merged_df$lat[i])), crs = st_crs(sub_2))
    
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
      
      to_bind <- data.frame("individ_id" = indiv_merged_df$individ_id[i], "lon" = indiv_merged_df$lon[i], "lat" = indiv_merged_df$lat[i], "Location" = j)
      long_results_4 <- rbind(long_results_4, to_bind)
      
    } # Second if loop ends
    
    else{ # First else loop starts
      long_results_4 <- long_results_4 
    } # First else loop ends
    
    
  } # Second for loop ends
} # First for loop ends


View(long_results_4)
nrow(long_results_4)
(long_results_4)[nrow(long_results_4),]
write.csv(long_results_4, "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/second_hyp_tillnow_4_df.csv", row.names = F)

indiv_merged_df[indiv_merged_df$lon == "-5.9478" & indiv_merged_df$lat == "59.4272",]
# long_results_2 started from 58676 row; long_results_3 starts from 89080 row then;
# long_results_4 starts from 159618 row

# Do it for the entire dataset at once instead of doing it line by line. Use st_join
# For Hyp 3, calculate kernels according to the nbs year

# New_script using st_join----

# new_result_1 <- data.frame(matrix(ncol = 5))
# colnames(new_result_1) <- c("individ_id", "lon", "lat", "timestamp", "Location")

# Creating an sf object of the main input dataframe
sf_df <- indiv_merged_df %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Making polygon geometry valid
valid_pame <- st_make_valid(PAME_shapefile)

# Finding which points lie inside which polygons
joined <- st_join(sf_df, valid_pame, join = st_intersects)
joined_df <- as.data.frame(joined)
condensed_joined_df <- joined_df %>% dplyr::select(c("individ_id", "timestamp", "LME_NAME"))
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/")
write.csv(condensed_joined_df, "condensed_joined_df_with_pame_regions.csv")

# Defining relevant regions for hypothesis
relevant_regions <- c("North Sea", "Faroe Plateau", "Iceland Shelf and Sea", "Canadian Eastern Arctic - West Greenland", "Barents Sea")
# Added this line to get a comprehensive list
relevant_regions <- unique(condensed_joined_df$LME_NAME)

# Pasting month_kernel script here----

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/"
mylocs <- readRDS(paste0(datadir, "SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"summaryTable.rds"))

# Exploring data
n_birds_tracked_per_year <- mylocs %>%
  dplyr::left_join(summary_info, join_by(ring)) %>%
  dplyr::mutate(year = year(timestamp)) %>%
  dplyr::group_by(colony,year) %>%
  dplyr::summarise(n_birds_tracked = n_distinct(ring))
write.csv(n_birds_tracked_per_year, 
    "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/n_birds_tracked_per_year.csv", row.names = F)

# Merging----
indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
names(indiv_merged_df)[1] <- "individ_id"
new_indiv_merged_df <- merge(indiv_merged_df, condensed_joined_df, by = c("individ_id", "timestamp"))
write.csv(new_indiv_merged_df, "final_df_for_hyp_2.csv")

nbs_mylocs <- new_indiv_merged_df %>% 
  # filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"

# Creating an sf object----
nbs_mylocs_sf <- nbs_mylocs %>% 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) 

# Assorting by months----

all_data <- nbs_mylocs
df <- all_data[!is.na(all_data$timestamp), ]
df <- df %>% dplyr::mutate(month = month(timestamp))
months <- sort(unique(df$month))

# Setting a directory for month-wise output rasters
dir_kernels <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs"

# Defining land----

land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# In case, one wants to exactly replicate what the petrels' paper did
# setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/")
# land <- sf::read_sf(dsn = "baselayer", layer = "ne_10m_land")
# plot(land)
# proj_wgs84 <- crs(land)

ss_for_kernel_cumulative_df <- data.frame(matrix(ncol = 3))
colnames(ss_for_kernel_cumulative_df) <- c("Region", "Month", "Sample_size_for_kde")

# for loop for computing input hr rasters for multiplication---- 
for(i in relevant_regions){ # First for loop begins
  sub <- df %>% filter(LME_NAME == i)
  
  for(month_number in 1:length(months)){ # Second for loop begins
    tracks_wgs <- sub[sub$month == months[month_number],] 
    
    sub <- subset(sub, lat < 90)  
    png(filename = paste0(dir_kernels, "/distributions/", i,"_",month_number,".png"))  
    plot(lat~lon, data = sub, type= "n", asp = 1, 
         xlim=c(min(all_data$lon) + 0.2,max(all_data$lon) - 0.2), 
         ylim=c(min(all_data$lat) + 0.5, max(all_data$lat) - 0.5), 
         main="", frame = T, xlab="", ylab="")
    plot(land, col='lightgrey', add= T)
    points(lat~lon, data = sub, pch = 16, cex = 0.5, col="blue")
    points(lat~lon, data = tracks_wgs, pch = 16, cex = 0.5, col="red")
    dev.off()
    
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
      tracks$month <- factor(tracks@data$month)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"month"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      
      vud <- adehabitatHR::getvolumeUD(kudl)
      
      fud <- vud[[1]]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- as.numeric(hr95 <= 95)
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
      
      raster::writeRaster(rast_mask_final, filename = paste0(datadir,"Second_hyp_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      mask_wgs84 <- projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
      
      KDE_ref <- paste0(i, "_", months[month_number])
      
      raster::writeRaster(mask_wgs84, filename = paste0(dir_kernels,"/unique_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      png(filename = paste0(dir_kernels,"/unique_distributions/",KDE_ref,".png"))
      plot(mask_wgs84, main = paste0("",i))
      plot(land, add = T, col = "#66000000")
      dev.off()
      
      # ss_for_kernel <- data.frame(Region = i, Month = month_number, Sample_size_for_kde = nrow(tracks_wgs))
      # ss_for_kernel_cumulative_df <- rbind(ss_for_kernel_cumulative_df, ss_for_kernel)
      
      
    } # First if loop ends
  } # Second for loop ends
} # First for loop ends

# Canadian Eastern Arctic- West Greenland is missing. Because there was just 1 observation in the year 2021 for March. 

View(ss_for_kernel_cumulative_df)
ss_for_kernel <- ss_for_kernel_cumulative_df[-1,]
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/")
write.csv(ss_for_kernel, "Sample_size_for_KDE.csv")



# Pasting month multiplication script now

# Resetting mapping parameters to default
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)

# Loading data----
dir_demClasses <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/Second_hyp_tifs/" # input directory
dir_1by1 <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/multiplication_rasters" # output directory
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
  
  if(sum(na_over@data@values)> 0){
    png(paste0(dir_1by1,"/na_maps/",name,".png"), width=1399,height=455)
    par(mfrow=c(1,2))
    plot(a_proj2,main=name,col=colsviri,legend=F)
    plot(na_over,main=paste0("Exposure Nas = ",sum(na_over@data@values)),
         col=colsinf,legend=F)
    dev.off()
  }
  
  name_split <- strsplit(name,"_")[[1]]
  month <- name_split[length(name_split)]
  population <- paste(name_split[1])
  
  check <- cbind(b,c)
  result <- cbind(population, month, check, exposure_score)
  dat <- rbind(dat, as.data.frame(result))
  print(i)
}

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/")
write.csv(dat, "correct_exposure_scores_by_month.csv",
          row.names = F)  
nas$percent_na <- nas$nas/nas$vals*100
View(nas)
write.csv(nas, "correct_nas_exposure_scores_by_month.csv",
          row.names = F)  


#2 don't exist ----------this isn't required in our case
nas_no_na <- subset(nas,vals > 0)
nas_no_na
mean(nas_no_na$percent_na)

exposure_score_csv <- read.csv("correct_exposure_scores_by_month.csv")
pop_exposure <- exposure_score_csv %>%
  group_by(population) %>%  
  summarise(population_exposure = round(median(exposure_score), 4)) %>%
  data.frame() 
write.csv(pop_exposure, "correct_exposure_scores_by_population.csv",
          row.names = F) 
Species_exposure_score <- mean(pop_exposure$population_exposure)
Species_exposure_score 

new_df_1 <- exposure_score_csv[,-c(3,4)]
colnames(new_df_1) <- c("Region", "Month", "pers")

new_df_2 <- read.csv("Sample_size_for_KDE.csv")
new_df_2$Region <- gsub(" ", ".", new_df_2$Region)
new_df_2 <- new_df_2[,-1]

new_df_2$Region[new_df_2$Region == "Celtic-Biscay.Shelf"] <- "Celtic.Biscay.Shelf" 
new_df_2$Region[new_df_2$Region == "Canadian.Eastern.Arctic.-.West.Greenland"] <- "Canadian.Eastern.Arctic...West.Greenland"
new_df_2$Region[new_df_2$Region == "Labrador.-.Newfoundland"] <- "Labrador...Newfoundland" 
new_df_2$Region[new_df_2$Region ==  "Northeast.U.S..Continental.Shelf"] <- "Northeast.U.S..Continental.Shelf" 

new_df_3 <- merge(new_df_1, new_df_2, by = c("Region", "Month"))
nas[c("Region", "Month")] <- do.call(rbind, strsplit(nas$name, "_", fixed = T))
new_df_4 <- nas[,-1]
new_df_5 <- merge(new_df_3, new_df_4, by = c("Region", "Month"))
write.csv(new_df_5, "Correct_base_df.csv")

correct_base_analysis_df <- new_df_5 %>% group_by(Region) %>% 
  summarise(median_percent_na = median(percent_na), median_pers = median(pers))
write.csv(correct_base_analysis_df, "correct_base_analysis_df.csv")

# Correlation test now!!!!----

# Creating a dataset for EcoQO values
relevant_regions <- c("North Sea", "Faroe Plateau", "Iceland Shelf and Sea", "Canadian Eastern Arctic - West Greenland", "Barents Sea")
ecoqo_df <- data.frame(population = relevant_regions, ecoqo = c(51, 40.5, 27.6, 14, 22.5))
ecoqo_df$population <- gsub(" ",".", ecoqo_df$population)

pop_exposure$population[2] <- "Canadian.Eastern.Arctic.-.West.Greenland"

normality_test_1 <- shapiro.test(pop_exposure$population_exposure)  
print(normality_test_1) # non-normal univariate distribution 
normality_test_2 <- shapiro.test(ecoqo_df$ecoqo)
print(normality_test_2) # non-normal univariate distribution 

analysis_df <- merge(pop_exposure, ecoqo_df, by = "population")
View(analysis_df)

corr_output <- cor.test(analysis_df$population_exposure, analysis_df$ecoqo, method = "kendall")
print(corr_output) # can be explained by percent_nas

# Controlling for percentage nas----
nas_df <- read.csv("nas_exposure_scores_by_month.csv")
nas_df[c("population", "month")] <-  do.call(rbind, strsplit(as.character(nas_df$name), "_")) 
View(nas_df)

median_nas_df <- nas_df %>% group_by(population) %>% summarize(median_percent_na <- median(percent_na))
median_nas_df$population[2] <- "Canadian.Eastern.Arctic.-.West.Greenland"
write.csv(median_nas_df, "median_nas_df.csv")
colnames(median_nas_df) <- c("Region","Median_percent_nas")
analysis_df_2 <- merge(analysis_df, median_nas_df, by = "population")
colnames(analysis_df_2) <- c("region", "pers", "ecoqo", "percent_nas")

corr_output_2 <- spcor.test(analysis_df_2$ecoqo, analysis_df_2$pers, analysis_df_2$percent_nas, method = "kendall")
corr_output_3 <- cor.test(analysis_df_2$ecoqo, analysis_df_2$pers, method = "kendall")

print(corr_output_2)
print(corr_output_3)

# New corrected correlation test

colnames(ecoqo_df) <- c("Region", "EcoQO value")
correct_base_analysis_df$Region[correct_base_analysis_df$Region ==  "Canadian.Eastern.Arctic...West.Greenland"] <- "Canadian.Eastern.Arctic.-.West.Greenland"
correct_hyp_2_analysis_df <- merge(ecoqo_df, correct_base_analysis_df, by = "Region")
write.csv(correct_hyp_2_analysis_df, "correct_hyp_2_analysis_df.csv")
correct_correlation_test <- cor.test(correct_hyp_2_analysis_df$`EcoQO value`, correct_hyp_2_analysis_df$median_pers, method = "kendall")

print(correct_correlation_test)

for_corr <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/correct_hyp_2_analysis_df.csv")
cor.test(for_corr$EcoQO.value, for_corr$median_pers, method = "kendall")
