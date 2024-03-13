# Hypothesis 3 but kernels calculated on a sampling year basis

# first testing the correlation between the values calculated on the basis of amt and adehabitathr kde overlap functions----

input_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv"

amt_overlap <- read.csv(paste0(input_dir, "/ba_third_hyp_corr_test.csv"))[,1:2]
colnames(amt_overlap) <- c("colony", "amt_overlap_score")
View(amt_overlap)

adehabitathr_overlap <- read.csv(paste0(input_dir, "/third_hyp_pure_ade_corr_test.csv"))[,1:2]
colnames(adehabitathr_overlap) <- c("colony", "adehabitathr_overlap_score")
View(adehabitathr_overlap)

analysis_df <- merge(amt_overlap, adehabitathr_overlap, by = "colony")

normality_test_1 <- shapiro.test(analysis_df$amt_overlap_score)
print(normality_test_1)
normality_test_2 <- shapiro.test(analysis_df$adehabitathr_overlap_score)
print(normality_test_2)

corr_output <- cor.test(analysis_df$amt_overlap_score, analysis_df$adehabitathr_overlap_score, method = "kendall")
print(corr_output)


# Loading essential packages----
library(lubridate)
library(tidyverse)
library(dplyr)
library(sf)
library(sp)
library(viridisLite)
library(viridis)
library(RColorBrewer)
library(raster)
library(spData)
library(adehabitatHR)

# Pasting month_kernel script----

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/"
mylocs <- readRDS(paste0(datadir, "SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"summaryTable.rds"))

# # Exploring data
# mylocs %>%
#   dplyr::left_join(summary_info, join_by(ring)) %>%
#   dplyr::mutate(year = year(timestamp)) %>%
#   dplyr::group_by(colony,year) %>%
#   dplyr::summarise(n_birds_tracked = n_distinct(ring))

# Creating an sf object----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
dim(indiv_merged_df)
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"
dim(nbs_mylocs)
nbs_mylocs_sf <- nbs_mylocs %>% 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) 

# Assorting by months----

all_data <- nbs_mylocs
df <- all_data[!is.na(all_data$timestamp), ]
df <- df %>% dplyr::mutate(month = month(timestamp))
months <- sort(unique(df$month))

# Adding a column for nbs tracking year----

df_mod <- df %>% mutate(tracking_year = ifelse(month < 7, year - 1, year))
nbs_years <- sort(unique(df_mod$tracking_year))
unique(df_mod$colony)
# Defining land----

land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# Doing a separate for loop for calculating kernel overlap on a tracking_year basis----

results_df <- data.frame(matrix(ncol = 3))
colnames(results_df) <- c("Colony", "udoi_Overlap_score", "ba_Overlap_score")

to_combine <- c("Skjalfandi", "Langanes")
df_mod$colony[df_mod$colony %in% to_combine] <- "Combined"

# Define custom function to divide each row by its maximum value
divide_by_max <- function(row) {
  max_value <- max(row)
  result <- row / max_value
  return(result)
}

for(i in unique(df_mod$colony)){ # First for loop begins
  sub <- df_mod %>% filter(colony == i)
  
  for(j in 1:length(unique(sub$tracking_year))){ # Second for loop begins
    tracks_wgs <- sub[sub$tracking_year == unique(sub$tracking_year)[j],]
    
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
      tracks$tracking_year <- factor(tracks@data$tracking_year)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"individ_id"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      
      
if(length(attributes(kudl)$names) > 1){ # Second if loop starts
  
  # Calculating home range overlap using kerneloverlaphr----
  ade_kde_overlap <- kerneloverlaphr(kudl, meth = "UDOI", percent = 95, conditional = F)
  # range(as.data.frame(ade_kde_overlap)) # 0 to 0.9998733
  df <- as.data.frame(ade_kde_overlap)
  df_2 <- apply(df, 1, divide_by_max)
  print(range(df_2))
  median_df <- apply(df_2, 1, median, na.rm = T)
  median_overlap_df <- as.data.frame(median_df)
  colnames(median_overlap_df) <- "overlap"
  colony_hroverlap <- median(median_overlap_df$overlap) 
  
  ade_kde_overlap_2 <- kerneloverlaphr(kudl, meth = "BA", percent = 95, conditional = F)
  df_ba <- as.data.frame(ade_kde_overlap_2)
  df_ba_2 <- apply(df_ba, 1, median, na.rm = T)
  median_overlap_df_ba <- as.data.frame(df_ba_2)
  colnames(median_overlap_df_ba) <- "overlap"
  colony_hroverlap_ba <- median(median_overlap_df_ba$overlap) 
  
  to_bind <- data.frame("Colony" = paste0(i,"_",(unique(sub$tracking_year))[j]), "udoi_Overlap_score" = colony_hroverlap, "ba_Overlap_score" = colony_hroverlap_ba)
  results_df <- rbind(results_df, to_bind)

  
} # Second if loop ends
} # First if loop ends
  } # Second for loop ends
} # First for loop ends
  
View(results_df)
kde_final_overlap <- results_df[-1,]
View(kde_final_overlap)
write.csv(kde_final_overlap, "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/final_new_third_hyp_3_with_ba.csv")

range(kde_final_overlap$udoi_Overlap_score) #  0.3009095 0.9449845
range(kde_final_overlap$ba_Overlap_score) # 0.6639227 0.9944815

install.packages("psych")
library(psych)

kde_final_overlap[c("Colony", "Tracking_year")] <- do.call(rbind, str_split(as.character(kde_final_overlap$Colony), pattern = "_"))
                                                        
harmonic_kde_final_overlap <- kde_final_overlap %>% group_by(Colony) %>%
  summarise(udoi = harmonic.mean(udoi_Overlap_score), ba = harmonic.mean(ba_Overlap_score))
harmonic_kde_final_overlap$final_score <- apply(harmonic_kde_final_overlap[, c(2, 3)], 1, harmonic.mean)

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/")
write.csv(harmonic_kde_final_overlap, "harmonic_hr_overlap.csv")

hyp_1_input <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/input_from_hyp_1/correct_ind_pers_by_population.csv")
names(hyp_1_input)[1] <- "Colony"
harmonic_kde_final_overlap$Colony <- gsub(" ", ".", harmonic_kde_final_overlap$Colony)
merged_df <- merge(hyp_1_input[-c(1,2),], harmonic_kde_final_overlap, by = "Colony")
write.csv(merged_df, "harmonic_overlap_complete_df.csv")
View(merged_df)
cor.test(merged_df$pers, merged_df$udoi, method = "kendal")

# Just checking a hunch----

names(final_df)[1] <- "Colony" # final_df borrowed by lat_exposure script
lat_merged_df <- merge(merged_df, final_df[,-2], by = "Colony")
View(lat_merged_df)
cor.test(lat_merged_df$lat, lat_merged_df$pers, method = "kendal")
lat_merged_df



merged_df_2 <- merge(hyp_1_input, harmonic_kde_final_overlap, by = "Colony")
merged_df_with_Alk_Bjo <- merge(merged_df_2, final_df[,-2], by = "Colony")
View(merged_df_with_Alk_Bjo)
write.csv(merged_df_with_Alk_Bjo, "complete_df_hyp_3_with_lat.csv")
cor.test(merged_df_with_Alk_Bjo$lat, merged_df_with_Alk_Bjo$pers, method = "kendal") # significant
cor.test(merged_df_with_Alk_Bjo[-c(1,2),]$lat, merged_df_with_Alk_Bjo[-c(1,2),]$final_score, method = "kendal") # close 










# This reminds me to add one more column for sample size of each estimation 

kde_final_overlap[c("Just_colony", "Tracking_year")] <- do.call(rbind, strsplit(as.character(kde_final_overlap$Colony), "_", fixed = T))
to_save <- kde_final_overlap[,-1] 

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/")
write.csv(to_save, "final_hyp_3_part_1.csv")

# Calculating a median according to just colony
new_df <- kde_final_overlap %>% group_by(Just_colony) %>% summarize(median_overlap_score <- median(Overlap_score))
colnames(new_df) <- c("Colony", "Overlap_score")
# new_df$Colony <- gsub(" ", ".", new_df$population)
write.csv(new_df, "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/latest_final_third_hyp_3_part_1.csv")
  
# Moment of truth
final_final_analysis_df <- merge(new_df, tracking_year_pers_var, by = "population")
cor.test(final_final_analysis_df$Overlap_score, final_final_analysis_df$variance_pers, method = "kendall")

# Incorporating sample-size in this:
# Adding a segment for sample size----

ss_df <- data.frame(matrix(ncol = 2))
colnames(ss_df) <- c("Colony", "Number_of_individuals")

for(i in unique(df_mod$colony)){ # First for loop begins
  
  sub <- df_mod %>% filter(colony == i)
  
  for(j in 1:length(unique(sub$tracking_year))){ # Second for loop begins
    tracks_wgs <- sub[sub$tracking_year == unique(sub$tracking_year)[j],]
    
    if(nrow(tracks_wgs) > 4){ # First if loop begins
        
        sample_size <- length(unique(tracks_wgs$individ_id))
        ss_vector_to_add <- data.frame("Colony" = paste0(i,"_",(unique(sub$tracking_year))[j]), "Number_of_individuals" = sample_size)
        ss_df <- rbind(ss_df, ss_vector_to_add)
      

} # First if loop ends
  } # Second for loop ends
} # First for loop ends

View(ss_df)
ss_df <- ss_df[-1,]
ss_df_plural <- ss_df[ss_df$Number_of_individuals != 1,]

analysis_df_part_1 <- merge(kde_final_overlap, ss_df_plural, by = "Colony")
View(analysis_df_part_1)

analysis_df_part_1[c("Just_colony", "Tracking_year")] <-  do.call(rbind, strsplit(as.character(analysis_df_part_1$Colony), "_")) 
condensed_analysis_df_part_1 <- analysis_df_part_1 %>% group_by(Just_colony) %>% summarise(median_number_of_individuals = median(Number_of_individuals), median_overlap_score = median(Overlap_score),
                                                                                           number_of_years_tracked = n())

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/")
View(analysis_df_part_1)
write.csv(analysis_df_part_1, "analysis_part_1.csv")
write.csv(condensed_analysis_df_part_1, "condensed_analysis_part_1.csv")



# Separate script for calculating pers on a tracking year basis and then calculating the variance for each tracking year and then considering this final output table 
# as the input for correlation script 

# ty_vector <- sort(unique(df_mod$tracking_year)) # for kernel script
datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/" # for kernel script

# # Function to adjust latitude values----
# adjust_latitudes <- function(raster_data) {
#   # Get latitude values
#   lat_values <- raster::getValues(raster_data)
#   
#   # Get the raster extent
#   extent <- raster::extent(raster_data)
#   
#   # Mask out invalid latitude values
#   lat_values[lat_values < extent@ymin] <- NA
#   lat_values[lat_values > extent@ymax] <- NA
#   
#   # Set adjusted latitude values back to the raster
#   raster::values(raster_data) <- lat_values
#   
#   return(raster_data)
# }

# # Function to project raster and adjust latitudes with error handling----
# project_and_adjust <- function(raster_data, crs) {
#   tryCatch(
#     {
#       raster_data <- projectRaster(raster_data, crs = crs, over = FALSE)
#       raster_data <- adjust_latitudes(raster_data)
#     },
#     warning = function(w) {
#       # Handle warnings (e.g., print or log the warning message)
#       print(paste("Warning:", w))
#     },
#     error = function(e) {
#       # Handle errors (e.g., print or log the error message)
#       print(paste("Error:", e))
#     }
#   )
#   return(raster_data)
# }
 
# Filtering out Alkefjellet 2020 from df_mod----
df_mod_2 <- df_mod[!(df_mod$colony == "Alkefjellet" & df_mod$tracking_year == 2020),]
df_mod_3 <- df_mod[!(df_mod$colony == "Alkefjellet" & df_mod$month == 10),]

df_mod_4 <- df_mod[df_mod$colony != "Alkefjellet",]


# for loop for KDE----

for(i in unique(df_mod_4$colony)){ # First for loop begins
  sub <- df_mod_4 %>% filter(colony == i)
  
 
  for(j in 1:length(unique(sub$tracking_year))){ # Second for loop begins
    tracks_nbs <- sub[sub$tracking_year == unique(sub$tracking_year)[j],]
    
    # Setting a null grid for kernel estimation 
    if(nrow(tracks_nbs) > 4){ # First if loop starts
      
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
      
      # Setting a colony-centered crs 
      median_loc <- cbind(median(tracks_nbs$lon), median(tracks_nbs$lat))
      DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",median_loc[1]," +lat_0=",median_loc[2]," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m ", sep="")) 
      
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
      mask_proj_pol
      
      rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = TRUE)
      rast_mask <- rast_mask_na
      rast_mask[is.na(rast_mask)] <- 0
      rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
      
      rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE)
      rast_mask_final2 <- rast_mask_final 
      
      mask_wgs84 <- projectRaster(rast_mask_final2, crs = proj_wgs84, over = FALSE)
     
      
      
      
      KDE_ref <- paste0(i, "_", ty_vector[j])
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
  
  
 
  # I'm in an attempt right now to change precisely this
  variance_df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/exposure_scores_by_individual.csv")
  to_combine <- c("Skjalfandi", "Langanes")
  variance_df$population[variance_df$population %in% to_combine] <- "Combined"
  
  variance_output <- variance_df %>% group_by(population) %>% summarize(variance = var(exposure_score)) 
  colnames(variance_output)[colnames(variance_output) == "population"] <- "Colony"
  
  # Correlation test----
  corr_df <- merge(new_df, variance_output, by = "Colony")
  no_alk_corr_df <- corr_df[-1,]
  
  bleak_correlation_1 <- cor.test(corr_df$Overlap_score, corr_df$variance, method = "kendall")
  bleak_correlation_2 <- cor.test(no_alk_corr_df$Overlap_score, no_alk_corr_df$variance, method = "kendall")

  print(bleak_correlation_1)  
  print(bleak_correlation_2)  

  # one suggestion is to calculate pers according to tracking years; and then run a correlation test
  # on the entire new df instead of taking a median value of variance and median value of overlap score
  
  
  
  
  
  
  
  # New multiplication script----
  
  # Prep for multiplication script----
  dir_demClasses <- paste0(datadir, "/input_rasters/tifs/") # input directory
  dir_1by1 <- paste0(datadir, "output/multiplication_rasters") # output directory
  
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
    tracking_year <- name_split[length(name_split)]
    population <- paste(name_split[1])
    
    check <- cbind(b,c)
    result <- cbind(population, tracking_year, check, exposure_score)
    dat <- rbind(dat, as.data.frame(result))
    print(i)
  }
  
  setwd(paste0(datadir, "output/csv/"))
  write.csv(dat, "exposure_scores_by_tracking_year.csv",
            row.names = F)  
  nas$percent_na <- nas$nas/nas$vals*100
  head(nas)
  write.csv(nas, "nas_exposure_scores_by_tracking_year.csv",
            row.names = F)  
  
  
  
  #2 don't exist ----------this isn't required in our case
  nas_no_na <- subset(nas,vals > 0)
  nas_no_na
  mean(nas_no_na$percent_na)
  
  setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/")
  ty_exposure_score_csv <- read.csv("exposure_scores_by_tracking_year.csv")
  normality_test_1 <- shapiro.test(ty_exposure_score_csv$exposure_score)
  print(normality_test_1) # p value 0.2989
  
  tracking_year_pers_var <- ty_exposure_score_csv %>%
    group_by(population) %>%  
    summarise(variance_pers = var(exposure_score)) %>%
    data.frame() 
  tracking_year_pers_var$population <- gsub("\\."," ", tracking_year_pers_var$population)

  
  write.csv(tracking_year_pop_exposure, "ty_exposure_scores_by_population.csv",
            row.names = F) 
  Species_exposure_score <- mean(pop_exposure$population_exposure)
  Species_exposure_score 

  # Moment of truth----
  analysis_part_1 <- read.csv("condensed_analysis_part_1.csv")
  names(analysis_part_1)[2] <- "population"
  
    
  analysis_df_2 <- merge(tracking_year_pers_var, analysis_part_1, by = "population")
  
  final_cor_test_hyp_3 <- cor.test(analysis_df_2$variance_pers, analysis_df_2$median_overlap_score, method = "kendall")
  print(final_cor_test_hyp_3)
  
  
  # Final moment of truth----
  df_part_1 <- new_df
  df_part_1$Colony <- gsub(" ",".",df_part_1$Colony)
 
  df_part_2 <- read.csv(paste0("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/csv/Final_hypothesis_1_analysis_df.csv"))
  df_part_3 <- df_part_2 %>% group_by(Colony) %>% summarise(variance_pers = var(pers), 
                                                            median_percent_na = median(percent_na))
 
  write.csv(df_part_3, "Variance_with_median_percent_na.csv")
  df_part_4 <- merge(df_part_1, df_part_3, by = "Colony")
 

  write.csv(df_part_4, "Base_analysis_df_for_hyp_3.csv")
  
  analysis_df_hyp_3 <- df_part_4[-c(1,2,8),]
  analysis_df_hyp_3[order(analysis_df_hyp_3$variance_pers),]
  cor.test(analysis_df_hyp_3$Overlap_score, analysis_df_hyp_3$variance_pers, method = "kendall") # still needs refining
  