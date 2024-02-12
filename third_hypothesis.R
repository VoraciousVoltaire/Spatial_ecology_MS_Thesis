# Home range overlap script using the amt package

# Clearing environment----
remove(list = ls())

# Loading essential packages----
library(amt)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(sp)
library(raster)
library(spData)
library(raster)
library(terra)
library(sf)
library(cowplot)
library(stringr)
library(RColorBrewer)
library(viridisLite)
library(viridis)
library(ppcor) # for partial correlation test

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
mylocs <- readRDS(paste0(datadir, "/SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"/summaryTable.rds"))

# Merging and curating for NBS (Non-Breeding Season)----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"

# Creating a new combined dataset which merges tracks from Skalfandi and Langanes (distance between colonies: 130.429 km)
colonies_to_combine <- c("Skjalfandi", "Langanes")
nbs_mylocs$colony[nbs_mylocs$colony %in% colonies_to_combine] <- "Combined"

# Defining land----
land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# Calculating individual home ranges using the amt package----

results_df <- data.frame(matrix(ncol = 2))
colnames(results_df) <- c("Colony", "Overlap_score")

for(i in unique(nbs_mylocs$colony)){ # First for loop start
  sub <- as.data.frame(nbs_mylocs) %>% filter(colony == i) 
  
  # Setting a colony-centered crs 
  median_loc <- cbind(median(sub$lon), median(sub$lat))
  DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",median_loc[1]," +lat_0=",median_loc[2]," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m ", sep="")) 
  
NF_tracks <- sub %>% 
  make_track(lon, lat, timestamp, id = individ_id, 
  crs = 4326)
attr(NF_tracks, "crs") <- DgProj
NF_tracks <- NF_tracks %>% nest(data = -"id") %>% 
  arrange(id)



# # Checking mcp overlap using the hr method
# 
# NF_tracks_mcp <- NF_tracks %>% 
#   mutate(mcp = map(data, function(x) 
#     x %>% hr_mcp(levels = c(1.0))))
# 
# NF_tracks_mcp_overlap <- hr_overlap(NF_tracks_mcp$mcp,
#                                 labels = NF_tracks$id, 
#                                 which = "all", 
#                                 # alternative which = "consecutive",
#                                 # "one_to_all"
#                                 conditional = FALSE)
# # the only catch is that the values are directional in nature,
# # that is they depend on which home range is mapped onto which 

# Using kernel density estimation----

# Make a trast that isn't nested
NF_nonest_tracks <- sub %>% 
  make_track(lon, lat, timestamp, id = individ_id, 
             crs = 4326)
attr(NF_nonest_tracks, "crs") <- DgProj


# Making a base raster because each Kernel Density Estimator will likely use different pixel sizes and 
# output of hr_overlap will otherwise be blank

base_trast <- make_trast(NF_nonest_tracks, res = 50)

NF_tracks_KDE <- NF_tracks %>% 
  mutate(kde = map(data, function(x) 
    x %>% hr_kde(trast = base_trast, levels = 0.95)))

NF_tracks_kde_overlap <- hr_overlap(NF_tracks_KDE$kde,
                                    labels = NF_tracks$id,
                                    which = "all",
                                    type = "ba",
                                    conditional = F)
range(NF_tracks_kde_overlap$overlap)

output_df <- NF_tracks_kde_overlap %>% group_by(from) %>%
  summarize(mean_overlap = mean(overlap))

colony_hr_kde_overlap_index <- median(output_df$mean_overlap) 
to_bind <- data.frame("Colony" = i, "Overlap_score" = colony_hr_kde_overlap_index)
results_df <- rbind(results_df, to_bind)

} # First for loop ends

View(results_df)
kde_final_overlap <- results_df[-1,]
kde_final_overlap$Colony <- gsub(" ",".",kde_final_overlap$Colony)

# Saving individual home range rasters for raster multiplication with the plastics raster----

# You're screwing up the crs part that's why the rasters aren't getting trimmed 

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/")

for(i in unique(nbs_mylocs$colony)){ # First for loop start
  sub <- as.data.frame(nbs_mylocs) %>% filter(colony == i) 
  
  # Setting a colony-centered crs 
  median_loc <- cbind(median(sub$lon), median(sub$lat))
  DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",median_loc[1]," +lat_0=",median_loc[2]," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m ", sep="")) 
  
  NF_tracks <- sub %>% 
    make_track(lon, lat, timestamp, id = individ_id, 
               crs = DgProj) %>% nest(data = -"id") %>% arrange(id)
  # attr(NF_tracks, "crs") <- DgProj # didn't work

  
  # crs couldn't be set here- error lies in this step. attr(NF_tracks, "crs") yields a null response
  # attr(NF_tracks, "crs") <- DgProj
  
  NF_nonest_tracks <- sub %>% 
    make_track(lon, lat, timestamp, id = individ_id, 
               crs = DgProj)
  # attr(NF_nonest_tracks, "crs") <- DgProj
  
  base_trast <- make_trast(NF_nonest_tracks, res = 50) # Spending too much time here; shall look into the res argument followed by the crs thingy and getting the masked plot appear
  
  NF_tracks_KDE <- NF_tracks %>% 
    mutate(kde = map(data, function(x) 
      x %>% hr_kde(trast = base_trast, levels = 0.95)))
  
  for(j in 1:length(unique(sub$individ_id))){ # Second for loop starts
  
  rast <- raster(NF_tracks_KDE$kde[[j]]$ud) # Doubtful step; because I'm getting empty rasters as outputs, messing up here
  
  # Proportionate scaling of home range rasters
  rast[is.na(rast)] <- 0
  rast <- rast/sum(raster::getValues(rast))
  rast[rast == 0] <- NA
  
  # Cropped extent
  x.matrix <- is.na(as.matrix(rast))
  colNotNA <- which(colSums(x.matrix) != nrow(rast))
  rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
  croppedExtent <- raster::extent(rast, 
                                  r1 = rowNotNA[1], 
                                  r2 = rowNotNA[length(rowNotNA)],
                                  c1 = colNotNA[1], 
                                  c2 = colNotNA[length(colNotNA)])
  cropped <- raster::crop(rast, croppedExtent)
  cropped[is.na(rast)] <- 0
  
 #  # Cropping extent
 #  rast[rast == 0] <- NA
 #  x.matrix <- is.na(as.matrix(rast))
 #  colNotNA <- which(colSums(x.matrix) != nrow(rast))
 #  rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
 # rowNotNA
 #  croppedExtent <- raster::extent(rast, 
 #                                  rowNotNA[1] - 2, 
 #                                  rowNotNA[length(rowNotNA)] + 2,
 #                                  colNotNA[1] - 2, 
 #                                  colNotNA[length(colNotNA)] + 2)
 #  cropped <- raster::crop(rast, croppedExtent)
 #  cropped[is.na(rast)] <- 0
  
  # Changing land's projection
  mask_proj <- sp::spTransform(land, DgProj) # changing projection to EPSG 3035
  mask_proj_pol <- as(mask_proj, "SpatialPolygons") 
  
  ## set to NA cells that overlap mask (land)
  rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = T)
  trial <- projectRaster(rast, crs = proj_wgs84)
  
  
  rast_mask <- rast_mask_na
  rast_mask[is.na(rast_mask)] <- 0
  rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
  rast_mask[rast_mask == 0] <- NA
  rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE) # this is to be used as the input whilst doing raster multiplication
  rast_mask_final
  
  #PLOT & SAVE ####
  mask_wgs84 <- raster::projectRaster(rast_mask_final,  crs = proj_wgs84, over = F)
  raster::writeRaster(mask_wgs84, filename = paste0("input_rasters/",i,"_", j ,".tif"),
                      format = "GTiff", overwrite = T)      
  
  # PLOT
  png(filename = paste0("output/kde_plots/",i,"_",j,".png"))
  plot(mask_wgs84, main = paste0(i,"_",j))
  plot(land, add = T, col = "#66000000")
  dev.off()
  
  } # Second for loop ends
  } # First for loop ends

# Now, just the multplication script has to be recreated by simply changing the working directory




# Multiplication script----

# Defining and loading data

wd <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp" # setting a working directory
dir_1by1 <- paste0(wd, "/output/multiplication_rasters/") # output directory
dir_demClasses <- paste0(wd, "/input_rasters/") # input directory which contains 10km^2 resolution rasters
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data/") # directory where plastics raster is stored
plastics <- raster::raster("00_PlasticsRaster.tif")

r <- raster()

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

# For loop starts here

dat <- data.frame()
result <- c()

for (i in 1:length(files)){
  
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
  
  raster::writeRaster(a_proj2, filename = raster_name_2, format = "GTiff", overwrite = TRUE)
  
  over <- a_proj2 * p_sum1
  
  over_score <- over
  summary(over_score@data@values)
  over_score[is.na(over_score)] <- 0
  exposure_score <- round(sum(raster::getValues(over_score))*1000000,4)
  
  png(paste0(dir_1by1,"maps/",name,".png"), width=1399,height=455)
  par(mfrow=c(1,2))
  plot(a_proj2, main=paste0(name," distribution"), col = colsviri,legend = F)
  plot(over, main = paste0("Exposure score = ", exposure_score),
       col = colsinf,legend = F)
  
  dev.off()
  
  if(sum(na_over@data@values)> 0){
    png(paste0(dir_1by1,"na_maps/",name,".png"), width = 1399, height = 455)
    par(mfrow=c(1,2))
    plot(a_proj2, main = name, col = colsviri, legend = F)
    plot(na_over, main = paste0("Exposure Nas = ", sum(na_over@data@values)),
         col = colsinf, legend = F)
    dev.off()
  }
  
  name_split <- strsplit(name,"_")[[1]]
  name_split
  individual <- name_split[2]
  population <- paste(name_split[1])
  
  check <- cbind(b,c)
  result <- cbind(population, individual, check, exposure_score)
  dat <- rbind(dat, as.data.frame(result))
  print(i)
} 

setwd(paste0(wd, "/output/csv/"))
write.csv(dat, "exposure_scores_by_individual.csv",
          row.names = F)  
nas$percent_na <- nas$nas/nas$vals*100
head(nas)

exposure_score_csv <- read.csv(paste0(wd, "/output/csv/exposure_scores_by_individual.csv"))
pop_exposure <- exposure_score_csv %>%
  group_by(population) %>%  
  summarise(population_exposure = round(mean(exposure_score), 4)) %>%
  data.frame() 
setwd(paste0(wd,"/csv/"))
write.csv(pop_exposure, "exposure_scores_by_population.csv",
          row.names = F) 
Species_exposure_score <- mean(pop_exposure$population_exposure)
Species_exposure_score 


# The final step----
variance_df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/exposure_scores_by_individual.csv")
to_combine <- c("Skjalfandi", "Langanes")
variance_df$population[variance_df$population %in% to_combine] <- "Combined"

variance_output <- variance_df %>% group_by(population) %>% summarize(variance = var(exposure_score)) 
colnames(variance_output)[colnames(variance_output) == "population"] <- "Colony"

# Correlation test----
corr_df <- merge(kde_final_overlap, variance_output, by = "Colony")
write.csv(corr_df, "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/ba_third_hyp_corr_test.csv", row.names = F)

correlation <- cor.test(corr_df$Overlap_score, corr_df$variance, method = "kendall")
print(correlation)

# Trying without Alkefjellet
corr_df_2 <- corr_df[-1,]
View(corr_df_2)
correlation_2 <- cor.test(corr_df_2$Overlap_score, corr_df_2$variance, method = "kendall")
print(correlation_2)

# Trying to control for sample size and number of individual tracks----
ss_by_colony <- nbs_mylocs %>% group_by(colony) %>% summarize(ss = length(unique(individ_id)))
colnames(ss_by_colony) <- c("Colony", "Sample_size")
ss_by_colony$Colony <- gsub(" ",".",ss_by_colony$Colony)

ss_by_ind <- nbs_mylocs %>% group_by(colony) %>% summarize(n_tracks = length(individ_id))
colnames(ss_by_ind) <- c("Colony", "Total_number_of_tracks")
ss_by_ind$Colony <- gsub(" ",".",ss_by_ind$Colony)

par_list <- list(kde_final_overlap, variance_output, ss_by_colony, ss_by_ind)
par_df <- par_list %>% reduce(full_join, by = "Colony")
View(par_df)

partial_correlation <- pcor.test(par_df$Overlap_score, par_df$variance, c(par_df[,c(4:5)]), method = "kendall") # controlling for both, sample size and number of total tracks
print(partial_correlation)

# par_corr_1 <- pcor(par_df[,-1], method = "kendall") # gives pair-wise partial correlation p values
# print(par_corr_1)

par_corr_2 <- pcor.test(par_df$Overlap_score, par_df$variance, par_df$Sample_size, method = "kendall") # controlling just for sample size
print(par_corr_2)

par_corr_3 <- pcor.test(par_df$Overlap_score, par_df$variance, par_df$Total_number_of_tracks, method = "kendall") # controlling just for number of tracks
print(par_corr_3)

# Creating a new dataframe with mean number of tracks instead of total number of tracks
to_join <- nbs_mylocs %>% group_by(individ_id) %>% summarise(n = n())
merged_n_tracks <- merge(nbs_mylocs, to_join, by = "individ_id")
df <- merged_n_tracks %>% group_by(individ_id) %>% summarise(number_of_tracks = unique(n))
merged_df <- merge(df, nbs_mylocs[,c(1,7)], by = "individ_id")
unique_merged_df <- unique(merged_df)
mean_n_tracks <- unique_merged_df %>% group_by(colony) %>% summarise(mean = mean(number_of_tracks))
colnames(mean_n_tracks) <- c("Colony", "Mean_total_number_of_tracks")
mean_n_tracks$Colony <- gsub(" ",".",mean_n_tracks$Colony)

new_par_list <- list(kde_final_overlap, variance_output, ss_by_colony, mean_n_tracks)
new_par_df <- new_par_list %>% reduce(full_join, by = "Colony")
View(new_par_df)

new_partial_correlation <- spcor.test(new_par_df$Overlap_score, new_par_df$variance, new_par_df$Mean_total_number_of_tracks, method = "kendall")
print(new_partial_correlation) # controlling for mean_total_number of tracks had no effect

# No significant relationship 

# This concludes the script for obtaining individual rasters and multiplying them with the plastics raster using the amt package----

