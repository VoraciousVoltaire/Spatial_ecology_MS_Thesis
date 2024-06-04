# Unique individual multiplication script for calculating exposure risk scores on an individual basis
# and then taking a mean of those values to arrive at one single value for the colony

# Clearing environment
rm(list = ls())

# Resetting mapping parameters to default
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)

# Loading essential packages----
library(adehabitatHR)
library(sf)
library(sp)
library(raster)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(viridisLite)
library(viridis)

# Defining and loading data----

wd <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/Latest_right_attempt/" # setting a working directory
dir_1by1 <- paste0(wd, "/outputs/multiplication_rasters/") # output directory
dir_demClasses <- paste0(wd, "/outputs/input_tifs") # input directory which contains 10km^2 resolution rasters
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

dat <- data.frame()
result <- c()

# For loop starts here----

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

  if(sum(na_over@data@values)> 0){ # if loop starts
    png(paste0(dir_1by1,"na_maps/",name,".png"), width = 1399, height = 455)
    par(mfrow=c(1,2))
    plot(a_proj2, main = name, col = colsviri, legend = F)
    plot(na_over, main = paste0("Exposure Nas = ", sum(na_over@data@values)),
         col = colsinf, legend = F)
    dev.off()
  } # if loop ends
  
  name_split <- strsplit(name,"_")[[1]]
  name_split
  individual <- name_split[2]
  population <- paste(name_split[1])
  tracking_year <- paste(name_split[3])
  
  check <- cbind(b,c)
  result <- cbind(population, individual, tracking_year, check, exposure_score)
  dat <- rbind(dat, as.data.frame(result))
  print(i)
} 

setwd(paste0(wd, "/outputs/csv/"))
write.csv(dat, "correct_exposure_scores_by_individual_according_to_tracking_yr.csv",
          row.names = F)  
nas$percent_na <- nas$nas/nas$vals*100
head(nas)

write.csv(nas, "correct_nas_ind.csv", row.names = F)
trial_nas <- nas
trial_nas[c("population", "individual", "tracking_year")] <- do.call(rbind, str_split(as.character(nas$name), pattern = "_"))
exposure_score_csv <- read.csv(paste0(wd, "/outputs/csv/correct_exposure_scores_by_individual_according_to_tracking_yr.csv"))

merged_df <- merge(exposure_score_csv, trial_nas, by = c("population", "individual", "tracking_year"))
View(merged_df)


variance_pers_ind_tracking_yr <- merged_df %>%
  group_by(population, individual) %>%
  summarise(number_of_years_tracked = n(), intra_ind_variance_pers = var(exposure_score),
            intra_ind_median_pers = median(exposure_score),
            intra_ind_percent_nas = median(percent_na)) %>%
  data.frame()
write.csv(variance_pers_ind_tracking_yr, "var_pers_by_ind_tracking_yr.csv",
          row.names = F) 

ind_pop_exposure <- variance_pers_ind_tracking_yr %>%
  group_by(population) %>%  
  summarise(pers = median(intra_ind_median_pers), 
           median_intra_ind_variance_pers = median(na.omit(intra_ind_variance_pers)),
           variance_pers = var(intra_ind_median_pers),
           percent_nas = median(intra_ind_percent_nas)) %>%
  data.frame() 

write.csv(ind_pop_exposure, "correct_ind_pers_by_population.csv",
          row.names = F) 
Species_exposure_score <- mean(ind_pop_exposure$pers)
Species_exposure_score 

# Creating a new (final) analysis dataframe
df_1 <- ind_pop_exposure
names(df_1)[1] <- "Colony"
df_2 <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/input_3h_csv/overlap_scores.csv")
df_2$Colony <- gsub(" ",".", df_2$Colony)
df_3 <- merge(df_1, df_2, by = "Colony")
analysis_df_for_hyp_3 <- df_3[,-6]

second_analysis_df_for_hyp_3 <- analysis_df_for_hyp_3[-c(1,2),]
cor.test(second_analysis_df_for_hyp_3$variance_pers, 
          second_analysis_df_for_hyp_3$Overlap_score,
          method = "kendall")
View(second_analysis_df_for_hyp_3)

# nas[c("population", "individual")] <- do.call(rbind, strsplit(as.character(nas$name), "_", fixed = T))
# nas_2 <- nas[,-c(1)]
# end_analysis_df <- merge(df_2, nas_2, by = c("population", "individual"))
# write.csv(end_analysis_df, "Analysis_dataframe_hyp_1.csv")
