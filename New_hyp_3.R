# Hypothesis 3 but kernels calculated on a sampling year basis

# first testing the correlation between the values calculated on the basis of amt and adehabitathr kde overlap functions----

getwd()
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
nbs_years <- sort(unique(df_mod$tracking_year))

# Defining land----

land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# Now calculating kernel hr overlap on a tracking_year basis----

results_df <- data.frame(matrix(ncol = 2))
colnames(results_df) <- c("Colony", "Overlap_score")

for(i in unique(df_mod$colony)){ # First for loop begins
  sub <- df_mod %>% filter(colony == i)
 
  for(j in unique(sub$tracking_year)){ # Second for loop begins
    tracks_nbs <- sub[sub$tracking_year == j,] 
    
    # Setting a colony-centered crs 
    median_loc <- cbind(median(tracks_nbs$lon), median(tracks_nbs$lat))
    tracks_nbs
    DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",median_loc[1]," +lat_0=",median_loc[2]," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m ", sep="")) 
    
    # Creating a colony-specific spatial points data frame
    sp::coordinates(tracks_nbs) <- ~ lon + lat
    sp::proj4string(tracks_nbs) <- sp::proj4string(land)
    
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
      
      tracks <- sp::spTransform(tracks_nbs, CRS = DgProj)
      tracks$individ_id <- factor(tracks$individ_id)
      
      # Kernel estimation
      kudl <- kernelUD(tracks[,"individ_id"], grid = null.grid, h = 200000) ## smoothing factor equals 200 km for GLS data
      length(attributes(kudl)$names)
      
      # Second if loop to eliminate datasets with just 1 animal
      
      if(length(attributes(kudl)$names) > 1){ # Second if loop starts
        
        # Calculating home range overlap using kerneloverlaphr
        ade_kde_overlap <- kerneloverlaphr(kudl, meth = "BA", percent = 95, conditional = F)
        range(as.data.frame(ade_kde_overlap)) # 0 to 0.9998733
        df <- as.data.frame(ade_kde_overlap)
        median_df <- apply(df, 1, median, na.rm = T)
        median_overlap_df <- as.data.frame(median_df)
        colnames(median_overlap_df) <- "overlap"
        colony_hroverlap <- median(median_overlap_df$overlap) 
        to_bind <- data.frame("Colony" = paste0(i,"_",j), "Overlap_score" = colony_hroverlap)
        results_df <- rbind(results_df, to_bind)
        
        
      } # Second if loop ends
      
    } # First if loop ends
  } # Second for loop ends
  } # First for loop ends
  
  View(results_df)
  kde_final_overlap <- results_df[-1,]
  View(kde_final_overlap)
  write.csv(kde_final_overlap, "/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/new_third_hyp.csv")
  
  # String spliting
  kde_final_overlap[c("Just_colony", "Tracking_year")] <-  do.call(rbind, strsplit(as.character(kde_final_overlap$Colony), "_"))  
  View(kde_final_overlap)
  
  # Calculating a median accoridng to just colony
  new_df <- kde_final_overlap %>% group_by(Just_colony) %>% summarize(median_overlap_score <- median(Overlap_score))
  colnames(new_df) <- c("Colony", "Overlap_score")
  new_df$Colony <- gsub(" ", ".", new_df$Colony)
  
 
  
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
  
  