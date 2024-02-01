# Home range overlap script using the amt package

# Loading essential packages----
library(amt)
library(dplyr)
library(ggplot2)
library(sp)

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
  crs = 4326) %>%
  nest(data = -"id") %>% 
  arrange(id)
attr(NF_tracks, "crs") <- DgProj

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
                                    type = "vi",
                                    conditional = F)

output_df <- NF_tracks_kde_overlap %>% group_by(from) %>%
  summarize(mean_overlap = mean(overlap))

colony_hr_kde_overlap_index <- median(output_df$mean_overlap) 
to_bind <- data.frame("Colony" = i, "Overlap_score" = colony_hr_kde_overlap_index)
results_df <- rbind(results_df, to_bind)

} # First for loop ends

kde_final_overlap <- results_df[-(1:10),]
kde_final_overlap$Colony <- gsub(" ",".",kde_final_overlap$Colony)

variance_df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/exposure_scores_by_individual.csv")
to_combine <- c("Skjalfandi", "Langanes")
variance_df$population[variance_df$population %in% to_combine] <- "Combined"

variance_output <- variance_df %>% group_by(population) %>% summarize(variance = var(exposure_score)) 
colnames(variance_output)[colnames(variance_output) == "population"] <- "Colony"

# correlation test
corr_df <- merge(kde_final_overlap, variance_output, by = "Colony")
View(corr_df)

correlation <- cor.test(corr_df$Overlap_score, corr_df$variance, method = "kendall")
print(correlation)
