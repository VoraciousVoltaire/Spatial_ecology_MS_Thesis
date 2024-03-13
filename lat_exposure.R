# Script for arranging colonies according to increasing latitudes 

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

desc_nbs_mylocs <- nbs_mylocs[order(-nbs_mylocs$col_lat),]
order_to_match <- unique(desc_nbs_mylocs$colony) 
order_to_match

df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/exposure_scores_by_individual.csv")
to_combine <- c("Skjalfandi", "Langanes")
df$population[df$population %in% to_combine] <- "Combined"

pop_exposure_by_lat <- df %>%
  group_by(population) %>%  
  summarise(population_exposure = round(mean(exposure_score), 4))

combined_lat <- nbs_mylocs
combined_lat[combined_lat$colony %in% to_combine] <- "Combined"
combined_lat$col_lat[combined_lat$colony == "Combined"] <- median(unique(combined_lat$col_lat[combined_lat$colony == "Combined"]))

lat <- sort(unique(combined_lat$col_lat), decreasing = T)
final_df <- cbind(pop_exposure_by_lat, lat)
final_df

# Checking normality of variables
shapiro.test(final_df$population_exposure)
shapiro.test(final_df$lat)

correlation <- cor.test(final_df$population_exposure, final_df$lat, method = "spearman")
print(correlation)


