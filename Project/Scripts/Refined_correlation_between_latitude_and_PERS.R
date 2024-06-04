# Refined correlation test between increasing latitudes and PERS

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
nbs_mylocs$colony[nbs_mylocs$colony %in% colonies_to_combine] <- "Iceland"

desc_nbs_mylocs <- nbs_mylocs[order(-nbs_mylocs$col_lat),]
order_to_match <- unique(desc_nbs_mylocs$colony) 
order_to_match
lat_df <- nbs_mylocs %>% dplyr::select(c("colony", "col_lat")) %>% distinct()


df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/latest_right_attempt/outputs/csv/correct_ind_pers_by_population.csv")
df$population[df$population == "Combined"] <- "Iceland"
df$population <- gsub("\\.", " ", df$population)
df$population <- factor(df$population, levels = order_to_match)
names(df)[1] <- "colony"

final_df <- merge(df, lat_df, by = "colony")
View(final_df)

correlation <- cor.test(final_df$pers, final_df$col_lat)
print(correlation)

final_df_sorted <- final_df %>% arrange(desc(col_lat))
View(final_df_sorted)

final_df_wo <- final_df_sorted[-2,]
correlation_wo <- cor.test(final_df_wo$pers, final_df_wo$col_lat)
print(correlation_wo)

# Unrelated: Calculating 95% CI for %NA----
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/latest_right_attempt/outputs/csv/")
nas_df <- read.csv("correct_nas_ind.csv")
nas_df[c("Colony", "Individ_id", "Tracking_year")] <- do.call(rbind, str_split(as.character(nas_df$name), pattern = "_"))
View(nas_df)
df <- nas_df %>%
  group_by(Colony) %>%
  mutate(
    CI_lower = quantile(percent_na, 0.025),
    CI_upper = quantile(percent_na, 0.975)
  )
View(df)
write.csv(df, "Hyp1_NA_df_with_CI.csv")
