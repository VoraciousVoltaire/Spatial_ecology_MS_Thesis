# For Supplementary info----

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/"
mylocs <- readRDS(paste0(datadir, "SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"summaryTable.rds"))

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
input_df <- indiv_merged_df %>% 
  dplyr::mutate(year = year(timestamp)) %>%
  dplyr::mutate(month = month(timestamp)) %>%
  mutate(tracking_year = ifelse(month < 7, year - 1, year))
names(input_df)[names(input_df) == "ring"] <- "individ_id"

# Exploring data----

# 2 dataframes for tracking year
csv_1 <- input_df %>% 
dplyr::group_by(colony,tracking_year) %>%
dplyr::summarise(ty_number_birds_tracked = n_distinct(individ_id))

csv_2 <- input_df %>% 
dplyr::group_by(colony,tracking_year) %>%
dplyr::summarise(ty_number_of_total_tracks = n()) 

# 2 dataframes for month
csv_3 <- input_df %>% 
  dplyr::group_by(colony,month) %>%
  dplyr::summarise(month_number_birds_tracked = n_distinct(individ_id))

csv_4 <- input_df %>% 
  dplyr::group_by(colony,month) %>%
  dplyr::summarise(month_number_of_total_tracks = n()) 


# 1 dataframe for individuals
csv_5 <- input_df %>% 
  dplyr::group_by(colony,individ_id) %>%
  dplyr::summarise(ind_number_of_total_tracks = n())

dim(csv_5[!(csv_5$ind_number_of_total_tracks < 50 | csv_5$ind_number_of_total_tracks == 50),]) == dim(csv_5) # All sample sizes for hre > 50 
range(csv_5$ind_number_of_total_tracks) #  101 6103

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Supplementary_info/csv/")

output_1 <- merge(csv_1, csv_2, by = c("colony", "tracking_year"))
dim(output_1[!(output_1$ty_number_of_total_tracks < 50 | output_1$ty_number_of_total_tracks == 50),]) == dim(output_1) # All sample sizes for hre > 50 
range(output_1$ty_number_of_total_tracks) #  73 31286
write.csv(output_1, "tracking_year_number_of_tracks.csv")

output_2 <- merge(csv_3, csv_4, by = c("colony", "month"))
dim(output_2[!(output_2$month_number_of_total_tracks < 50 | output_2$month_number_of_total_tracks == 50),]) == dim(output_2) # All sample sizes for hre > 50 
range(output_2$month_number_of_total_tracks) # 120 23180
write.csv(output_2, "month_number_of_tracks.csv")

write.csv(csv_5, "ind_number_of_tracks.csv")
