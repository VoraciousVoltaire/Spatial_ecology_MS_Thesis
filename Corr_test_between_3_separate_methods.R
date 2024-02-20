# Checking whether pers values positively correlate with each other by three different methods: by ind, by month, by tracking year
# by remembering the fact that 

df_month <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Month/outputs/csv/month_exposure_scores_by_population.csv")
View(df_month)

df_tracking_year <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/ty_exposure_scores_by_population.csv")
View(df_tracking_year)

df_ind <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/ind_exposure_scores_by_population.csv")
View(df_ind)

Species_exposure_score_ind <- mean(df_ind$population_exposure)
Species_exposure_score_ind

Species_exposure_score_month <- mean(df_month$population_exposure)
Species_exposure_score_month

Species_exposure_score_ty <- mean(df_tracking_year$population_exposure)
Species_exposure_score_ty

colnames(df_ind) <- c("Colony", "ind_pers")
colnames(df_month) <- c("Colony", "month_pers")
colnames(df_tracking_year) <- c("Colony", "ty_pers")

three_dfs_merged <- merge(merge(df_ind, df_month, by = "Colony"), df_tracking_year, by = "Colony")

# Correlation test----

# Normality test first

# Function to perform Shapiro-Wilk test on each column
shapiro_multiple <- function(df) {
  results <- lapply(df, function(col) shapiro.test(col))
  return(results)
}

# Perform Shapiro-Wilk test on each column
shapiro_results <- shapiro_multiple(three_dfs_merged[,-1])

print(shapiro_results) # p = 0.9875 (ind), 0.1648 (month), 0.2803 (tracking year)

three_dfs_cor_test_1 <- cor(three_dfs_merged[,-1], method = "kendall")
print(three_dfs_cor_test_1)

sep_cor_test_1 <- cor.test(three_dfs_merged$ind_pers, three_dfs_merged$month_pers, method = "kendall")
sep_cor_test_2 <- cor.test(three_dfs_merged$ind_pers, three_dfs_merged$ty_pers, method = "kendall")
sep_cor_test_3 <- cor.test(three_dfs_merged$month_pers, three_dfs_merged$ty_pers, method = "kendall")

print(sep_cor_test_1)
print(sep_cor_test_2)
print(sep_cor_test_3)

