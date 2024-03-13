# Exploring hyp 2 dataset----
getwd()

df_hyp2 <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/final_df_for_hyp_2.csv")

regions_visited <- data.frame(matrix(ncol = 2))
colnames(regions_visited) <- c("Colony", "Regions")

for(i in unique(df_hyp2$colony)){ # first for loop starts
  sub <- df_hyp2[df_hyp2$colony == i,]
  unique(sub$LME_NAME)
  print(paste0("Individuals from colony ",i," visit the following ",length(unique(sub$LME_NAME))," regions: ",unique(sub$LME_NAME)))
  to_add <- data.frame(Colony = i, Regions = unique(sub$LME_NAME))
  regions_visited <- rbind(regions_visited, to_add)
} # first for loop ends

View(regions_visited_2)
regions_visited_2 <- na.omit(regions_visited)
regions_visited_2 %>% group_by(Colony) %>% summarise(n = n())
length(unique(df_hyp2$LME_NAME))
