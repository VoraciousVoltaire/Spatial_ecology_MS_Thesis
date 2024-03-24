# Code for jitter plot

# Jitter plot----

# Loading data: 

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/latest_right_attempt/outputs/csv")
df <- read.csv("correct_exposure_scores_by_individual_according_to_tracking_yr.csv")
corrected_df <- df[df$exposure_score > 10,]
analysis_df <- corrected_df[!(corrected_df$population == "Jan.Mayen" & 
                                corrected_df$individual == "NOS.4181597" & 
                                corrected_df$tracking_year == 2018), ]
analysis_df <- analysis_df[,-c(4,5)]
colnames(analysis_df) <- c("Colony", "Individ_id", "Tracking_year", "pers")
analysis_df$Colony[analysis_df$Colony == "Combined"] <- "Iceland"
analysis_df$Colony <- gsub("\\.", " ", analysis_df$Colony)
analysis_df$Colony <- as.factor(analysis_df$Colony)
order_to_match <- c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Iceland", "Faroe Islands", "Jarsteinen", "Eynhallow", "Isle of Canna", "Inishkea", "Little Saltee")
reversed_order_to_match <- rev(order_to_match)
analysis_df$Colony <- factor(analysis_df$Colony, levels = reversed_order_to_match)

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/")
color_vector <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF",
                  "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")

jitter_plot <- ggplot(analysis_df, aes(x = pers, y = Colony, col = Colony)) +
  geom_jitter(width = 0, height = 0.2) +  # Add jitter for better visualization
  scale_color_manual(values = color_vector) +  # Specify colors manually
  labs(x = "Plastic Exposure Risk Score (PERS)", y = "Colony") +
  theme_minimal() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  ggtitle("Tracking-year wise individual-specific PERS across Northern fulmar colonies") +
  guides(col = guide_legend(reverse = TRUE))  # Reverse order of legend

jitter_plot

ggsave("Jitter_plot_for_PERS.png", height = 8, width = 14, unit = "in", dpi = 900, plot = jitter_plot, bg = "white")
