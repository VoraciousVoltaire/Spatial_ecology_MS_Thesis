# Bar plot for PERS 

hyp1_df <- read.csv("correct_exposure_scores_by_individual_according_to_tracking_yr.csv")
hyp1_df$population <- gsub("\\.", " ", hyp1_df$population)
colonies_to_combine <- "Combined"
hyp1_df$population[hyp1_df$population %in% colonies_to_combine] <- "Iceland"

# Calculate 95% confidence intervals for pers

df <- hyp1_df %>%
  group_by(population) %>%
  mutate(
    CI_lower = quantile(exposure_score, 0.025),
    CI_upper = quantile(exposure_score, 0.975)
  )
write.csv(df, "Hyp1_df_with_CI.csv")

# Create the boxplot with error bars
library(ggplot2)

# Assuming df is your data frame

CI_boxplot <- ggplot(data = df, aes(y = population, x = exposure_score)) +
  geom_boxplot(fill = "lightblue", color = "black", alpha = 0.7) + # Adjust fill and outline colors
  geom_errorbar(
    aes(xmin = CI_lower, xmax = CI_upper),
    width = 0.1,
    position = position_dodge(width = 0.75),
    color = "black" # Adjust color of error bars
  ) +
  ylab("Colony") + xlab("PERS") +
  ggtitle("Colony vs PERS") +
  theme_minimal() + # Remove background grid
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center title and make it bold
    axis.line = element_line(color = "black"), # Adjust color of axis lines
    axis.title = element_text(face = "bold") # Make axis titles bold
  )

ggsave("Boxplot_Colony_vs_PERS.png", plot = CI_boxplot, height = 8, width = 12, dpi = 900, unit = "in", bg = "white")

