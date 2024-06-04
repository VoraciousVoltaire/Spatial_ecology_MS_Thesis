# Doing a simple dot plot for Hypothesis 2

# Loading data----
hyp2_analysis_df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/correct_hyp_2_analysis_df.csv")
hyp2_analysis_df$Region <- gsub("\\.", " ", hyp2_analysis_df$Region)


simple_plot <- ggplot(data = hyp2_analysis_df, aes(y = EcoQO.value, x = median_pers , fill = Region)) +
  geom_point(size = 6, shape = 21) +  # Increase point size for better visibility, use shape 21 for filled circles
  scale_fill_manual(values = rev(viridis_pal(option = "D")(5))) + 
  ylab("EcoQO") + xlab("PERS") +  # Label axes
  ggtitle("EcoQO vs PERS") +  # Title
  theme_minimal() +  # Minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title and make it bold
    axis.line = element_line(color = "black"),  # Adjust color of axis lines
    axis.title = element_text(face = "bold")  # Make axis titles bold
  )
simple_plot

# Adding a straight line----

simple_plot <- ggplot(data = hyp2_analysis_df, aes(y = EcoQO.value, x = median_pers , fill = Region)) +
  geom_point(size = 6, shape = 21) +  # Increase point size for better visibility, use shape 21 for filled circles
  geom_abline(intercept = coef(lm(EcoQO.value ~ median_pers, data = hyp2_analysis_df))[1], 
              slope = coef(lm(EcoQO.value ~ median_pers, data = hyp2_analysis_df))[2], 
              color = "black") +  # Add linear regression line
  scale_fill_manual(values = rev(viridis_pal(option = "D")(5))) + 
  ylab("EcoQO") + xlab("PERS") +  # Label axes
  ggtitle("EcoQO vs PERS") +  # Title
  theme_minimal() +  # Minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title and make it bold
    axis.line = element_line(color = "black"),  # Adjust color of axis lines
    axis.title = element_text(face = "bold")  # Make axis titles bold
  ) +
  labs(fill = "PAME Region")
simple_plot
ggsave(plot = simple_plot, "EcoQO_vs_PERS_v2.png", bg = "white", height = 14, width = 10, dpi = 900)


