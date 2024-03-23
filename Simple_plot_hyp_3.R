# Image for simple plot between PERS variance and overlap index

# Loading essential packages----
library(psych)

# Loading data----
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Third_hyp/output/csv/")
hyp3_df <- read.csv("final_new_third_hyp_3_with_ba.csv")
hyp3_df <- hyp3_df[,-1]
hyp3_df <- hyp3_df %>% rowwise() %>% mutate(harmonic_mean = harmonic.mean(c(udoi_Overlap_score, ba_Overlap_score)))
hyp3_df[,c("Colony", "Tracking_year")] <- do.call(rbind, str_split(hyp3_df$Colony, pattern = "_"))
df <- hyp3_df %>%
  group_by(Colony) %>%
  mutate(
    CI_lower = quantile(harmonic_mean, 0.025),
    CI_upper = quantile(harmonic_mean, 0.975)
  )
write.csv(df, "Hyp3_df_with_CI.csv")

df_var <- read.csv("Variance_with_median_percent_na.csv")
df_var <- df_var[,-1]
df_var$Colony <- gsub("\\.", " ", df_var$Colony)

comprehensive_df <- merge(df, df_var, by = "Colony")
write.csv(df, "Comprehensive_hyp3_df.csv")

df <- comprehensive_df %>% group_by(Colony) %>% summarise(Overlap_score = harmonic.mean(harmonic_mean))
df_final <- merge(comprehensive_df, df, by = "Colony")
df_final$Colony[df_final$Colony == "Combined"] <- "Iceland"
df_final$Colony <- as.factor(df_final$Colony)
order_to_match <- c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Iceland", "Faroe Islands", "Jarsteinen", "Eynhallow", "Isle of Canna", "Inishkea", "Little Saltee")
reversed_order_to_match <- rev(order_to_match)
reversed_order_to_match
df_final$Colony <- factor(df_final$Colony, levels = reversed_order_to_match)
unique(df_final$Colony)

simple_plot <- ggplot(data = df_final, aes(y = variance_pers, x = Overlap_score)) +
  geom_point() +  # Add points
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.1) +  # Add horizontal error bars
  ylab("Variance in PERS") + xlab("Home range overlap score") +  # Label axes
  ggtitle("Variance in PERS vs Home range overlap score") +  # Title
  theme_minimal() +  # Minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title and make it bold
    axis.line = element_line(color = "black"),  # Adjust color of axis lines
    axis.title = element_text(face = "bold")  # Make axis titles bold
  )
simple_plot

# Better readability----
# Create the plot with error bars representing confidence intervals
simple_plot <- ggplot(data = df_final, aes(y = variance_pers, x = Overlap_score)) +
  geom_point(size = 3) +  # Increase point size for better visibility
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.1) +  # Add horizontal error bars
  ylab("Variance in PERS") + xlab("Home range overlap score") +  # Label axes
  ggtitle("Variance in PERS vs Home range overlap score") +  # Title
  ylim(min(df_final$variance_pers) - 0.1*abs(min(df_final$variance_pers)), max(df_final$variance_pers) + 0.1*abs(max(df_final$variance_pers))) +  # Expand Y axis
  theme_minimal() +  # Minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title and make it bold
    axis.line = element_line(color = "black"),  # Adjust color of axis lines
    axis.title = element_text(face = "bold")  # Make axis titles bold
  )

simple_plot


# Create the plot with error bars representing confidence intervals
simple_plot <- ggplot(data = df_final, aes(y = variance_pers, x = Overlap_score, fill = Colony)) +
  geom_point(size = 2) +  # Increase point size for better visibility
  scale_fill_manual(values = viridis_pal(option = "D")(10)) + 
  ylab("Variance in PERS") + xlab("Home range overlap score") +  # Label axes
  ggtitle("Variance in PERS vs Home range overlap score") +  # Title
  theme_minimal() +  # Minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title and make it bold
    axis.line = element_line(color = "black"),  # Adjust color of axis lines
    axis.title = element_text(face = "bold")  # Make axis titles bold
  )
simple_plot

# Plot with the reversed color scheme
simple_plot <- ggplot(data = df_final, aes(y = variance_pers, x = Overlap_score, fill = Colony)) +
  geom_point(size = 6, shape = 21) +  # Increase point size for better visibility, use shape 21 for filled circles
  scale_fill_manual(values = viridis_pal(option = "D")(10)) + 
  ylab("Variance in PERS") + xlab("Home range overlap score") +  # Label axes
  ggtitle("Variance in PERS vs Home range overlap score") +  # Title
  theme_minimal() +  # Minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title and make it bold
    axis.line = element_line(color = "black"),  # Adjust color of axis lines
    axis.title = element_text(face = "bold")  # Make axis titles bold
  )
simple_plot
ggsave(plot = simple_plot, "Variance_in_PERS_vs_Home_range_overlap_score.png", bg = "white", height = 8, width = 14, dpi = 900)

# Colour palette code:

# Define the colony list
colonies <- c("Little Saltee", "Inishkea", "Isle of Canna", "Eynhallow", "Jarsteinen",
              "Faroe Islands", "Iceland", "Jan Mayen", "Bjørnøya", "Alkefjellet")

# Define the corresponding color codes (in reverse order)
color_codes <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF",
                 "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")

# Define the corresponding color names
color_names <- c("Dark Purple", "Purple", "Indigo", "Navy Blue", "Blue",
                 "Teal", "Turquoise Green", "Green", "Lime Green", "Yellow")

# Create a data frame
colony_colors_df <- data.frame(Colony = colonies, Color = color_codes, Color_Name = color_names)

# Write the data frame to a CSV file
write.csv(colony_colors_df, file = "colony_colors.csv", row.names = FALSE)





