# Bar plot for PERS 

# Loading essential packages----
library(tidyverse)
library(ggplot2)
library(viridisLite)
library(viridis)

# Loading data----
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/latest_right_attempt/outputs/csv/")
hyp1_df <- read.csv("correct_exposure_scores_by_individual_according_to_tracking_yr.csv")
hyp1_df$population <- gsub("\\.", " ", hyp1_df$population)
colonies_to_combine <- "Combined"
hyp1_df$population[hyp1_df$population %in% colonies_to_combine] <- "Iceland"
hyp1_df$population <- as.factor(hyp1_df$population)
order_to_match <- c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Iceland", "Faroe Islands", "Jarsteinen", "Eynhallow", "Isle of Canna", "Inishkea", "Little Saltee")
reversed_order_to_match <- rev(order_to_match)
hyp1_df$population <- factor(hyp1_df$population, levels = reversed_order_to_match)

# Calculate 95% confidence intervals for pers

df <- hyp1_df %>%
  group_by(population) %>%
  mutate(
    CI_lower = quantile(exposure_score, 0.025),
    CI_upper = quantile(exposure_score, 0.975)
  )
write.csv(df, "Hyp1_df_with_CI.csv")

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
CI_boxplot
ggsave("Boxplot_Colony_vs_PERS.png", plot = CI_boxplot, height = 8, width = 12, dpi = 900, unit = "in", bg = "white")

# Friendly colour scheme

CI_boxplot <- ggplot(data = df, aes(y = population, x = exposure_score)) +
  geom_boxplot(fill = viridis_pal(option = "D")(10), color = "black", alpha = 0.7) + # Use viridis palette for fill color
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
CI_boxplot
ggsave("Boxplot_Colony_vs_PERS.png", plot = CI_boxplot, height = 8, width = 12, dpi = 900, unit = "in", bg = "white")
summary(posthoc)
# Boxplot for Colony vs PERS but with significance codes----

# Obtain the summary of the final_model
model_summary <- summary(final_model)
model_summary

# Extract coefficient estimates, standard errors, and degrees of freedom
coefficients <- model_summary$coefficients[, "Estimate"]
standard_errors <- model_summary$coefficients[, "Std. Error"]

# Calculate z-values
z_values <- coefficients / standard_errors

# Calculate two-tailed p-values
p_values <- 2 * (1 - pnorm(abs(z_values)))

data.frame(Colony = rownames(model_summary$coefficients), Estimate = coefficients, 
           Std.Error = standard_errors, z_value = z_values, p_value = p_values)


# Create the data frame
summary_data <- data.frame(
  Colony = rownames(model_summary$coefficients),
  Estimate = coefficients,
  Std.Error = standard_errors,
  z_value = z_values,
  p_value = p_values
)
summary_data
write.csv(file = "Significance_codes_for_colonies.csv", summary_data, row.names = F)

# Reverse the order of factor levels in the population variable
df$population <- factor(df$population, levels = rev(levels(df$population)))

# Plot with reversed x-axis labels
CI_boxplot <- ggplot(data = df, aes(x = population, y = exposure_score)) +
  geom_boxplot(fill = rev(viridis_pal(option = "D")(10)), color = "black", alpha = 0.7) + # Use viridis palette for fill color
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.1,
    position = position_dodge(width = 0.75),
    color = "black" # Adjust color of error bars
  ) +
  xlab("Colony") + ylab("PERS") + # Corrected axis labels
  ggtitle("Colony vs PERS") +
  theme_minimal() + # Remove background grid
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center title and make it bold
    axis.line = element_line(color = "black"), # Adjust color of axis lines
    axis.title = element_text(face = "bold") # Make axis titles bold
  ) 

CI_boxplot

# Remove the intercept from the summary data
summary_data <- summary_data[-1, ]

# Extract colony names by removing the word "Colony" from the first column
summary_data$Colony <- sub("^Colony", "", summary_data$Colony)

# Define a function to convert p-values to stars
significance_stars <- function(p_values) {
  stars <- ifelse(p_values < 0.001, "***",
                  ifelse(p_values < 0.01, "**",
                         ifelse(p_values < 0.05, "*",
                                ifelse(p_values < 0.1, ".", ""))))
  return(stars)
}

citation("DHARMa")

# Apply the significance_stars function to the p-values
summary_data$stars <- significance_stars(summary_data$p_value)

# Plot with significance level on top of each box
CI_boxplot <- ggplot(data = df, aes(x = population, y = exposure_score)) +
  geom_boxplot(fill = rev(viridis_pal(option = "D")(10)), color = "black", alpha = 0.7) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.1,
    position = position_dodge(width = 0.75),
    color = "black"
  ) +
  geom_text(data = summary_data, aes(label = stars, x = Colony, y = max(df$exposure_score) * 1.05),
            size = 4, vjust = -0.5) +  # Adjust position and size of text
  xlab("Colony") + ylab("PERS") +
  ggtitle("Colony vs PERS") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black"),
    axis.title = element_text(face = "bold")
  )

CI_boxplot
ggsave("Boxplot_Colony_vs_PERS_with_significance.png", plot = CI_boxplot, height = 8, width = 12, dpi = 900, unit = "in", bg = "white")

# Removing Alkefjellet from the dataset----
df_wo_Alk <- df[!(df$population == "Alkefjellet"),]
# Removing Eynhallow outlier
df_wo_Eyn_outlier <- df_wo_Alk[!(df_wo_Alk$exposure_score < 10),]
# Reverse the color codes
color_codes <- rev(c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF",
                     "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF"))

CI_boxplot <- ggplot(data = df_wo_Eyn_outlier, aes(x = population, y = exposure_score)) +
  geom_boxplot(fill = color_codes, color = "black", alpha = 0.7) + # Use reversed color codes here
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.1,
    position = position_dodge(width = 0.75),
    color = "black"
  ) +
  geom_text(data = summary_data, aes(label = stars, x = Colony, y = max(df$exposure_score) * 1.05),
            size = 4, vjust = -0.5) +  # Adjust position and size of text
  xlab("Colony") + ylab("PERS") +
  ggtitle("Colony vs PERS") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black"),
    axis.title = element_text(face = "bold")
  )

CI_boxplot




# Doing a similar boxplot but with PAME regions----

# Loading data----
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/")
hyp2_df <- read.csv("correct_exposure_scores_by_month.csv")
hyp2_df$population <- gsub("\\.", " ", hyp2_df$population)
hyp2_df$population[hyp2_df$population == "Canadian Eastern Arctic   West Greenland"] <- "Canadian Eastern Arctic - West Greenland"
relevant_regions <- c("North Sea", "Faroe Plateau", "Iceland Shelf and Sea", "Canadian Eastern Arctic - West Greenland", "Barents Sea")
hyp2_df <- hyp2_df %>% filter(population %in% relevant_regions,)
hyp2_df$population <- as.factor(hyp2_df$population)
order_to_match <- c("Barents Sea", "Canadian Eastern Arctic - West Greenland", "Iceland Shelf and Sea", "Faroe Plateau", "North Sea")
reversed_order_to_match <- rev(order_to_match)
hyp2_df$population <- factor(hyp2_df$population, levels = reversed_order_to_match)

df <- hyp2_df %>%
  group_by(population) %>%
  mutate(
    CI_lower = quantile(exposure_score, 0.025),
    CI_upper = quantile(exposure_score, 0.975)
  )
write.csv(df, "Hyp2_df_with_CI.csv")

CI_boxplot <- ggplot(data = df, aes(y = population, x = exposure_score)) +
  geom_boxplot(fill = "lightblue", color = "black", alpha = 0.7) + # Adjust fill and outline colors
  geom_errorbar(
    aes(xmin = CI_lower, xmax = CI_upper),
    width = 0.1,
    position = position_dodge(width = 0.75),
    color = "black" # Adjust color of error bars
  ) +
  ylab("PAME region") + xlab("PERS") +
  ggtitle("PAME region vs PERS") +
  theme_minimal() + # Remove background grid
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center title and make it bold
    axis.line = element_line(color = "black"), # Adjust color of axis lines
    axis.title = element_text(face = "bold") # Make axis titles bold
  )
CI_boxplot
ggsave("Boxplot_Colony_vs_PERS.png", plot = CI_boxplot, height = 8, width = 12, dpi = 900, unit = "in", bg = "white")

# Friendly colour filling

CI_boxplot <- ggplot(data = df, aes(x = population, y = exposure_score)) +
  geom_boxplot(fill = viridis_pal(option = "D")(5), color = "black", alpha = 0.7) + # Use viridis palette for fill color
  geom_errorbar(
    aes(xmin = CI_lower, xmax = CI_upper),
    width = 0.1,
    position = position_dodge(width = 0.75),
    color = "black" # Adjust color of error bars
  ) +
  ylab("PAME region") + xlab("PERS") +
  ggtitle("PAME region vs PERS") +
  theme_minimal() + # Remove background grid
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center title and make it bold
    axis.line = element_line(color = "black"), # Adjust color of axis lines
    axis.title = element_text(face = "bold") # Make axis titles bold
  )
CI_boxplot
ggsave("Boxplot_PAME_region_vs_PERS.png", plot = CI_boxplot, height = 8, width = 12, dpi = 900, unit = "in", bg = "white")

# Suggestions

CI_boxplot <- ggplot(data = df, aes(x = population, y = exposure_score)) +
  geom_boxplot(fill = viridis_pal(option = "D")(5), color = "black", alpha = 0.7, position = position_dodge(width = 0.75)) + # Use viridis palette for fill color
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.1,
    position = position_dodge(width = 0.75),
    color = "black" # Adjust color of error bars
  ) +
  xlab("PAME region") + ylab("PERS") +
  ggtitle("PERS vs PAME region") +
  theme_minimal() + # Remove background grid
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center title and make it bold
    axis.line = element_line(color = "black"), # Adjust color of axis lines
    axis.title = element_text(face = "bold") # Make axis titles bold
  )
CI_boxplot
ggsave("Boxplot_PERS_vs_PAME_region.png", plot = CI_boxplot, height = 8, width = 12, dpi = 900, unit = "in", bg = "white")
