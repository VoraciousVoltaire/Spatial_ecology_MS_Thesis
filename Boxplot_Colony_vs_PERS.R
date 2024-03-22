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

# Doing a similar boxplot but with PAME regions 

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

CI_boxplot <- ggplot(data = df, aes(y = population, x = exposure_score)) +
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
