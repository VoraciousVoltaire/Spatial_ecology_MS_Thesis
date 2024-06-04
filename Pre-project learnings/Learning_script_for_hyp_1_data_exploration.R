# Learning script for data exploration----

# Loading essential packages----
library(ggplot2)

exp_df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/csv/Final_hypothesis_1_analysis_df.csv")
p1 <- ggplot(data = exp_df, aes(x = Colony, y = pers)) + geom_boxplot() + theme_minimal()
p1
p2 <- ggplot(exp_df, aes(x = pers, y = percent_na)) + 
geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, binwidth = 1)
summary(exp_df)
p2

# Checking for outliers----
# Set up the plotting area
par(mfrow = c(2, 3))

# Create dotplots for each column
for (col in c("pers", "Number_of_tracks", "vals", "nas", "percent_na")) {
  dotchart(exp_df[[col]], main = col, xlab = col, pch = 19)
}
 
# # ggplot----
# library(ggplot2)
# 
# # Create a long-form dataframe for ggplot
# long_df <- reshape2::melt(exp_df[, c("pers", "Number_of_tracks", "vals", "nas", "percent_na")])
# 
# # Create multipanel dotplots using ggplot
# ggplot(long_df, aes(x = variable, y = value)) +
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
#   facet_wrap(~variable, scales = "free") +
#   labs(title = "Multipanel Cleveland Dotplots")

# Checking for homogenity of variance----
library(ggplot2)

# Create a scatterplot to check homogeneity of variance
ggplot(exp_df, aes(x = Colony, y = pers)) +
  geom_point() +
  labs(title = "Scatterplot of pers by Colony",
       x = "Colony",
       y = "pers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

library(car)

# Perform Levene's test to assess homogeneity of variance
leveneTest(pers ~ Colony, data = exp_df)

# Load the car package for the oneway.test function
library(car)

# Perform Welch's ANOVA
welch_anova <- oneway.test(pers ~ Colony, data = exp_df)

# Print the results
print(welch_anova)

# Post-hoc test----
# Install and load the emmeans package
install.packages("emmeans")
library(emmeans)

# Perform Welch's ANOVA
welch_anova <- aov(pers ~ Colony, data = exp_df)
welch_anova_emm <- emmeans(welch_anova, "Colony")

# Perform pairwise comparisons with Bonferroni correction
pairwise_comp <- pairs(welch_anova_emm, adjust = "bonferroni")

# Print the results
print(pairwise_comp)

# Result----
| Comparison                        | Estimate | SE    | df  | t-ratio | p-value | Significant |
  |----------------------------------|----------|-------|-----|---------|---------|-------------|
  | Alkefjellet - Bjørnøya           | -3.1642  | 0.390 | 378 | -8.115  | < 0.0001| Yes         |
  | Alkefjellet - Combined           | -1.6602  | 0.367 | 378 | -4.519  | 0.0004  | Yes         |
  | Alkefjellet - Eynhallow          | -4.5958  | 0.348 | 378 | -13.210 | < 0.0001| Yes         |
  | Alkefjellet - Faroe Islands      | -3.5187  | 0.428 | 378 | -8.218  | < 0.0001| Yes         |
  | Alkefjellet - Inishkea           | -4.8375  | 0.720 | 378 | -6.723  | < 0.0001| Yes         |
  | Alkefjellet - Isle of Canna      | -4.1699  | 0.560 | 378 | -7.449  | < 0.0001| Yes         |
  | Alkefjellet - Jan Mayen          | -3.2224  | 0.372 | 378 | -8.651  | < 0.0001| Yes         |
  | Alkefjellet - Jarsteinen         | -4.9554  | 0.467 | 378 | -10.613 | < 0.0001| Yes         |
  | Alkefjellet - Little Saltee      | -6.4048  | 0.400 | 378 | -16.006 | < 0.0001| Yes         |
  | Bjørnøya - Combined              | 1.5040   | 0.263 | 378 | 5.726   | < 0.0001| Yes         |
  | Bjørnøya - Eynhallow             | -1.4316  | 0.235 | 378 | -6.101  | < 0.0001| Yes         |
  | Bjørnøya - Faroe Islands         | -0.3545  | 0.343 | 378 | -1.035  | 1.0000  | No          |
  | Bjørnøya - Inishkea              | -1.6733  | 0.672 | 378 | -2.489  | 0.5952  | No          |
  | Bjørnøya - Isle of Canna         | -1.0057  | 0.497 | 378 | -2.022  | 1.0000  | No          |
  | Bjørnøya - Jan Mayen             | -0.0582  | 0.270 | 378 | -0.216  | 1.0000  | No          |
  | Bjørnøya - Jarsteinen            | -1.7913  | 0.390 | 378 | -4.594  | 0.0003  | Yes         |
  | Bjørnøya - Little Saltee         | -3.2406  | 0.307 | 378 | -10.562 | < 0.0001| Yes         |
  | Combined - Eynhallow             | -2.9356  | 0.195 | 378 | -15.064 | < 0.0001| Yes         |
  | Combined - Faroe Islands         | -1.8585  | 0.317 | 378 | -5.869  | < 0.0001| Yes         |
  | Combined - Inishkea              | -3.1773  | 0.659 | 378 | -4.819  | 0.0001  | Yes         |
  | Combined - Isle of Canna         | -2.5097  | 0.480 | 378 | -5.229  | < 0.0001| Yes         |
  | Combined - Jan Mayen             | -1.5622  | 0.236 | 378 | -6.620  | < 0.0001| Yes         |
  | Combined - Jarsteinen            | -3.2953  | 0.367 | 378 | -8.970  | < 0.0001| Yes         |
  | Combined - Little Saltee         | -4.7446  | 0.278 | 378 | -17.092 | < 0.0001| Yes         |
  | Eynhallow - Faroe Islands        | 1.0771   | 0.294 | 378 | 3.665   | 0.0127  | Yes         |
  | Eynhallow - Inishkea             | -0.2417  | 0.649 | 378 | -0.373  | 1.0000  | No          |
  | Eynhallow - Isle of Canna        | 0.4259   | 0.465 | 378 | 0.915   | 1.0000  | No          |
  | Eynhallow - Jan Mayen            | 1.3734   | 0.204 | 378 | 6.721   | < 0.0001| Yes         |
  | Eynhallow - Jarsteinen           | -0.3597  | 0.348 | 378 | -1.034  | 1.0000  | No          |
  | Eynhallow - Little Saltee        | -1.8090  | 0.251 | 378 | -7.200  | < 0.0001| Yes         |
  | Faroe Islands - Inishkea         | -1.3188  | 0.695 | 378 | -1.897  | 1.0000  | No          |
  | Faroe Islands - Isle of Canna    | -0.6512  | 0.528 | 378 | -1.233  | 1.0000  | No          |
  | Faroe Islands - Jan Mayen        | 0.2963   | 0.323 | 378 | 0.919   | 1.0000  | No          |
  | Faroe Islands - Jarsteinen       | -1.4367  | 0.428 | 378 | -3.355  | 0.0393  | Yes         |
  | Faroe Islands - Little Saltee    | -2.8861  | 0.354 | 378 | -8.149  | < 0.0001| Yes         |
  | Inishkea - Isle of Canna         | 0.6676   | 0.783 | 378 | 0.853   | 1.0000  | No          |
  | Inishkea - Jan Mayen             | 1.6151   | 0.662 | 378 | 2.439   | 0.6835  | No          |
  | Inishkea - Jarsteinen            | -0.1180  | 0.720 | 378 | -0.164  | 1.0000  | No          |
  | Inishkea - Little Saltee         | -1.5673  | 0.678 | 378 | -2.311  | 0.9613  | No          |
  | Isle of Canna - Jan Mayen        | 0.9475   | 0.484 | 378 | 1.958   | 1.0000  | No          |
  | Isle of Canna - Jarsteinen       | -0.7855  | 0.560 | 378 | -1.403  | 1.0000  | No          |
  | Isle of Canna - Little Saltee    | -2.2349  | 0.505 | 378 | -4.421  | 0.0006  | Yes         |
  | Jan Mayen - Jarsteinen           | -1.7331  | 0.372 | 378 | -4.653  | 0.0002  | Yes         |
  | Jan Mayen - Little Saltee        | -3.1824  | 0.284 | 378 | -11.194 | < 0.0001| Yes         |
  | Jarsteinen - Little Saltee       | -1.4494  | 0.400 | 378 | -3.622  | 0.0149  | Yes         |
  
  P-value adjustment: Bonferroni method for 45 tests

# Using robust linear methods to do the same thing----
# Load the MASS package for the rlm function
library(MASS)

# Perform robust regression
robust_lm <- rlm(pers ~ Colony, data = exp_df)

# Print the results
print(summary(robust_lm))

# Post-hoc test now----
# Load necessary libraries
library("WRS2")  # For robust post hoc tests
install.packages("WRS2")
# Load necessary libraries
library("WRS2")  # For robust post hoc tests

# Extract residuals from the robust linear regression model
residuals <- residuals(robust_lm)

# Conduct robust post hoc tests with Bonferroni correction
# Load necessary libraries
library("WRS2")  # For robust post hoc tests

# Conduct robust post hoc tests with Bonferroni correction
pairwise_comparisons <- pairwise.wilcox.test(exp_df$pers, exp_df$Colony, 
                                             p.adjust.method = "bonferroni",
                                             alternative = "two.sided")  # Specify two-sided test
print(pairwise_comparisons)

# Result
|          Comparison         | Estimate | Std. Error |  t value |  p value | Significant |
  |-----------------------------|----------|------------|----------|----------|-------------|
  | Alkefjellet - Bjørnøya      |  -3.1642 |     0.390  |   -8.115 |   <.0001 |     Yes     |
  | Alkefjellet - Combined      |  -1.6602 |     0.367  |   -4.519 |   0.0004 |     Yes     |
  | Alkefjellet - Eynhallow     |  -4.5958 |     0.348  |  -13.210 |   <.0001 |     Yes     |
  | Alkefjellet - Faroe.Islands |  -3.5187 |     0.428  |   -8.218 |   <.0001 |     Yes     |
  | Alkefjellet - Inishkea      |  -4.8375 |     0.720  |   -6.723 |   <.0001 |     Yes     |
  | Alkefjellet - Isle.of.Canna |  -4.1699 |     0.560  |   -7.449 |   <.0001 |     Yes     |
  | Alkefjellet - Jan.Mayen     |  -3.2224 |     0.372  |   -8.651 |   <.0001 |     Yes     |
  | Alkefjellet - Jarsteinen    |  -4.9554 |     0.467  |  -10.613 |   <.0001 |     Yes     |
  | Alkefjellet - Little.Saltee |  -6.4048 |     0.400  |  -16.006 |   <.0001 |     Yes     |
  | Bjørnøya - Combined         |   1.5040 |     0.263  |    5.726 |   <.0001 |     Yes     |
  | Bjørnøya - Eynhallow        |  -1.4316 |     0.235  |   -6.101 |   <.0001 |     Yes     |
  | Bjørnøya - Faroe.Islands    |  -0.3545 |     0.343  |   -1.035 |    1.000 |      No     |
  | Bjørnøya - Inishkea         |  -1.6733 |     0.672  |   -2.489 |    0.595 |      No     |
  | Bjørnøya - Isle.of.Canna    |  -1.0057 |     0.497  |   -2.022 |    1.000 |      No     |
  | Bjørnøya - Jan.Mayen        |  -0.0582 |     0.270  |   -0.216 |    1.000 |      No     |
  | Bjørnøya - Jarsteinen       |  -1.7913 |     0.390  |   -4.594 |   0.0003 |     Yes     |
  | Bjørnøya - Little.Saltee    |  -3.2406 |     0.307  |  -10.562 |   <.0001 |     Yes     |
  | Combined - Eynhallow        |  -2.9356 |     0.195  |  -15.064 |   <.0001 |     Yes     |
  | Combined - Faroe.Islands    |  -1.8585 |     0.317  |   -5.869 |   <.0001 |     Yes     |
  | Combined - Inishkea         |  -3.1773 |     0.659  |   -4.819 |   0.0001 |     Yes     |
  | Combined - Isle.of.Canna    |  -2.5097 |     0.480  |   -5.229 |   <.0001 |     Yes     |
  | Combined - Jan.Mayen        |  -1.5622 |     0.236  |   -6.620 |   <.0001 |     Yes     |
  | Combined - Jarsteinen       |  -3.2953 |     0.367  |   -8.970 |   <.0001 |     Yes     |
  | Combined - Little.Saltee    |  -4.7446 |     0.278  |  -17.092 |   <.0001 |     Yes     |
  | Eynhallow - Faroe.Islands   |   1.0771 |     0.294  |    3.665 |   0.0127 |     Yes     |
  | Eynhallow - Inishkea        |  -0.2417 |     0.649  |   -0.373 |    1.000 |      No     |
  | Eynhallow - Isle.of.Canna   |   0.4259 |     0.465  |    0.915 |    1.000 |      No     |
  | Eynhallow - Jan.Mayen       |   1.3734 |     0.204  |    6.721 |   <.0001 |     Yes     |
  | Eynhallow - Jarsteinen      |  -0.3597 |     0.348  |   -1.034 |    1.000 |      No     |
  | Eynhallow - Little.Saltee   |  -1.8090 |     0.251  |   -7.200 |   <.0001 |     Yes     |
  | Faroe.Islands - Inishkea    |  -1.3188 |     0.695  |   -1.897 |    1.000 |      No     |
  | Faroe.Islands - Isle.of.Canna | -0.6512 |    0.528  |   -1.233 |    1.000 |      No     |
  | Faroe.Islands - Jan.Mayen   |   0.2963 |     0.323  |    0.919 |    1.000 |      No     |
  | Faroe.Islands - Jarsteinen  |  -1.4367 |     0.428  |   -3.355 |   0.0393 |     Yes     |
  | Faroe.Islands - Little.Saltee | -2.8861 |    0.354  |   -8.149 |   <.0001 |     Yes     |
  | Inishkea - Isle.of.Canna    |   0.6676 |     0.783  |    0.853 |    1.000 |      No     |
  | Inishkea - Jan.Mayen        |   1.6151 |     0.662  |    2.439 |   0.6835 |      No     |
  | Inishkea - Jarsteinen       |  -0.1180 |     0.720  |   -0.164 |    1.000 |      No     |
  | Inishkea - Little.Saltee    |  -1.5673 |     0.678  |   -2.311 |   0.9613 |      No     |
  | Isle.of.Canna - Jan.Mayen   |   0.9475 |     0.484  |    1.958 |    1.000 |      No     |
  | Isle.of.Canna - Jarsteinen  |  -0.7855 |     0.560  |   -1.403 |    1.000 |      No     |
  | Isle.of.Canna - Little.Saltee | -2.2349 |    0.505  |   -4.421 |   0.0006 |     Yes     |
  | Jan.Mayen - Jarsteinen      |  -1.7331 |     0.372  |   -4.653 |   0.0002 |     Yes     |
  | Jan.Mayen - Little.Saltee   |  -3.1824 |     0.284  |  -11.194 |   <.0001 |     Yes     |
  | Jarsteinen - Little.Saltee  |  -1.4494 |     0.400  |   -3.622 |   0.0149 |     Yes     |
  
# Do a conditional boxplot for different months for each colony as in the months 10,11,12,1,2,3 on
# the X-axis with pers on the Y-axis



