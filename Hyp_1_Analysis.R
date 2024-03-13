# Analysis script for hypothesis 1: trying out a generalized mixed effects model (GLMM) using the lme4 package with a gamma / quasi-poisson/ gamma 
# distribution of course after checking if the requirements of each and whether the pers data fits that criteria. After selecting the distribution family
# also, choose the link function carefully. 
# I have decided to use the following syntax: glmer(pers ~ Colony + (1 | individ_id))

# One of the assumptions of the models is homogenity of variance. This has to come prior to model-assumption because if the 
# data is non-normal and doesn't have homogenity of variance, only then glmm over glm is decided. 
# Then check if the model residuals are normally distrbuted or not. If not, then Kruskal-Wallis test (if yes then anova), followed by appropraite post-hoc
# tests like Mann-Whitney u-test for parametric or dunn test for non-parametric. 

# Loading relevant packages----
library(dplyr)
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(lme4) # for glmm

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
mylocs <- readRDS(paste0(datadir, "/SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"/summaryTable.rds"))

# Merging and curating for NBS (Non-Breeding Season)----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"

# Creating a new combined dataset which merges tracks from Skalfandi and Langanes (distance between colonies: 130.429 km)
colonies_to_combine <- c("Skjalfandi", "Langanes")
nbs_mylocs$colony[nbs_mylocs$colony %in% colonies_to_combine] <- "Combined"

unique_colony_list <- sort(unique(nbs_mylocs$colony), decreasing = F)
# unique_colony_list_2 <- gsub(" ", ".", unique_colony_list)

# Preparing dataframe----

n_tracks <- nbs_mylocs %>% group_by(colony, individ_id) %>% summarise(Number_of_tracks = n())
colnames(n_tracks) <- c("Colony", "Individ_id", "Number_of_tracks")
n_tracks$Colony <- gsub(" ",".",n_tracks$Colony)

# sample_size <- nbs_mylocs %>% group_by(colony) %>% summzarise(Sample_size = length(unique(individ_id)))
# colnames(sample_size) <- c("Colony", "Sample_size")
# sample_size$Colony <- gsub(" ",".",sample_size$Colony)

# df <- merge(mean_n_tracks, sample_size, by = "Colony")

pers_df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Ind/outputs/csv/exposure_scores_by_individual.csv")
pers_df$population[pers_df$population %in% colonies_to_combine] <- "Combined"
pers_df

# sorted_pers_df <- pers_df[order(pers_df$population),]
# sorted_pers_df_2 <- sorted_pers_df %>% arrange(population, individual)
# sorted_pers_df_2
# unique_colony_list_2 == unique(sorted_pers_df_2$population)

individ_id_df <- as.data.frame(matrix(ncol = 3))
colnames(individ_id_df) <- c("population", "Individ_id", "individual")

for(i in unique_colony_list){ # First for loop starts
  sub <- nbs_mylocs[nbs_mylocs$colony == i,]
  individ_id_vector <- unique(sub$individ_id)
  sub_df <- data.frame(population = i, Individ_id = individ_id_vector, individual = 1:length(unique(sub$individ_id)))
  individ_id_df <- rbind(individ_id_df, sub_df)
} # First for loop ends

individ_id_df <- individ_id_df[-1,]
individ_id_df$population <- gsub(" ", ".", individ_id_df$population)
individ_id_df

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/csv/")
hyp_1_df <- merge(pers_df, individ_id_df, by = c("population", "individual"))
hyp_1_df <- hyp_1_df[,-c(2:4)]
colnames(hyp_1_df) <- c("Colony", "pers", "Individ_id")

# Now just have to add no of tracks dataframe to this
final_hyp_1_df <- merge(hyp_1_df, n_tracks, by = c("Colony", "Individ_id"))
df <- final_hyp_1_df

# Saving new input dataframes----
write.csv(final_hyp_1_df, "ntracks_glmer_df.csv")
write.csv(hyp_1_df, "glmer_df.csv")

# # Final df----
# df <- merge(nbs_mylocs, df_transient, by = c("colony", "individ_id"))
# df

# colnames(pers_df) <- c("Colony", "pers")


# df <- merge(df, pers_df, by = "Colony")

# GLMM----

# For checking if Poisson's regression should be used, check pers mean and variance. If they aren't same, poisson's won't be a good fit.
# Most ecological data is overdispersed, i.e. variance is greater than mean. 

# mean_pers <- df %>% group_by(population) %>% summarise(mean_exposure_score = mean(exposure_score))
# variance_pers <- df %>% group_by(population) %>% summarise(variance_exposure_score = var(exposure_score))
# poisson_testing <- merge(mean_pers, variance_pers, by = "population")
# poisson_testing

# Homoskedasticity assumption testing----
skedasticity_test_model <- lm(pers ~ Colony, data = df)
residuals <- resid(skedasticity_test_model)
fitted_values <- fitted(skedasticity_test_model)
plot(fitted_values, residuals, main = "Residuals vs Fitted",
     xlab = "Fitted values", ylab = "Residuals")

# Function to perform Breusch-Pagan test for heteroscedasticity

breusch_pagan_test <- function(model) {
  library(lmtest)
  bptest(model)
}

# Function to perform White test for heteroscedasticity
white_test <- function(model) {
  library(lmtest)
  bptest(model, ~ fitted(model) + I(fitted(model)^2))
}

# Perform Breusch-Pagan test
breusch_pagan_result <- breusch_pagan_test(skedasticity_test_model)
print(breusch_pagan_result)

# Perform White test
white_test_result <- white_test(skedasticity_test_model)
print(white_test_result)

# Baseline model----
baseline_model <- lmer(pers ~ 1 + (1 | Colony), data = df, ) # No predictors in the baseline model, this is an intercept-only model
# meaning we're assuming only the intercept to change and the mean to remain the same across individuals; random effects are to be stated 
# on the left side of the vertical bar; by this line we say that each individ_id's pers differs from the overall mean of all individual pers;
# the repeated measurements are within the individ_id variable (cluster variable)
print(baseline_model)
summary(baseline_model) # gives coefficients, standard errors, significance tests, and goodness-of-fit measures.
# plot(simulateResiduals(baseline_model)) # used to check goodness of fit
plot(residuals(baseline_model))
drop1(baseline_model, test = "Chisq") # use for other models with predictors

hist(df$pers, breaks = 30)
# Fit models with different distribution families----

# model_gaussian <- glmer(pers ~ Colony + (1 | Colony / Individ_id), data = df, family = gaussian(link = "identity"))
View(df)
model_gamma <- glmer(pers ~ 1 + (1 | Colony / Individ_id), data = df, family = Gamma(link = "identity"))
model_gamma_2 <- glmer(pers ~ Colony + (1 | Individ_id / Number_of_tracks) + (1 | Colony / Individ_id), data = df, family = Gamma(link = "identity"))

plot(model_gamma)

model_inverse_gaussian <- glmer(pers ~ 1 + (1 | Individ_id / Number_of_tracks) + (1 | Colony / Individ_id), data = df, family = inverse.gaussian(link = "identity"))
plot(model_inverse_gaussian)

# Compute AIC for each model
AIC_values <- c(AIC(model_gamma), AIC(model_inverse_gaussian))
BIC_values <- c(BIC(model_gamma), BIC(model_inverse_gaussian))


# Compare AIC and BIC values
models <- c("Gamma", "Inverse Gaussian")
comparison <- data.frame(models, AIC_values, BIC_values)
print(comparison)

print(model_inverse_gaussian)

# Use confint(model_name) for calculating confidence intervals based on likelihood-ratio based tests (LRT) for a model; methods include method = "Wald"
# but this method doesn't have enough power compard to lrt but it's faster to compute. Bootstrap intervals: method = "boot", nsim = 1000 (up to you)

# Need to change things around here since sample size isn't a necessary covariate. 
# Do you expect the slope to change with the #tracks too? From what I 
# know, as #tracks increases, home range volume decreases so pers should decrease in a way. So, slope
# remains the same I think, I mean it's so white and black but for now, the slope can be assumed to be same
# across diff #tracks. 




# Residual analysis
residuals <- resid(model_inverse_gaussian)
fitted_values <- fitted(model_inverse_gaussian)
residual_plot <- ggplot(data.frame(residuals, fitted_values), aes(x = fitted_values, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Fitted Values") +
  ylab("Residuals") +
  ggtitle("Residuals vs Fitted Values")

print(residual_plot)

# Calculate deviance
deviance <- deviance(model_inverse_gaussian)

# Print AIC, BIC, and deviance
cat("AIC:", AIC(model_inverse_gaussian), "\n")
cat("BIC:", BIC(model_inverse_gaussian), "\n")
cat("Deviance:", deviance, "\n")

# Check the random effects
ranef(model_inverse_gaussian)

# Check fixed effects
fixef(model_inverse_gaussian)

# Extracting residuals from the model
residuals <- resid(model_inverse_gaussian)

# Creating a residual plot
plot(fitted(model_inverse_gaussian), residuals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residual Plot")

# Adding a horizontal line at y = 0 for reference
abline(h = 0, col = "red")

# # Adding a lowess smooth line to identify any patterns
# lines(lowess(fitted(model_inverse_gaussian), residuals), col = "blue")


trial_model <- glmer(pers ~ Colony + (1 | Colony / Individ_id), family = inverse.gaussian(link = "identity"), data = df)
ggplot(df, aes(x = Colony, y = pers)) +
  geom_boxplot() +
  labs(x = "Colony", y = "pers") +
  ggtitle("Boxplot of pers by Colony")

ggplot(df, aes(x = Colony, y = pers)) +
  geom_point() +
  labs(x = "Colony", y = "pers") +
  ggtitle("Scatterplot of pers by Colony")

summary_df <- df %>%
  group_by(Colony) %>%
  summarise(mean_pers = mean(pers), 
            se_pers = sd(pers) / sqrt(n()))

# Create scatterplot with standard errors
ggplot() +
  geom_point(data = df, aes(x = Colony, y = pers)) +
  geom_errorbar(data = summary_df, aes(x = Colony, ymin = mean_pers - se_pers, ymax = mean_pers + se_pers), 
                width = 0.2, color = "blue") +
  labs(x = "Colony", y = "pers") +
  ggtitle("Scatterplot of pers by Colony with Standard Errors of Mean")

shapiro.test(residuals)
kruskal.test(df$pers ~ df$Colony)

library(dunn.test)

# Perform post hoc pairwise comparisons using Dunn test with Bonferroni correction
posthoc_dunn <- dunn.test(df$pers, df$Colony, method = "bonferroni")

# Print the post hoc test results
print(posthoc_dunn)


# 






# Starting analysis of hypothesis 1 anew----

library(lme4)

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/csv/")
df_1 <- read.csv("ntracks_glmer_df.csv")
df_1 <- df_1[,-1]
colonies_to_combine <- c("Skjalfandi", "Langanes")
df_1$Colony[df_1$Colony %in% colonies_to_combine] <- "Combined"
df_2 <- read.csv("glmer_df.csv")
df_2 <- df_2[,-c(1,4,5)]
colnames(df_2) <- c("Colony", "Individual", "pers", "Individ_id")
df_3 <- merge(df_1, df_2, by = c("Colony", "Individ_id", "pers"))

df_4 <- read.csv("nas_ind.csv")
df_4[c("Colony", "Individual")] <- do.call(rbind, strsplit(as.character(df_4$name), "_", fixed = T))
df_4 <- df_4[,-1]
df_4$Colony[df_4$Colony %in% colonies_to_combine] <- "Combined"

df_5 <- merge(df_3, df_4, by = c("Colony", "Individual"))
df <- df_5
write.csv(df_5, "Final_hypothesis_1_analysis_df.csv")

# New glmm model

full_model_gamma <- glmer(pers ~ Colony + (1 | Individ_id), data = df, family = Gamma(link = "identity"))
full_model_inverse_gaussian <- glmer(pers ~ Colony + (1 | Individ_id), data = df, family = inverse.gaussian(link = "identity"))

print(full_model_inverse_gaussian)


