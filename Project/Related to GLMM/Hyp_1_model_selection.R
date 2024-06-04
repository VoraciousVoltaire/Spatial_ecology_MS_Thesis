# From the basics---- 

# Installing relevant packages ----
library(tidyverse)
library(gamm4)
library(mgcv)
library(lme4)
library(gratia)
library(plotly)
library(ggeffects)
library(viridis)
library(viridisLite)
library(car) # for ANOVA
library(multcomp) # for glht
library(lsmeans) # for lsm
library(rcompanion) # for the function nagelkerke

# The code for gam is derived from Wood's book on GAMs and from the following website: https://m-clark.github.io/generalized-additive-models/application.html

# Loading data ----

# Loading correct data----

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
analysis_df_2 <- analysis_df[!(analysis_df$Colony == "Alkefjellet"),] # removing Alkefjellet
analysis_df_3 <- analysis_df_2 %>% filter(!(Colony == "Inishkea" | Colony == "Isle of Canna"),) # removing Alkefjellet, Inishkea and Isle of Canna
analysis_df$Colony <- as.factor(analysis_df$Colony)
analysis_df$Individ_id <- as.factor(analysis_df$Individ_id)
analysis_df_2$Colony <- as.factor(analysis_df_2$Colony)
analysis_df_2$Individ_id <- as.factor(analysis_df_2$Individ_id)
analysis_df_3$Colony <- as.factor(analysis_df_3$Colony)
analysis_df_3$Individ_id <- as.factor(analysis_df_3$Individ_id)

# Exploring data----

# Scatterplot of colony vs. pers ----


jitter_plot <- ggplot(analysis_df, aes(x = pers, y = Colony, col = Colony)) +
  geom_jitter(width = 0, height = 0.2) +  # Add jitter for better visualization
  scale_color_viridis_d() +  # Use the viridis color palette
  labs(x = "Plastic Exposure Risk Score (PERS)", y = "Colony") +
  theme_minimal()

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Images/")
ggsave()

linear_model <- lm(pers ~ Colony, data = analysis_df)
summary(linear_model)
AIC(linear_model) # 3500.205

library(mgcv)
gam_model <- gam(pers ~ Colony, data = analysis_df)
summary(gam_model)
AIC(gam_model) # 3500.205
summary(gam_model)$sp.criterion # GCV.Cp 
# 2.419914 
summary(gam_model)$r.sq # 0.4688855

analysis_df$Colony <- as.numeric(as.factor(analysis_df$Colony)) # this step makes this analysis seem wrong since the gap between orders is not consistent
gam_model_2 <- gam(pers ~ s(Colony, bs = "cr"), data = analysis_df)
summary(gam_model_2)
AIC(gam_model_2) # 3499.806
plot(gam_model_2, col = "blue", shade = TRUE, shade.col = "#CCCCCC", seWithMean = T)
summary(gam_model_2)$sp.criterion #  GCV.Cp 
# 2.418866 
summary(gam_model_2)$r.sq # 0.4689038
anova(gam_model, gam_model_2, test = "Chisq")
coef(gam_model_2)

# Visualisation ----
plot(ggeffects::ggpredict(gam_model_2), facets = TRUE)
gratia::draw(gam_model_2)

# Fit a GAM with dummy variables for Colony
analysis_df$Colony <- factor(analysis_df$Colony)
gam_model_3 <- gam(pers ~ s(Colony, bs = "re"), data = analysis_df)
AIC(gam_model_3)  # 3499.887
summary(gam_model_3)$sp.criterion # GCV.Cp 
# 2.419079 
summary(gam_model_3)$r.sq # 0.4688905
anova(gam_model, gam_model_3, test = "Chisq")
plot(ggeffects::ggpredict(gam_model_3), facets = TRUE)
gratia::draw(gam_model_3)

# GAMM finally----
analysis_df$Colony <- as.numeric(as.factor(analysis_df$Colony)) 
gamm_model <- gamm4(pers ~ s(Colony, bs = "cr"), family = Gamma(link = "identity"), data = analysis_df, random = ~ (1 | Individ_id))
summary(gamm_model)
plot(gamm_model$mer)
plot(gamm_model$gam)
AIC_mer <- AIC(gamm_model$mer)
gam_model_gam_only <- gam(pers ~ s(Colony, bs = "cr"), family = Gamma(link = "identity"), data = analysis_df)
AIC_gam <- AIC(gam_model_gam_only)
AIC_total <- AIC_mer + AIC_gam
AIC_total # don't know if this correct but even if it is, then doesn't matter since I'm choosing glmm over gamm

# Fitting a glmm finally----

control_params <- glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 100000))  # Maximum number of function evaluations

glmm_model <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df, control = control_params)
AIC(glmm_model) # 3453.074
summary(glmm_model)
qqnorm(resid(glmm_model))
qqline(resid(glmm_model))
plot(glmm_model)
hist(resid(glmm_model), breaks=90)
analysis_df$resid = resid(glmm_model)
boxplot(analysis_df$resid~analysis_df$Colony)
ggplot(analysis_df, aes(x = Colony, y = resid)) + geom_boxplot() 

# Trying out different link functions----
inverse_glmm <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "inverse"), data = analysis_df, control = control_params) # didn't converge
log_glmm <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "log"), data = analysis_df, control = control_params)
AIC(log_glmm) # 3383.849
plot(log_glmm)
hist(resid(log_glmm), breaks = 90)
qqnorm(resid(log_glmm))
qqline(resid(log_glmm))
analysis_df$resid = resid(log_glmm)
boxplot(analysis_df$resid~analysis_df$Colony)
ggplot(analysis_df, aes(x = Colony, y = resid)) + geom_boxplot() 

# Comparing it with other models----

# Choosing random effect
trial_mod_1 = gls(pers ~ 1 + Colony, data = analysis_df, method = "REML")
AIC(trial_mod_1) # 3974.622
trial_mod_2 = lmer(pers ~ Colony + (1 | Individ_id), data = analysis_df,REML = T)
AIC(trial_mod_2) # 3607.024; this tells us that (1 | Individ_id) should be added as a random effect


# Choosing fixed effect
trial_mod_3 = lmer(pers ~ 1 + (1 | Individ_id), data = analysis_df,REML = F)
AIC(trial_mod_3) # 3673.142
trial_mod_4 = lmer(pers ~ Colony + (1 | Individ_id), data = analysis_df,REML = F)
AIC(trial_mod_4) # 3599.283; thus Colony is a suitable fixed effect

# Reaffirming the significance of including colony as a fixed effect----
full_model <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df, control = control_params)
null_model <- glmer(pers ~ 1 + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df, control = control_params)
AIC(null_model)
summary(null_model) # Model 1.0 in results
lrt_result <- anova(full_model, null_model) # Likelihood Ratio Test
print(lrt_result)

# trying out variations of glmm_model to improve fit----

# Log-transforming pers
glmm_model_2 <- glmer(log(pers) ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df, control = control_params)
AIC(glmm_model_2) # -2597.023
qqnorm(resid(glmm_model_2))
qqline(resid(glmm_model_2)) # didn't improve normality to any significant extent
plot(glmm_model_2)
hist(resid(glmm_model_2), breaks=90)
analysis_df$resid = resid(glmm_model_2)
boxplot(analysis_df$resid~analysis_df$Colony)
ggplot(analysis_df, aes(x = Colony, y = resid)) + geom_boxplot() 

# Will skew interpretability but trying this out with diff link functions
inverse_glmm_model_2 <- glmer(log(pers) ~ Colony + (1 | Individ_id), family = Gamma(link = "inverse"), data = analysis_df, control = control_params) # didn't ccnverge
log_glmm_model_2 <- glmer(log(pers) ~ Colony + (1 | Individ_id), family = Gamma(link = "log"), data = analysis_df, control = control_params) # didn't converge

# Removing Alkefjellet from the dataset
View(analysis_df_2)
glmm_model_3 <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df_2, control = control_params)
AIC(glmm_model_3) # 2485; lower than glmm_model_1 
summary(glmm_model_3)
shapiro.test(resid(glmm_model_3))
qqnorm(resid(glmm_model_3))
qqline(resid(glmm_model_3)) # non-normality remains
hist(resid(glmm_model_3), breaks=90)
analysis_df_2$resid = resid(glmm_model_3)
boxplot(analysis_df_2$resid~analysis_df_2$Colony) # More or less homoskedastic
ggplot(analysis_df_2, aes(x = Colony, y = resid)) + geom_boxplot() 

glm <- glm(pers ~ Colony, family = Gamma(link = "identity"), data = analysis_df_2)
AIC(glm)
summary(glm)

# Will skew interpretability but trying this out with diff link functions
inverse_glmm_model_3 <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "inverse"), data = analysis_df_2, control = control_params) # didn't converge
log_glmm_model_3 <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "log"), data = analysis_df_2, control = control_params) # didn't converge

# Square-root trasnforming pers
glmm_model_4 <- glmer(sqrt(pers) ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df, control = control_params)
AIC(glmm_model_4) # -979.8162

# Incorporating both, log-transforming and removing Alkefjellet
glmm_model_5 <- glmer(log(pers) ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df_2, control = control_params)
AIC(glmm_model_5) # -3272.464; the most favourable model so far

Anova(log_glmm, type = 3)
plot(analysis_df$pers)
residuals_log_glmm <- residuals(log_glmm)

# Plot residuals against fitted values
plot(fitted(log_glmm), residuals_log_glmm, 
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")

# Add a horizontal line at y = 0 for reference
abline(h = 0, col = "red", lty = 2)

# Running a GEE
# Remove prefixes from Individ_id
analysis_df$Individ_id <- as.numeric(gsub("[^0-9]", "", analysis_df$Individ_id))
# Fit GEE model after converting Individ_id to numeric
gee_model <- geeglm(pers ~ Colony, 
                    id = Individ_id,
                    family = Gamma(link = "log"),
                    corstr = "exchangeable",
                    data = analysis_df)
summary(gee_model)
plot(gee_model, which = 1)  # Residuals vs. Fitted for homogeneity of variance
hist(resid(gee_model), breaks = 90)
qqnorm(resid(gee_model))
qqline(resid(gee_model))
analysis_df$resid = resid(gee_model)
boxplot(analysis_df$resid~analysis_df$Colony) 

# Final model----
# After careful consideration, I have decided to drop Alkefjellet from the dataset to maintain homogenity of variance across
# all predictor variables: glmm_model_3

final_model <- glmm_model_3

null_model <- glmer(pers ~ 1 + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df_2, control = control_params)
nagelkerke_results <- nagelkerke(glmm_model_3, null = null_model)
print(nagelkerke_results) # variance explained by fixed effects

# Parametric bootstrap to test whether the random effect has a significant effect
nBoot <- 1000
lrStat <- rep(NA, nBoot)

ft.null <- lm(pers~Colony, data = analysis_df) # null model
ft.alt <- lme(pers~Colony, random = ~1|Individ_id, data = analysis_df) # alternate model
ft.alt <- try(lme(pers ~ Colony, random = ~1 | Individ_id, data = analysis_df, method = "ML"))
if (inherits(ft.alt, "try-error")) {
  # Handle convergence issues
  warning("Convergence issue in fitting the alternate model with default options. Trying different optimization algorithm.")
  
  # Try fitting the alternate model with a different optimization algorithm
  ft.alt <- try(lme(pers ~ Colony, random = ~1 | Individ_id, data = analysis_df, method = "REML"))
  
  if (inherits(ft.alt, "try-error")) {
    # Handle convergence issues
    warning("Convergence issue in fitting the alternate model with different optimization algorithm. Results may not be reliable.")
  }
}
lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) # observed test stat

for(iBoot in 1:nBoot){
  analysis_df$pers <- unlist(simulate(ft.null)) # resampled data
  bNull <- lm(pers~Colony, data = analysis_df) # null model
  bAlt <- lme(pers~Colony, random = ~1|Individ_id, data = analysis_df) # alternate model
  lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull) # resampled test stat
}

mean(lrStat > lrObs) # p-value for test of the random effect
hist(lrStat, col = "blue")
abline(v = lrObs, col = "red", lwd = 3, lty = 3)


# Post-hoc test----

summary(final_model)
library(emmeans)
posthoc <- emmeans(final_model, pairwise ~ Colony, median = T)
summary(posthoc)
Tukey_HSD_results <- summary(pairs(posthoc))
Tukey_HSD_results

# Filter pairs with non-significant differences (p > 0.05)
non_sig_pairs <- Tukey_HSD_results[Tukey_HSD_results$p.value > 0.05, c("contrast", "p.value")]
nrow(Tukey_HSD_results)
# Print the table
non_sig_pairs
nrow(non_sig_pairs)


# Conf ints for hyp 1 and 2 PERS

conf_df <- read.csv("correct_exposure_scores_by_individual_according_to_tracking_yr.csv")
# Calculate median and 95% confidence interval for each colony
colony_medians <- tapply(conf_df$exposure_score, conf_df$population, function(x) {
  med <- median(x)
  ci <- quantile(x, probs = c(0.025, 0.975))
  return(c(median = med, lower = ci[1], upper = ci[2]))
})

# Convert the results to a dataframe
colony_medians_df <- do.call(rbind, colony_medians)

# Print the results
write.csv(colony_medians_df,"colony_PERS_confints.csv")

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/")
conf_df_PAME <- read.csv("correct_exposure_scores_by_month.csv")

# Calculate median and 95% confidence interval for each colony
PAME_medians <- tapply(conf_df_PAME$exposure_score, conf_df_PAME$population, function(x) {
  med <- median(x)
  ci <- quantile(x, probs = c(0.025, 0.975))
  return(c(median = med, lower = ci[1], upper = ci[2]))
})

# Convert the results to a dataframe
PAME_medians_df <- do.call(rbind, PAME_medians)
View(PAME_medians_df)

# Print the results
write.csv(PAME_medians_df,"PAME_PERS_confints.csv")

Anova(final_model, type = "2")






