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

# Exploring data----

# Scatterplot of colony vs. pers ----

ggplot(analysis_df, aes(x = pers, y = Colony, col = Colony)) +
  geom_jitter(width = 0, height = 0.2) +  # Add jitter for better visualization
  scale_color_viridis_d() +  # Use the viridis color palette
  labs(x = "Plastic Exposure Risk Score", y = "Colony") +
  theme_minimal()

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
glmm_model_3 <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df_2, control = control_params)
AIC(glmm_model_3) # 2485.363; lower than glmm_model_1 
summary(glmm_model_3)
qqnorm(resid(glmm_model_3))
qqline(resid(glmm_model_3)) # non-normality remains
hist(resid(glmm_model_3), breaks=90)
analysis_df_2$resid = resid(glmm_model_3)
boxplot(analysis_df_2$resid~analysis_df_2$Colony) # More or less homoskedastic
ggplot(analysis_df_2, aes(x = Colony, y = resid)) + geom_boxplot() 

# Will skew interpretability but trying this out with diff link functions
inverse_glmm_model_3 <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "inverse"), data = analysis_df_2, control = control_params) # didn't converge
log_glmm_model_3 <- glmer(pers ~ Colony + (1 | Individ_id), family = Gamma(link = "log"), data = analysis_df_2, control = control_params) # didn't converge

# Square-root trasnforming pers
glmm_model_4 <- glmer(sqrt(pers) ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df, control = control_params)
AIC(glmm_model_4) # -979.8162

# Incorporating both, log-transforming and removing Alkefjellet
glmm_model_5 <- glmer(log(pers) ~ Colony + (1 | Individ_id), family = Gamma(link = "identity"), data = analysis_df_2, control = control_params)
AIC(glmm_model_5) # -3272.464; the most favourable model so far

# Hence, if I want a log-link function, log_glmm is my safest bet since it doesn't remove Alkefjellet and models data correctly


