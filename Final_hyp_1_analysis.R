# New Analaysis for hypothesis 1: this reminds me, run the months script again with hr[>=95]----

# Prerequisites----
knitr::opts_chunk$set(echo = TRUE,
                      results = 'show',
                      message = FALSE,
                      warning = FALSE ) # this sets settings for all the script

# Loading essential packages----
library(nlme) # for gls function 
library(lme4) # for glmer function 
library(ggplot2)
library(dplyr)
library(pander)
library(car)
library(lmerTest)
library(MuMIn)
install.packages("emmeans")
library(emmeans) # to generate specific contrasts and tests for various comparisons over time and/or between groups

# Some exploratory data analysis first----

# df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/csv/Final_hypothesis_1_analysis_df.csv")
# df # has colony name, individual number, individual ID, pers, number of tracks per individual to 
# calculate its KDE, vals, nas and percent_na. 

# Loading correct data----

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/latest_right_attempt/outputs/csv")
df <- read.csv("correct_exposure_scores_by_individual_according_to_tracking_yr.csv")

nas_df <- read.csv("correct_nas_ind.csv")
View(nas_df)

View(nas_df[order((nas_df$percent_na), decreasing = T),])

corrected_df <- df[df$exposure_score > 10,]

# Eynhallow_outlier <- df_mod[df_mod$colony == "Eynhallow" & df_mod$individ_id == "GBT-FH89021" & df_mod$tracking_year == 2016,]
# plot(Eynhallow_outlier$lon, Eynhallow_outlier$lat)
# plot(land, add = T)
# plot(plastics)

# Creating a sample size dataframe first and removing colony_tracking year with less than 50 relocations
sample_size_df <- df_mod %>% group_by(colony, individ_id, tracking_year) %>% summarise(n = n())
sample_size_df[sample_size_df$n < 50,]
# Jan Mayen NOS-4181597          2018    46; so remove this and Eynhallow GBT-FH89021

analysis_df <- corrected_df[!(corrected_df$population == "Jan.Mayen" & 
                                corrected_df$individual == "NOS.4181597" & 
                                corrected_df$tracking_year == 2018), ]
View(analysis_df[analysis_df$population == "Jan.Mayen", ])

analysis_df <- analysis_df[,-c(4,5)]
colnames(analysis_df) <- c("Colony", "Individ_id", "Tracking_year", "pers")

# 1. Checking for outliers----
par(mfrow = c(1,1))
dotchart(df$pers) # a left tail, a few outliers
df[df$pers < 15,]

# 2. Checking for homogenity of variance----
ggplot(df, aes(x=Colony, y = pers)) + geom_boxplot() # + facet_grid(cols=vars(Individ_id))
# more or less sequal spread, some have higher disparity than others

# 3. Checking for normality----
ggplot(df, aes(x=pers)) + geom_histogram(binwidth = 0.2) # doesn't look like a bell shaped curve, has a left tail
shapiro.test(df$pers)
# 4. GLMM----

# Selecting random effect----
modGLS = gls(pers ~ 1 + Colony, data = df, method = "REML")
AIC(modGLS) # 1312.738

modA = lmer(pers ~ Colony + (1 | Individ_id), data = df,REML = T)
AIC(modA) # 1313.102

# this tells us that probably (1 | Individ_id) isn't that significant after all

# Selecting fixed effects----
mod0 = lmer(pers ~ 1 + (1 | Individ_id), data = df,REML = F)
AIC(mod0) # 1613.719

mod1 = lmer(pers ~ Colony + (1 | Individ_id), data = df,REML = F)
AIC(mod1) # 1304.388

# Colony as a fixed effect decreases the model AIC 

# To obtain paramters of the best model----

mod1 = lmer(pers ~ Colony + (1 | Individ_id), data = df,REML = T)
AIC(mod1) #1313.102
summary(mod1)

confint(mod1, oldNames=FALSE)

plot(mod1)
qqnorm(resid(mod1))
qqline(resid(mod1)) # not great but ok 


hist(resid(mod1),breaks=20)


df$resid = resid(mod1)

boxplot(df$resid~df$Colony)
ggplot(df, aes(x=Colony, y = resid)) + geom_boxplot() # + facet_grid(cols = vars(species))

vif(mod1) # Error in vif.merMod(mod1) : model contains fewer than 2 terms

r.squaredGLMM(mod1)

# Now trying this same thing with glmm----

ig_mod_0 <- glmer(pers ~ 1 + (1 | Individ_id), data = df, family = inverse.gaussian(link = "identity"))
AIC(ig_mod_0) # 749.1618

ig_mod_1 <- glmer(pers ~ Colony + (1 | Individ_id), data = df, family = inverse.gaussian(link = "identity"))
# Warning message:
  # In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
                 # Model failed to converge with max|grad| = 0.0268077 (tol = 0.002, component 1)

# Increasing number of itertaions

# Specify control parameters
control_params <- glmerControl(optimizer = "bobyqa",  # Choose optimizer (optional)
                               optCtrl = list(maxfun = 100000))  # Maximum number of function evaluations

# Fit the GLMM with increased max iterations
ig_mod_1 <- glmer(pers ~ Colony + (1 | Individ_id), 
                  data = df, 
                  family = inverse.gaussian(link = "identity"),
                  control = control_params)  # Specify control parameters
AIC(ig_mod_1) # 600.5235; best model 
BIC(ig_mod_1) # 648.0555

ig_mod_2 <- glmer(pers ~ Colony + (1 | Individ_id) + (1 | percent_na), data = df, family = inverse.gaussian(link = "identity"), control = control_params)
AIC(ig_mod_2) # 602.5235

ig_mod_3 <- glmer(pers ~ Colony + (1 | Individ_id) + (1 | Number_of_tracks), data = df, family = inverse.gaussian(link = "identity"), control = control_params)
AIC(ig_mod_3) # 602.5235

ig_mod_4 <- glmer(pers ~ Colony + (1 | Individ_id) + (1 | Number_of_tracks) + (1 | percent_na), data = df, family = inverse.gaussian(link = "identity"), control = control_params)
AIC(ig_mod_4) # 604.5235

ig_mod_5 <- glmer(pers ~ Colony + (1 | Individ_id) + (1 | vals), data = df, family = inverse.gaussian(link = "identity"), control = control_params)
AIC(ig_mod_5) # doesn't make sense

ig_mod_6 <- glmer(pers ~ Colony + (1 | percent_na), data = df, family = inverse.gaussian(link = "identity"), control = control_params)
AIC(ig_mod_6) # 756.2511

ig_mod_7 <- glmer(pers ~ Colony + (1 | Number_of_tracks), data = df, family = inverse.gaussian(link = "identity"), control = control_params)
AIC(ig_mod_7) # 1402.586

# To obtain paramters of the best model----

hist(df$pers, breaks = 30)
names(df)[1] <- "Colony"
names(df)[2] <- "Individ_id"
names(df)[6] <- "pers"
View(df)
mod1 = glmer(pers ~ Colony + (1 | Individ_id), data = df, family = Gamma(link = "identity"), control = control_params)
AIC(mod1) # 599.9103
BIC(mod1) # 647.4424
summary(mod1)


# confint(mod1, oldNames=FALSE)

plot(mod1)
qqnorm(resid(mod1))
qqline(resid(mod1)) # not great but ok 


hist(resid(mod1),breaks=90)
shapiro.test(resid(mod1))

df$resid = resid(mod1)

boxplot(df$resid~df$Colony)
ggplot(df, aes(x=Colony, y = resid)) + geom_boxplot() # + facet_grid(cols = vars(species))

vif(mod1) # Error in vif.merMod(mod1) : model contains fewer than 2 terms

r.squaredGLMM(mod1)

# Organizing findings----

# So, analysis_df is the df to use now

# 1. Checking for outliers----

par(mfrow = c(1,1))
dotchart(analysis_df$pers)


# 2. Checking for homogenity of variance----

ggplot(analysis_df, aes(x= pers, y = Colony)) + geom_boxplot()

# 3. Checking for normality----

ggplot(analysis_df, aes(x=pers)) + geom_histogram(binwidth = 0.2) # has a left tail; probably due to Alkefjellet
analysis_df_2 <- analysis_df[!(analysis_df$Colony == "Alkefjellet"),]
ggplot(analysis_df_2, aes(x=pers)) + geom_histogram(binwidth = 0.2) 

# GLMM----

# 1. Selecting random effect

modGLS = gls(pers ~ 1 + Colony, data = analysis_df, method = "REML")
AIC(modGLS) # 3512.466

modA = lmer(pers ~ Colony + (1 | Individ_id), data = analysis_df, REML = T)
AIC(modA) # 3371.769, so relatively lower AIC than modGLS

# 2. Selecting fixed effect

mod0 = lmer(pers ~ 1 + (1 | Individ_id), data = analysis_df, REML = F)
AIC(mod0) # 3673.142

mod1 = lmer(pers ~ Colony + (1 | Individ_id), data = analysis_df, REML = F)
AIC(mod1) # 3363.664; lower AIC than mod0

# Trying out lmer

# To obtain parameters of the best model (analysis_df first)----

mod1 = lmer(pers ~ Colony + (1 | Individ_id), data = analysis_df, REML = T)
AIC(mod1) # 3371.769
BIC(mod1) # 3429.919
summary(mod1)

confint(mod1, oldNames=FALSE)

plot(mod1)
qqnorm(resid(mod1))
qqline(resid(mod1)) 

hist(resid(mod1),breaks=20)

analysis_df$resid = resid(mod1)

boxplot(analysis_df$resid~analysis_df$Colony)
ggplot(analysis_df, aes(x=Colony, y = resid)) + geom_boxplot() # + facet_grid(cols = vars(species))

vif(mod1) # Error in vif.merMod(mod1) : model contains fewer than 2 terms

r.squaredGLMM(mod1)
# R2m       R2c
# [1,] 0.4800856 0.6897765

# Now with analysis_df_2
analysis_df_2
mod2 = lmer(pers ~ Colony + (1 | Individ_id), data = analysis_df_2, REML = T)
AIC(mod2) # 2716.119
summary(mod2)

confint(mod2, oldNames=FALSE)

plot(mod2)
qqnorm(resid(mod2))
qqline(resid(mod2)) # not great but ok 

hist(resid(mod2),breaks=20)
shapiro.test(resid(mod2))

analysis_df_2$resid = resid(mod2)

boxplot(analysis_df_2$resid~analysis_df_2$Colony)
ggplot(analysis_df_2, aes(x=Colony, y = resid)) + geom_boxplot() # + facet_grid(cols = vars(species))

vif(mod2) # Error in vif.merMod(mod1) : model contains fewer than 2 terms

r.squaredGLMM(mod2)
# R2m       R2c
# [1,] 0.5230997 0.7766151

# Now trying analysis_df_1 and analysis_df_2 with glmer

g_mod_1 <- glmer(pers ~ Colony + (1 | Individ_id), 
                  data = analysis_df, 
                  family = Gamma(link = "reciprocal"),
                 control = control_params)  # Specify control parameters
AIC(g_mod_1) # 3347.626
BIC(g_mod_1) # 3405.776
plot(g_mod_1)
qqnorm(resid(g_mod_1))
qqline(resid(g_mod_1)) # not great but ok 

# Using the log link function
g_mod_1 <- glmer(pers ~ Colony + (1 | Individ_id), 
                 data = analysis_df, 
                 family = Gamma(link = "log"),
                 control = control_params)
AIC(g_mod_1) # 3383.849
BIC(g_mod_1) # 3441.999

# Using the inverse link function (approximation)
g_mod_1 <- glmer(pers ~ Colony + (1 | Individ_id), 
                 data = analysis_df, 
                 family = Gamma(link = "inverse"),
                 control = control_params)

# Using the square-root link function
g_mod_1 <- glmer(pers ~ Colony + (1 | Individ_id), 
                 data = analysis_df, 
                 family = Gamma(link = "sqrt"),
                 control = control_params)

AIC(g_mod_1) # 3364.964
BIC(g_mod_1) # 3423.114


hist(resid(g_mod_1),breaks=40)

analysis_df$resid = resid(g_mod_1)

boxplot(analysis_df$resid~analysis_df$Colony)
ggplot(analysis_df, aes(x=Colony, y = resid)) + geom_boxplot() # + facet_grid(cols = vars(species))


ig_mod_1 <- glmer(pers ~ Colony + (1 | Individ_id), 
                  data = analysis_df, 
                  family = inverse.gaussian(link = "identity"),
                  control = new_control_params)  # Specify control parameters

g_mod_2 <- glmer(pers ~ Colony + (1 | Individ_id), 
                  data = analysis_df_2, 
                  family = Gamma(link = "identity"),
                 control = control_params)  # Specify control parameters
AIC(g_mod_2) # 2485.363
BIC(g_mod_2) # 2538.104

ig_mod_2 <- glmer(pers ~ Colony + (1 | Individ_id), 
                 data = analysis_df_2, 
                 family = inverse.gaussian(link = "identity"),
                 control = control_params)  # Specify control parameters
AIC(ig_mod_2) # 2528.119
BIC(ig_mod_2) # 2580.859

# g_mod_2: lowest AIC and BIC

summary(g_mod_2)
confint(g_mod_2, oldNames=FALSE) # Computing profile confidence intervals ...
# Error in profile.merMod(object, which = parm, signames = oldNames, ...) : 
#   can't (yet) profile GLMMs with non-fixed scale parameters

dotchart(analysis_df_2$pers)
plot(g_mod_2)
qqnorm(resid(g_mod_2))
qqline(resid(g_mod_2)) # not great but ok 

hist(resid(g_mod_2),breaks=40)

analysis_df_2$resid = resid(mod2)

boxplot(analysis_df_2$resid~analysis_df_2$Colony)
ggplot(analysis_df_2, aes(x=Colony, y = resid)) + geom_boxplot() # + facet_grid(cols = vars(species))

vif(mod2) # Error in vif.merMod(mod1) : model contains fewer than 2 terms

r.squaredGLMM(mod2)

# trying out GAMM----

# Install and load the gamm4 package
library(gamm4)
# Fit the GAMM











