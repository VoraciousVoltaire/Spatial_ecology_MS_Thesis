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

# Some exploratory data analysis first----

df <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/csv/Final_hypothesis_1_analysis_df.csv")
df # has colony name, individual number, individual ID, pers, number of tracks per individual to 
# calculate its KDE, vals, nas and percent_na. 

# 1. Checking for outliers----
par(mfrow = c(1,1))
dotchart(df$pers) # a left tail, a few outliers

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
mod1 = glmer(pers ~ Colony + (1 | Individ_id), data = df, family = Gamma(link = "identity"), control = control_params)
AIC(mod1) # 599.9103
BIC(mod1) # 647.4424
summary(mod1)


# confint(mod1, oldNames=FALSE)

plot(mod1)
qqnorm(resid(mod1))
qqline(resid(mod1)) # not great but ok 


hist(resid(mod1),breaks=20)
shapiro.test(resid(mod1))

df$resid = resid(mod1)

boxplot(df$resid~df$Colony)
ggplot(df, aes(x=Colony, y = resid)) + geom_boxplot() # + facet_grid(cols = vars(species))

vif(mod1) # Error in vif.merMod(mod1) : model contains fewer than 2 terms

r.squaredGLMM(mod1)


