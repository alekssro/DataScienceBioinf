## A walkthrough anovas in R using the "ABD" Book Examples (chapter 15)

# 1-Way Anova with fixed effect (Book example 15.1) ####

circadian <- read.csv(url("http://www.zoology.ubc.ca/~schluter/WhitlockSchluter/wp-content/data/chapter15/chap15e1KneesWhoSayNight.csv"))
head(circadian) # a glimpse of the data
circadian$treatment <- factor(circadian$treatment, 
                              levels = c("control", "knee", "eyes")) ##to order treatments

stripchart(shift ~ treatment, data = circadian, method = "jitter", 
           vertical = TRUE, pch=20)
## Do a better graph :)

circadianAnova <- lm(shift ~ treatment, data = circadian) #Fitting the anova model with treatment as fixed effect
circadianAnova # Not very informative 
anova(circadianAnova) ## The most important way of summarizing an ANOVA should match the book table
summary(circadianAnova) ## to get the R^2 and other summaries
names(circadianAnova) ## inspecting the R fitted object 


plot(circadianAnova) ## diagnostic plot of the ANOVA 
hist(circadianAnova$residuals) ## visual check on residuals 
shapiro.test(circadianAnova$residuals)  ## normality test of residuals here we do not reject normality --> Good



#1 way ANOVA random (example 15.6) ####
walkingstick <- read.csv(url("http://www.zoology.ubc.ca/~schluter/WhitlockSchluter/wp-content/data/chapter15/chap15e6WalkingStickFemurs.csv"))
head(walkingstick)
library(nlme)
#Fitting the model 
walkingstickAnova <- lme(fixed = femurLength ~ 1, 
                         random = ~ 1|specimen, data = walkingstick)
walkingstickVarcomp <- VarCorr(walkingstickAnova) #n# Getting the variance components
walkingstickVarcomp ## the matrix of variance covariance fitted to the model 

varAmong  <- as.numeric( walkingstickVarcomp[1,1] ) ## Extracting the variance among groups
varWithin <- as.numeric( walkingstickVarcomp[2,1] ) ## Extracting the variance within groups
repeatability <- varAmong / (varAmong + varWithin)  ## Matching the book 
repeatability

