## A walk through the naked mole rat example of the book ch 18. 

library(tidyverse)

library(abd) # The R package that comes along with the book


abdData(chapters=18)

data(MoleRats)
df = MoleRats %>% tbl_df()

summary(df)

ggplot(data = df, mapping = aes(x = ln.mass, y = ln.energy, symbol=caste)) + 
  geom_point()

ggplot(data = df, mapping = aes(x = ln.mass, y = ln.energy, colour=caste)) + 
  geom_point()+
  geom_smooth(method="lm")


fit.full = glm(formula = ln.energy ~ caste + ln.mass +  caste * ln.mass, data = df)  # Full model
summary(fit.full)
anova(fit.full, test = "F")  # Compare with Table 18.4.1

# Alternative way - compare models with and without the interaction variable
fit2 = glm(formula = ln.energy ~ caste + ln.mass, data = df)                    # caste + mass, No interaction
anova(fit2, fit.full, test = "F")

# Compare the two p values

#### Analyse if MASS is significant in the simple model?  ####

fit = glm(formula = ln.energy ~ ln.mass + caste, data = df)  # MASS + CASTE
anova(fit, test="F") # Table 18.4.2

###
### Notice the line in the output "Terms added sequenctially" !!!!!
###
### It looks like MASS is significant - but here we are comparing
### energy ~ constant vs
### energy ~ constant + mass
###
### What we should compare is
### energy ~ constant + caste vs.
### energy ~ constant + caste + mass
###

# We may save it by changing the order of predictors
fit2 = glm(formula = ln.energy ~ caste + ln.mass, data = df)  # CASTE + MASS
anova(fit2, test="F")

# Q: What is the difference between fit and fit2?
# Notice how the p values change

# To do it correctly: Compare models with and without the term
fit3 = glm(formula = ln.energy ~ caste, data = df)  # removing the mass from the model
anova(fit3,fit, test="F")

# Check that it is robust to order of variables
anova(fit3,fit2, test="F")

# What happens if we test two identical models but with different order of predictors?
anova(fit,fit2, test="F")

#
# install.packages("ggfortify")
library(ggfortify)

ggplot2::autoplot(fit)
