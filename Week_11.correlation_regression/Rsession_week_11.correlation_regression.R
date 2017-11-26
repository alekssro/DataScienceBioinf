library(tidyverse)

#### Distribution of levels of gene expression in Orangutans ####

# We will work (again) on the data from 
# Nguyen L-P, Galtier N, Nabholz B. 2015 Gene expression, chromosome heterogeneity and the fast-X effect
# in mammals. Biol. Lett.11 : 20150010.
# Data source:  
# http://dx.doi.org/10.1098/rsbl.2015.0010
# The data underlying this study are available on Dryad:
# doi:10.5061/dryad.qr20n


df = read.table("rsbl20150010supp1.csv", header = T, sep=",")

# As in week 3:
# Q: identify how each column of the data corresponds to the variables described in the methods section of the paper 
# If you already did this in week 3 - just be lazy ;)

summary(df)

#1.Expresion Level ---------> RPKM
#2.Expresion specificity ---> Tau
#3.GC content 3rd position -> GC3
#4.Species -----------------> Species
#5.dN/dS -------------------> dNdS
#6.Chromosome type ---------> chrMark

#### Overall questions ####

# Overall research question: is the dN/dS ratio correlated with expression level in Orangutan? Is it significant?



# How well does a linear model (linear regression) fit the data? What about residuals? Is it significant?

#### Correlation analysis ####

#?cor

# Since we expect the data to be very noise (week 3) we summarise it pr. chromosome.
# Q: Get median expression, median dnds and number of genes pr. chromosome pr. species


# Q: What is the correlation of dNdS and expression level for orangutan ? Present these results both visually as well as a number.
# Hint: ?cor


# Q: discuss and interpret the results. 



#### Testing the correlation ####

# Q: Is the observed correlation significant - how can you test it?

# Q: Try and calculate the correlation of 10000 permuted datasets - how many of these datasets have a more extreme correlation than the observed?
# Hint: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Using_a_permutation_test
# Compare your result with the output from cor.test ?


#### Bootstrapping confidence intervals ####

# Q: Bootstrap the confidence interval on the correlation coefficient? (Again look at wikipedia).
# Compare it with the output of cor.test


#### Regression analysis ####

# We use the lm function to build a linear model in R: 
df = data.frame(x=1:100, y = 5 + 0.05*(1:100) + rnorm(100))
ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point()
fit = lm(y ~ x, data = df) 
summary(fit)
names(fit)

# Q: What is the difference between lm and cor? 

# Q: Make a linear regression of dnds ~ expression 
# Q: inspect your regression object using names() and summary(): What coefficients are estimated for the regression? What is the R^2 of the model? 


# Q: Make a scatterplot of the median expression (x) vs. the dnds (y) and add a regression line? (Google for help).
# Also: highlight number of genes on chromosome (size)


#### Residuals ####

# Q: Make a residual plot (observed x vs. residuals)  - does it look ok?
# HINT: names(lm) and ?lm


#### Residuals when we KNOW the model is wrong ####

# Now we make a fake variable "fake" that has the values of expression^2 + some noise
pd$fake = 0.25 *  pd$expression^2 + rnorm(n = nrow(pd))

# Make a linear regression on fake ~ expression and plot the residuals plots again 
# Can you see that a linear fit is a bad predictor for these data?
# Also make a lm fit where you transform the expression to expression^2
# Can you see this is a better predictor?


#### Weighted regression ####

# Q: Do you see some data points with small sample size may be influencing the regression too much?
# Basically: do you see problems?
# Discuss: If chromosome with more genes have more precise dnds estimates - how can we improve the regression?
# Hint: ?lm

# Q: Perform and visualize a weighted regression and compare it with the normal unweighted regression
# HINT1: lm(x ~ y, weights = n, data = ...)
# HINT2: geom_smooth(data = ..., mapping = aes(..., weight = n), color="blue", se=F, method="lm")

# Discuss the results - what does it mean and what have you done?


#### Testing the regression  ####

# Read "The ANOVA approach" page 554 in the book and anova test the unweighted regression

# Q: Compare the anova result with the cor.test result.

myregression1 = lm(dnds ~ expression, data = pd)  # I called my summarised data for pd
summary(myregression1)
cor.test(pd$expression, pd$dnds, alternative = "two.sided", method = "pearson", conf.level = 0.95)



