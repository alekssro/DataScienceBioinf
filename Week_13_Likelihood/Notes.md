## Class notes 28/11/2017
### Likelihood example

For the example about father and mother mutation rate passed to the offspring: Poisson regression model.

Also, we can assume that the response variable is normally distributed so we could use Generalized Linear Models with Likelihood as the way of fitting the model.  

glm(data, response_variable ~ covariable, family = "poisson")

anova(glm, test = "LRT") # Likelihood ratio test (Am I fitting the data better?) 

Yes, adding the slope gives a (massive) increase of the fit. Very likely.

