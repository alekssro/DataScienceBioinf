library(tidyverse)

df = data.frame(x=1:100, y = 5 + 0.05*(1:100) + rnorm(100))

# cor() is used for calculating correlations?
# cor.test() for testing
?cor
cor(x = df$x, y = df$y, method = "pearson")

# We use the lm function to build a linear model in R: 
?lm

fit = lm(y ~ x, data = df) 
summary(fit)
names(fit)

ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point() + geom_smooth(method = "lm")

# US states only
df = read_delim("firearms2014.csv",delim = ";" )
summary(df)


# USA vs. rest of the world
df = read_csv("country-gun-stats.csv")

# Get industrialized countries only
df %>% filter(GDP > 20000)





