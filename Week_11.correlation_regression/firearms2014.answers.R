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
df = read_delim("us_states_guns_and_deaths.csv",delim = ";" )
summary(df)

#!begin

# Working on overall death rate

library(ggrepel)
source("nature_theme.R")

cor(df$gun_ownership, df$rate_all, method = "pearson")

fit = lm(rate_all ~ gun_ownership, data = df) 
summary(fit)

plot1 = ggplot(data = df, mapping = aes(x = gun_ownership, y = rate_all, color=region, label=state, group=1)) +
  geom_point(aes(size=n)) + 
  geom_smooth(method="lm")+
  geom_text_repel()+
  ggtitle("Gun ownership vs. death rate in US states")

plot1 = plot1 + nature_theme()

plot(plot1)

# Working on non suicides only

cor(df$gun_ownership, df$rate, method = "pearson")

fit = lm(rate ~ gun_ownership, data = df) 
summary(fit)

plot1 = ggplot(data = df, mapping = aes(x = gun_ownership, y = rate_homicide, color=region, label=state, group=1)) +
  geom_point(aes(size=n)) + 
  geom_smooth(method="lm")

plot(plot1)



#!end

# USA vs. rest of the world
df = read_csv("country-gun-stats.csv")
summary(df)

# Get industrialized countries only
df %>% filter(GDP > 20000)

#!begin

pd = df %>% filter(GDP > 20000)

cor(pd$homicides, pd$guns_pr_capita)

fit = lm(homicides ~ guns_pr_capita, data = pd) 
summary(fit)

ggplot(data = pd, mapping = aes(x = guns_pr_capita, y = homicides, label=Country)) +
  geom_smooth(method="lm", se=TRUE, lty="dashed", color="#808080", fill="#d0d0d0")+
  geom_point(size=3)+
  geom_text_repel()+
  ggtitle("Gun ownership vs. homicide rate in industrialized countries")

plot1 = plot1 + nature_theme()
plot(plot1)
ggsave(plot = plot1, filename = "plot.png")


#!end




