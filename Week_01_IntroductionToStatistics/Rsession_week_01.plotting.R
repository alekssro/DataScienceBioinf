#### Introduction ####

#
# Datascience in Bioinformatics (an extended version of an earlier course called LEARNING FROM GENOME DATA I
#
# The questions and exercises are identified as # Q:

#### This is a section header ####

# Section headers make it a lot easier to navigate your script
# Please check the drop down menu at the bottom of this window.

#### R for data science exercises: Data visualization ####

# First you should work through two online tutorials that will introduce you to ggplot2
# We advice you to make a new file (R script) and use that for doing the tutorial.
# When done - save the R script with an appropriate filename and in an appropriate location
#
#
# The tutorial is from the free book "R for data science" written by the R Overlord Hadley Wickham
# URL: http://r4ds.had.co.nz/index.html
#
#
# Read Welcome
# Read 1. Introduction
# Read 2. Introduction
# Go through 3. Data visualization (this takes some time but will introduce you to ggplot2)
#
# If you need help have a look at the different cheatsheets
# And if you are stuck - use the online forum!
#


#### Use you new knowledge to Work on real data ####7
# install.packages(c("tidyverse", "ggplot2"))

library(tidyverse)
library(dplyr)

# By using the knowledge from above we will now work on a dataset
# Read the data 

microbes = read_tsv(file="taxontable.tsv") # read data
# Notice the warning - let's tell it that it may read more rows to guess the correct type of each row
microbes = read_tsv(file="taxontable.tsv", guess_max=5000) # read data

microbes
dim(microbes)
names(microbes)
summary(microbes)

#### Variation of genome size and number of genes ####

# # Hint for keeping and saving a plot
# plotdata = data.frame(x=rnorm(100), y=rnorm(100), type=sample(x = 1:3, size = 100, replace = T))
# plot1 = ggplot(data = plotdata) + geom_point(mapping = aes(x = x, y = y, color=type))
# plot(plot1)
# ggsave(filename = "test.png", plot = plot1)


# Q: Make a scatterplot of Genome.Size(x) vs. Gene.Count(y)
plot1 <- (
    ggplot(data = microbes) + 
    geom_point(mapping = aes(x = Genome.Size, y = Gene.Count, color = "red"))
)

plot(plot1)
ggsave(filename = "RsessionWeek01/Q1.png", plot = plot1)

# Q: Make a scatterplot of Genome.Size(x) vs. Gene.Count(y) on log-log scale
# Hint: ?scale_x_log10() (read examples)
plot2 <- (
    ggplot(data = microbes) + 
        geom_point(mapping = aes(x = Genome.Size, y = Gene.Count, color = "red")) +
        scale_x_log10() +
        scale_y_log10()
)


plot(plot2)
ggsave(filename = "RsessionWeek01/Q2.png", plot = plot2)

# Q: Make a log-log scatterplot so you compare the different Domains (different colors)
plot3 <- (
    ggplot(data = microbes) +
        geom_point(mapping = aes(x = Genome.Size, y = Gene.Count, color = Domain)) +
        scale_x_log10() +
        scale_y_log10()
)

plot(plot3)
ggsave(filename = "RsessionWeek01/Q3.png", plot = plot3)

# Q: Make a facetted log-log scatterplot so you compare the different Domains 
# HINT: ..facet_wrap()
plot4 <- (
    ggplot(data = microbes) +
        geom_point(mapping = aes(x = Genome.Size, y = Gene.Count, color = Domain)) +
        scale_x_log10("Gene Size") +
        scale_y_log10("Gene Count") +
        facet_wrap(~Domain, nrow = 1)
)

print(plot4)
ggsave(filename = "RsessionWeek01/Q4.png", plot = plot4)

# Q: What can you say about the relationship between genome size and gene count.
# You will present this and discuss it at the exercises

#### The relationship between genome size and gene density ####

# We define a new variable called "Gene.density" in the dataset (number of genes per Mb)
microbes = mutate(microbes, Gene.Density = Gene.Count/10^6)

# Q: How is gene density distributed ?  
# Hint: +geom_histogram()
plot5 <- (
    ggplot(data = microbes, aes(x = Gene.Density, fill = Domain)) +
        scale_x_log10("Gene Density") +
        geom_histogram(bins = 200)
)

print(plot5)
ggsave(filename = "RsessionWeek01/Q5.png", plot = plot5)

# Q: Make the histogram but with subplots for each Domain (facets)
plot6 <- (
    ggplot(data = microbes, aes(x = Gene.Density, fill = Domain)) +
        scale_x_log10("Gene Density") +
        facet_wrap(~Domain, nrow = 3) +
        geom_histogram(bins = 200)
)

print(plot6)
ggsave(filename = "RsessionWeek01/Q6.png", plot = plot6)

# Q: Plot the gene density as a function of genome size on a log-log scale
plot7 <- (
    ggplot(data = microbes) +
        geom_point(mapping = aes(x = Genome.Size, y = Gene.Density)) +
        scale_x_log10("Gene Size") +
        scale_y_log10("Gene.Density")
)

print(plot7)
ggsave("RsessionWeek01/Q7", plot = plot7)

# Q: Produce the same figure, but each domain in a different color
plot8 <- (
    ggplot(data = microbes) +
        geom_point(mapping = aes(x = Genome.Size, y = Gene.Density, color = Domain)) +
        scale_x_log10("Gene Size") +
        scale_y_log10("Gene Count") 
)

print(plot8)
ggsave("RsessionWeek01/Q8", plot = plot8)

# Q: Produce the same figure, but each domain in a different subplot (facets)
plot9 <- (
    ggplot(data = microbes, aes(x = Genome.Size, y = Gene.Density, color = Domain)) +
        geom_point() +
        geom_smooth() +
        scale_x_log10("Gene Size") +
        scale_y_log10("Gene Count") +
        facet_wrap(~Domain, nrow = 1)
)

print(plot9)
ggsave("RsessionWeek01/Q9", plot = plot9)
# Q: Present your results at the exercises.

# You are now done: https://goo.gl/3XLFxY

