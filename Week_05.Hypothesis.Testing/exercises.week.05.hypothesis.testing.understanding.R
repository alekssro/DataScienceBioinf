library(tidyverse)


#### Hypothesis testing ####

# This week we will focus on a simpler dataset (chromosome 6)
# But we will only look at the european population
# Our variable of interest will be the minor allele frequency of all variable SNPs

# Our null hypothesis is: minor allele frequencies are the same everywhere on chromosome 6
# Our alternative hypothesis: some places have higher minor allele frequencies than expected (balancing selection)


# Load data and calculate maf

d = read_delim(file="1000g.allele.frequencies.tsv.gz",delim=c("\t"))

d = d %>%
  filter(population=="EUR") %>%
  mutate(maf=ifelse(frequency>0.5, 1-frequency, frequency)) %>%
  filter(maf > 0) %>%
  mutate(population=factor(population), 
         reference_allele=factor(reference_allele), 
         alternative_allele=factor(alternative_allele))

summary(d)

# For all snps calculate a "bin25k" variable based on position, each bin should be 25 kb


#Checkpoint
names(d)
#[1] "position"           "reference_allele"   "alternative_allele" "population"         "frequency"          "maf"                "bin25k"            
d %>% head(3)
#1    63979                C                  T        EUR    0.1044 0.1044 2
#2    63980                A                  G        EUR    0.1044 0.1044 2
#3    73938                A                  G        EUR    0.0010 0.0010 2

# For each bin calculate number of snps with maf > 0.2 (n20) , number of snps (n) , a test statistic (ts) = n20/n
# Also calculate the position of the bin midpoint (x)
# Call the binned results (1 row pr bin) for "br" (binned results)


# Checkpoint
names(br)
#[1] "bin25k" "n20"    "n"      "ts"     "x"     
br %>% head(5)
#1      2     0     3 0.0000000  68958.5
#2      3     0     2 0.0000000  88089.0
#3      4     0     7 0.0000000 111903.5
#4      5     1     7 0.1428571 146700.5
#5      6     2    25 0.0800000 162538.5

# Q: Plot this teststatistic along the chromosome and also visualise the number of snps in each bin



# Discuss: is this test statistic influenced by the number of snps?

# Discuss: the plot - do you see a signal?


#### Using probability ####


# Q: What is the observed p20 = n20/n for the entire chromosome

# Q: So of 1000 random snps how many do you expect to have maf > 0.2 ?

#### Imagine three bins... ####

# Imagine you have a bin with 99 SNPs, 30 of them have maf > 0.2
# Imagine another bin, with 100 SNPS, 22 of them have maf > 0.2
# Imagine another bin, with 1000 SNPS, 220 of them have maf > 0.2
# You have estimated an expected frequency of SNPs with maf > 0.20 to be 
p20 = 0.18 # 18 % this is NOT what you have calculated for the chrosome

#### Use simulation to estimate the null distribution and the p value ####

# Use 1000000 simulations to estimate the p value for the observation (30 / 99)
# H0: The frequency of maf>0.20 snps in this sample is p20 
# HA: The frequency of maf>0.20 snps in this sample is HIGHER than p20
# So is 30/99 significantly higher than p20 when it is based on 100 snps only?

set.seed(1)
x = sample(x = c(0,1), size = 99, prob = c(1-p20,p20), replace = TRUE)
# ^^^
# Hint!
# This is one way of getting 0 and 1s with frequency 1-p20 and p20.
x
y = sum(x)
y
y >= 30

set.seed(8603)
x = sample(x = c(0,1), size = 99, prob = c(1-p20,p20), replace = TRUE)
x
y = sum(x)
y
y >= 30

# Another way

set.seed(1) 
x = runif(n = 99, min = 0, max = 1)
# Hint: this is another way of getting random numbers between 0 and 1
x
x <= p20
y = sum(x <= p20)
y
y >= 30


set.seed(6853)
x = runif(n = 99, min = 0, max = 1)
x
x <= p20 
y = sum(x <= p20)
y
y >= 30


# Now do 1 million simulations and plot and save make a histogram of your simulated number of snps with maf>0.20
# Add the estimated pvalue to the plot as title ( plot1 + ggtitle("Hello world") )
# Save the plot

# Do this for both simulation strategies hinted above
# Q: Is one of them faster? Why?

# Method 1
set.seed(0) # This will make all all our results comparable


# Method 2
set.seed(0) # This will make all all our results comparable


#### Using binomial simulations to estimate the null distribution and the p value ####

n = 105
n20 = 29

set.seed(0)
x = rbinom(n = 1, size = n, prob = p20)
x
x >= n20
x = rbinom(n = 1000, size = n, prob = p20)
x
sum(x >= n20)
?rbinom

# Q: Use rbinom to simulate the first bin (33/99)  1 million times
# make a histogram of the number of snps in the bin with maf > 0.20 and add the estimated p value as title... and save it


# Compare the three different null distributions you have obtained.
# Q: How similar are the three different null distributions?

# Q: Select the method you like most from above and test all 3 bins using 1 million simulations

# In all cases we use the HA: 
# HA: The frequency of maf>0.20 snps in this sample is HIGHER than p20


# Compare the three results - what can you save about sample size and the null diributions?
# Q: What happens with increasing sample size? 

# Discuss : How many different null distributions exist?
# When we have different p20 and number of snps...

#
# 
#
# So instead of estimating p values from our own null distribution we
# may use a function that calculates the area under the curve...
#
#
#

#### Parametric test - use the binomial distribution as null distribution ####

n = 99
n20 = 30
x    = 1:n
y   = dbinom(x = x, size=n, prob=p20)
pd  = data.frame(x=x, y=y)
ggplot(pd, aes(x=x, y=y)) + geom_col(width = 0.5) + geom_vline(xintercept = n20, lty=2)

# Q: What is it we are calculating above?

# You can calculate the area under the curve as 
sum(dbinom(x = n20:n, size = n, prob=p20))

# Q: Is observing 3 out of 10 SNPs significant? (what is the pvalue?)


# Q: Why not? 3/10 >> p20 ?

# Q: What about 220/1000 ?

# Q: Does the test take sample size into account?

# Q: Also make a histogram using the binomial distribution of the 3 imagined bins
# Add dbinom calculated pvalues as titles


# Now you have seen 4 different ways of getting p values
# If you were to do the test for several thousand bins, which method would you use?

# Can you imagine cases were you would use another method?

# What are the downsides of taking a complex problem (distribution of maf) 
# and turning it into a simple problem (frequency of maf > 0.20) ?

# Next week we will test all bins on the chromosome
# And we will also try and use permutations to estimate the null distribution


