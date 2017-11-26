library(tidyverse)


#### Hypothesis testing ####

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
d = d %>%                                                     
    mutate(bin25k = position %/% 25000) %>%
    ungroup()


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



# Q: What is the observed p20 = n20/n for the entire chromosome


#### Continuing from last week ####

# Q: select the best/easiest way of tesing if the observed proportion of SNPs with maf > 0.2 in the bin is higher than expected



# Q: List the 10 most significant bins

# What is the lowest pvalue?

# How many of your bins have p < 0.001?

# If all bins followed H0: how many would you expect to have p > 0.001

# Q: Plot the p value as function of bin position

# Q: It is really difficult to see the small pvalues - try to mutate a new pvalue2 = -log10(pvalue) and plot it
# This is called a Manhattan plot - strong signals will be skycrapers of significance


# Q: Do you see a skyscraper?

#### Dividing the chromosome into bins with same number of snps ####

# Q: Do the same analysis as before but using bins of size 500 snps instead
# Basically we just want the Manhattan plot
# HINT: 
d %>% 
  arrange(position) %>% 
  mutate(SNPnumber = row_number())



#### Testing a specific hypothesis ####

# Assume that I speculate that the overall frequency of SNPs with maf > 0.05 is really 50% in humans.
# NOTE: maf > 0.05 - I call these "high maf snps"
# Can you test if some bins of size 100 kb have significantly more high maf SNPs?
# Visualize the test results so it is easy to see where significant bins
# HINT: Make manhattan plots but also try and remove all bins with less than 50% of the snps having maf > 0.05 


#### Bonus question ####

# Can you repeat the final test for all three populations?
# Do you see a difference between the populations?

