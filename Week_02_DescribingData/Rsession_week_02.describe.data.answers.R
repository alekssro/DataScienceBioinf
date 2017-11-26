#!
#### Week 02 ####

#### R for data science exercises: describing data ####

# The next part of the exercises will learn you to manipulate data in the R tidyverse
# Like last week you will work through online tutorials that will introduce you to data manipulation
# We advice you to make a new file (R script) and use that for doing the tutorial.
# When done - save the R script with an appropriate filename and in an appropriate location
#
#
# The tutorial is from the free book "R for data science" written by the R Overlord Hadley Wickham
# URL: http://r4ds.had.co.nz/index.html
#
#
# Read Welcome
# Read and do 4. and 5. 
# For general introduction to Rstudio: read 6.
#
# If you need help have a look at the different cheatsheets
# And if you are stuck - use the online forum!
#


#### Use your new knowledge to Work on real data ####

library(tidyverse)

# By using the knowledge from above we will now work on a dataset
# The data is the SNP allele frequencies from the 1000 genomes project
# from chromosome 6 (you can see the README file to see how I downloaded and formatted this dataset)

# We can read zipped files (saves transfer time and space on disk)
d = read_delim(file="1000g.allele.frequencies.tsv.gz",delim=c("\t"))

#### Counting SNPs ####

# Q: What is the distribution of the variable "frequency"?

ggplot(data = d) + geom_histogram(aes(x=frequency)) #!

# Q: Is it different for the different populations?

# Notice the few SNPs with high frequency (> 0.5)
# Q: How many SNPs pr population have frequency > 0.5?
# Q: What does this mean? (a frequency > 0.5)

d %>% #! 
  filter(frequency > 0.5) %>%  #!
  group_by(population) %>%     #!
  summarize(n=n())             #!

# We calculate minor allele frequency
d = d %>% mutate(maf=ifelse(frequency>0.5, 1-frequency, frequency))

# Q: How many SNPs have maf > 0 (are polymorphic) for each population?
d %>% filter(maf > 0) %>% group_by(population) %>% summarize(n=n())  #!

# Q: How many SNPs have maf > 0.05 (common polymorphism) for each population?
d %>% filter(maf > 0.05) %>% group_by(population) %>% summarize(n=n())  #!

#### Interquartiles ####

# Q: Calculate the mean, sd and median and the 0.25, 0.5 and 0.75 quantiles of the maf for each population
# HINT: ?quantile

#!begin

d %>% 
  filter(maf > 0) %>%
  group_by(population) %>%
  summarise(mean_maf = mean(maf),
            sd_maf=sd(maf),
            median_maf=median(maf), 
            q25 = quantile(maf, probs=0.25),
            q50 = quantile(maf, probs=0.50),
            q75 = quantile(maf, probs=0.75))

#!end

# Q: How many SNPs are fixed for the reference allele in AFR but fixed for the alternative allele in EUR and EAS respectively?
# HINT: the original frequency is of the alternative allele

# A fixed allele means all individuals have the same allele (e.g. G/G = homozygote)
# For the computer scientists: you inherit an allele from you mother and your father.
# So basically your genotype can be:
# REF/REF = homozygotic reference allele
# REF/ ALT = heterozygote
# ALT/ALT = homozygotic alternative allele
# If the frequency==0.0 it means the reference allele is fixed (and 1.0 corresponds to the alternative allele being fixed)

#!begin

d %>% 
  group_by(population) %>% 
  summarize(fixed_ref = sum(frequency == 0.0),
            fixed_alt = sum(frequency == 1.0))
  

#!end

#### Cumulative  frequencies ####

# Q: Make a cumulative frequency plot of the minor allele frequency for the different populations (color)
# Zoom in to maf 0-0.05
# HINT: ?stat_ecdf()

#!begin

ggplot(d, aes(x=maf, colour = population)) + stat_ecdf(pad=F) + xlim(c(0,0.05)) #!

#!end

# Q: Explain the plot - what do you see?
# Q: Which population has the most "rare" alleles?

#! EAS followed by EUR and then Africa.


#### Looking for interesting patterns ####

# Now we will try and see if there is somewhere on the chromosome with a funny allele frequency structure.
# Basically we ask - can you find something interesting?

# Q: Try and plot a scatterplot of x=position and y=maf 

#ggplot(d) + geom_point(aes(x=position, y=maf))                #!

# Sorry - it will likely crash your machine
# HINT: ... %>% sample_n(100) will select 100 random rows from your data

# Q: Are you able to detect a pattern or not?

# Q: What kind of pattern do you see?

# Q: Instead of only plotting some of the data we want to summarise the data.
# Add a "bin" variable that is the binned position, in bins of size 25000
# and summarise the mean and sd pr. bin (of the minor allele frequency)
# Also get the number of observations for each bin
# HINT: 1223 %/% 100 = 12,  12 * 100.00  + 50 = 1250 = all values 1200-1299 will become 1250 
# ? %/% and ?%% - I think you used them in the R for datascience exercises as well.

#!begin

pd = d %>%                                                     #!
  mutate(bin = ((position %/% 25000) * 25000) + 12500) %>%     #!
  group_by(population, bin) %>%                                #!
  summarize(maf_mean  = mean(maf),                             #!
            maf_sd    = sd(maf),                               #!
            maf_count = n())                                   #!
pd #!

#!end

# Q: Plot the number of SNPs pr. bin - a subplot for each population - do you see something?
ggplot(pd,aes(x=bin, y=maf_count)) + geom_point() + facet_wrap( ~population, ncol=1) #!

# Q: What do you see - is there a pattern in snp density and/or allele frequency along the chromosome?

# Q: Plot the mean MAF pr. bin - a subplot for each population (facet) and also illustrate first sd then number of SNPs in bin (color)

#!begin

ggplot(pd,aes(x=bin, y=maf_mean, color=maf_sd)) + geom_point() + facet_wrap( ~population, ncol=1) #!
plot1 = ggplot(pd,aes(x=bin, y=maf_mean, color=maf_count)) + geom_point() + facet_wrap( ~population, ncol=1)      #!
plot1 #!

#!end

# Q: Do you see a peak in maf somewhere?

#! Yes - around 30 mb

# Q: Use the ggplot function coordinate_cartesian(xlim=c(,)) to zoom in on this region
plot1 + coord_cartesian(xlim=c(25000000,35000000)) #!

# Q: What kind of natural selection do you think is present here? (purifying or balancing selection?)

#! Balancing selection

# Q: What kind of genes do you think are present here?
# Google or use a genome browser to help you

#! Immune genes - the HLA



#### Transitions and transvertions ####

# A-G  and C-T are transitions (purine - purine and pyrimidine - pyrimidine)
# A-T and A-C are transversions, i.e. when REF allele is purine and ALT allele is a pyrimidine or the other way around

purines     = c("A", "G")
pyrimidines = c("C", "T")

# Q: If all mutations were random: what is the expected transition:transversion ratio?

# Q: What is the observed transition:transversion ratio for each population?

#!begin

d %>% filter(maf > 0) %>% 
  mutate(type=ifelse( (reference_allele %in% purines & alternative_allele %in% purines) 
                      |      
                        (reference_allele %in% pyrimidines & alternative_allele %in% pyrimidines), 
                      "transition", "transversion")) %>%
  group_by(population,type) %>% 
  summarize(n=n()) %>% 
  spread(type,n) %>%
  mutate(ratio=transition/transversion)

#!end

# End of exercise
