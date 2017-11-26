#### Week 04. PROBABLITY ####

#### R for data science exercises: Sampling from data & Sampling theoretical distributions  ####

# In week 1 and 2, you have spent some time learning how to manipulate data in the R tidyverse
# In week 3 you have explored what samplng variation means by sampling different types of distribution


# The goals of this week R session is to :
# Get familiarized with a few mathematical probability distributions and with the R functions that allow you to: 
# A. Sample random number from a probability distribution
# B. Calculate probabilties of an individual event or several events.
# These are the basic ingredient that we will then use in the coming weeks to do among other things hypothesis testing (Week 5)

# If you want to sample from a normal distribution you use rnorm()
#Q0: Look up the function dnorm() and qnorm() and examine what they are doing reagrding A and  B above . 
# 
# Q1:Make a graph of the probability density of the normal distrbution with mean =0 and sd=1.
# Q2:Sumperimpose and compare with a histogram of 10^3 draws from a normal.  
# (Think about what should be the scale of the Y axis to make a valid comparison)
# 
# Q3: Calculate the probability of drawing a random number from this normal distribution that exceeds 1, 2, 4 and 10 
# Q4: Calculate the probability of drawing a random number from this normal distribution that lies between 0, and 10. 

# 
# Q5 Approximate the probabilties above by calculating the proportions of random numbers that 
# exceed 1, 2, 4, 10  and and lie in the [0,10] in the histogram from above
# 
#### These are "pencil and paper question ####

#Q6 What is the probabilty that a sample of 10 individuals will contain a rare mutation in frequency 2.5% in the CEU population. 
# If we were planning a big sequencing studies, what number of individuals should be examined in order to ensure that we get at least one copy of the mutation with probabilty 50%, 90% and 99%?

# Q7 Developmental biologist have discovered that a set of 8 genes - Hox genes- control anterior-posterior identiy of regions of developping embryos. 
# We have recently discovered via GENOMIC & transcriptomic studies THAT THESE 8 genes are sequentially orderedon the chromosome in 
# the same ORDER AS THEY ARE EXPRESSED IN THE BODY FROM HEAD TO TAIL !
# HOW SURPRISING IS THAT FINDING? CAlculate the probability of that event in a model where order on the chromosome is independent of pattern of expression
