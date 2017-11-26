#### Week 04. PROBABLITY ANSWERS ####

#### R for data science exercises: Sampling from data & Sampling theoretical distributions  ####
library(ggplot2)
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
# Q1
df <- data.frame(x=rnorm(100, 0, 1))
head(df)
ggplot(df, aes(x)) +
  geom_histogram(binwidth = 0.2) +
  coord_cartesian( xlim= c(-3,3, by = 0.5)) +
  labs(x ="Random observation from a Normal", y = "Observed Counts per class")

ggplot(df, aes(x)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.2) +     #NOW A DENSITY                   
  coord_cartesian( xlim= c(-3,3, by = 0.5)) +
  labs(x ="Random obs from a Normal", y = "NOT COUNTS BUT A DENSITY = PROB (class)*binwidth")


# overlay histogram and "THE" normal density and a "fitted" Normal
ggplot(df, aes(x)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.2) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(df$x), sd = sd(df$x)), ###THIS IS A NORMAL FITTING THE DATA by using data to get mean and sd
                lwd = 2, 
                col = 'red')+
  stat_function(fun = dnorm, 
              args = list(mean = 0, sd = 1), ### THIS IS THE CANONICAL N(0,1) normal 
              lwd = 2, 
              col = 'blue')

# (Think about what should be the scale of the Y axis to make a valid comparison)
# 
# Q3: Calculate the probability of drawing a random number from this normal distribution that exceeds 1, 2, 4 and 10 
sum(df$x>2)/length(df$x)
sum(df$x>2)/length(df$x)
1- pnorm(q =  0,mean = 0,sd = 1)

# Just an old fashion plot of the cumulative distribution function 
Xs=seq(-5,5,by = 0.05)
Y1s=pnorm(q =  Xs,mean = 0,sd = 1) ## NOTE we again give a vector "Xs" to pnorm --> returns a vector
Y2s=pnorm(q =  Xs,mean = 0,sd = .5) ## a normal density with less variance --> faster transition from 0 to 1

plot(Xs,Y1s, type="l", xlab="t", ylab="P(X<= t)", col="cornflowerblue", lwd=3)
lines(Xs,Y2s) # LESS VARIANCE --> faster transition from 0 to 1

# Q4: Calculate the probability of drawing a random number from this normal distribution that lies between 0, and 10. 
probO_10=pnorm(q =  10,mean = 0,sd = 1) - pnorm(q =  0,mean = 0,sd = 1)

# 
# Q5 Approximate the probabilties above by calculating the proportions of random numbers that 
# exceed 1, 2, 4, 10  and and lie in the [0,10] in the histogram from above
# 
#Two possible ways ... 

sum(df$x>1)/length(df$x)
length(which(df$x>1))/length(df$x)

#Answer Compare these two possibilities using the logical comparison ">" as a function or using which 
df$x>1
length(df$x>1)
sum(df$x>1)
which(df$x>1)


sum(df$x>2)/length(df$x)

sum((df$x>=0)*(df$x <= 10))/length(df$x)

### NB now you can also quantif the uncertainty around these probabilities estimated by proportions of numbers in a certain probability range:
# for a given range that cover a true probability pR the number of random simulations that will fall in that range out of S simulations is 
 # >>>>> binomial with S trial and poproability pR!! ---> you can get "almost" automatically the precision on your estimation of probabilities from the binomial variance




##Mini example on Medelian segreagation in a set of recombinant inbred lines
# NG is the number of green ancestral bacground alleles
MyNG=rbinom(n=1000,size=100,prob=0.5)
hist((MyNG), xlab="Number of GREEN out of 100", col="green", main="Mendelian??")
rug(65, lwd=2, col="red") ##superimposing an aobservation for a SNP where 65 individuals out of 100 are bearing the "green" ancestral background




