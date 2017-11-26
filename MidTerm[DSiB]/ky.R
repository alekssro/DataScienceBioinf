#####startup and load#####
library(tidyverse)
library(dplyr)

setwd("\\Users\\Ky\\Downloads")

d = read_delim(file="dgrp2.tgeno", col_names = T,delim=c(" "), na = "-")
head(d, 10)

#memory.limit(size=56000)

#>  Describe briefly the data? how many different chr?

####Q1####

#filter out chromosome 3L minus anything not strictly an SNP
bases <- c("A","T","C","G")
m <- d %>% 
  filter(chr == "3L", ref %in% bases & alt %in% bases) %>% 
  head(10^5)
#get positions
q1_first_pos <- min(m$pos)
q1_last_pos <- max(m$pos)

##Figure out genotyped %

#Add column with number of NA values in that row
m$num_na <- rowSums(is.na(m[,10:214]))
m[1:10,c(1,215)]

#Exclude rows based on num_na/204 (number of lines) 
#<75%

m_75 <- m %>% 
  filter(num_na/204 < 0.25)
#number left: 98135
nrow(m_75)

#<80%

m_80 <- m %>% 
  filter(num_na/204 < 0.20)
#number left: 97588
nrow(m_80)

m_90 <- m %>% 
  filter(num_na/204 < 0.1)
#number left: 94944
nrow(m_90)

m_95 <- m %>% 
  filter(num_na/204 < 0.05)
#number left: 90717
nrow(m_95)

####Q2####

##the distribution of SNP allele frequencies
ggplot(m_95, aes(x = refc))+
  geom_histogram()

ggplot(m_95, aes(x = altc))+
  geom_histogram()

#>maybe use maf?

##the distribution of coverage of SNPs
ggplot(m_95, aes(x  = cov))+
  geom_histogram(bindwidth = 1)


##the distribution of number of lines genotyped for each snp
ggplot(m_95, aes(x  = (204-num_na)))+
  geom_histogram()

##Is there a statistical association between the different types of SNPs and the allele frequency of SNPs?
## H0: there is no association

AT <- m_95 %>% 
  filter(ref == "A" & alt == "T")

AT <- AT[1:6519,7] 
#>`AC <- AC[1:nrow(AC),7]` or `AC <- AC$altc`-->(!! you get a integer vector)
#>independent of number of rows

AC <- m_95 %>% 
  filter(ref == "A" & alt == "T")

AC <- AC[1:nrow(AC),7]

GC <- m_95 %>% 
  filter(ref == "G" & alt == "C")

GC <- GC[1:nrow(GC),7]  

ATAC <- data.frame(ac = AC, at = AT)
ATGC <- data.frame(gc = GC, at = AT)
#>bellow line gives me error (different nº of rows)
#>create data.frame after bin calculation?

AT <- unlist(AT)
AT <- sort(AT, decreasing = FALSE)
at_bin <- tapply(AT, cut(AT, 10), median)

AC <- unlist(AC)
AC <- sort(AC, decreasing = FALSE)
ac_bin <- tapply(AC, cut(AC, 10), median)

GC <- unlist(GC)
GC <- sort(GC, decreasing = FALSE)
gc_bin <- tapply(GC, cut(GC, 10), median)

df1 <- data.frame(gc_bin,at_bin)
chisq.test(df1, correct = FALSE)
##p-value = AC/AT = 1, AT/GC = 1 =>no bias

#>another way to check statistical association visually could be to plot
#>AT vs GC counts by bins and you can see it fits to a straight line
# ggplot(data = df1, mapping = aes(x = gc_bin, y = at_bin)) +
#     geom_point() +  #Each point represent counts in range
#     geom_smooth(method = lm)

####Q3####
# If coverage was homogenous throughout the genome,
# (by that we mean that on average the coverage is the same for any given position),
# what probability distribution is expected to capture well the coverage?
#Poisson around mean(cov)

####Q4####
tbl2 <- table(rpois(nrow(m_95),mean(m_95$cov)), m_95$cov)
chisq.test(tbl2, correct = FALSE)
# p-value varies, but mostly a good fit

#>one way to take into account variations on comparison with random generated distributions
#>is to repeat the test multiple times and take the mean p-value for this tests
#>Something similar to:
# pvals <- replicate(1000, {
#     chisq.test(data.frame(x = rpois(nrow(m_95), lambda = mean(m_95$cov)), y = m_95$cov))$p.value
# })
# 
# mean(pvals)

# poistest <- rpois(nrow(m_95), mean(m_95$cov))
# 
# ggplot(melt(poistest), aes(value)) +
#   geom_histogram(binwidth = 1)

####Q5####
m_kb <- d %>% 
  filter(chr == "3L", ref %in% bases & alt %in% bases) %>% 
  slice(c((3*10^5):(3*10^5+10^5-1)))

ggplot(m_kb, aes(x = refc))+
  geom_histogram()

#>Doesn't it have to be reference allele frequency by its position?
#>scattterplot: x = pos | y = ref frequency of the SNP

LD <- c(0,0)
LD[1:100] <- 0
dist <- c(0,0)
dist[1:100] <- 0

for(i in 1:100) {
  tbl3 <- c(0, 0, 0, 0)
  
  for (n in 10:214) {
    if (!is.na(m_kb[i, n]) & !is.na(m_kb[i+1, n]) ){
      
      if ((m_kb[i, n] == 0) & (m_kb[i+1, n] == 0)){
        tbl3[4] <- tbl3[4] +  1
      }
      if ((m_kb[i, n] == 2) & (m_kb[i+1, n] == 0)){
        tbl3[2] <- tbl3[2] +  1
      }
      if ((m_kb[i, n] == 0) & (m_kb[i+1, n] == 2)){
        tbl3[3] <- tbl3[3] +  1
      }
      if ((m_kb[i, n] == 2) & (m_kb[i+1, n] == 2)){
        tbl3[1] <- tbl3[1] +  1
      }
    }
  }
  N <- sum(tbl3)
  LD[i] <- (tbl3[1]-((tbl3[1]+tbl3[2])*(tbl3[1]+tbl3[3])))/N
  dist[i] <- sqrt((m_kb[i,2]-m_kb[i+1,2])^2)    #>abs() to save operations (you should use unlist() afterwards)

  # N <- sum(tbl3)
  # tbl3 <- tbl3 / N
  # LD[i] <- tbl3[1] - (tbl3[1] + tbl3[2]) * (tbl3[1] + tbl3[3])
  # dist[i] <- sqrt((m_kb[i,2]-m_kb[i+1,2])^2)
}

LD_df <- data.frame(LD = LD, distance = dist)
#>LD should be: -1 <= LD <= 1
#>I think the problem is that you have to divide by N to obtain the freqs before calculating LD
#>I have commented what I mean

ggplot(LD_df,aes(x = distance, y = LD))+geom_point()
#looks like fig 1, but upside down? I only chose the first 100 pairs because of computing time, 100 pairs took ages
#>try representing abs(LD) or LD^2 to avoid negative values?:
# ggplot(LD_df,aes(x = distance, y = abs(LD)))+geom_point()

m_kb[10:20,c(3,5)]
