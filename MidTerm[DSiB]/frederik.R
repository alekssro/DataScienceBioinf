## Q0. Import the data

# (space delimited, 0 = ref allele, 2 = alt allele (not necessarily minor), - = missing)
data <- read.delim("dgrp2.tgeno", header = T, sep = " ", nrows = -100, na = '-')
#>read_delim() is faster than read.delim
data <- read_delim(file="dgrp2.tgeno", col_names = T,delim=c(" "), na = "-")

# Data contains gene changes recorded for 205 different fruit flies.
# Each datapoint contains a chromosome, the position on the chromosome, type of change, the reference gene expression and the alternative, the count of each, coverage and qual, and for each of the 205 flies whether they express the reference or alternative for that particular change.
# The datapoints all comes from the chromosomes 2L 2R 3L 3R 4 X. 



## Q1. Extract the first 100,000 SNPs located on chromosome 3L, 
# Extract SNPS, i.e. all changes contain only a single nucleotide, and change chr ind id to characters instead of factors
snps <- data %>%
  filter(nchar(as.character(alt)) == 1 &
         nchar(as.character(ref)) == 1) %>%
  transform(chr = as.character(chr),
            id = as.character(id))

# Drop the now unused factor levels for ref and alt
snps[,'ref'] = droplevels(snps[,'ref'])
snps[,'alt'] = droplevels(snps[,'alt'])

# Only look at the first 100000 SNPs which are located on the left arm of chromosome 3 and are actually polymorphic
snps3L <- snps %>%
  filter(chr == '3L', refc > 0, altc > 0) %>%
  head(100000)


# Locate the genomic position of the first and last SNP
min(snps3L$pos)
max(snps3L$pos)
# filter out all SNPs for which less than 75%, 80%, 90%, 95%  of the individuals could be genotyped, report how many SNPs are left based on the filter you apply

proportion_sequenced <- function(dataset, prop) {
  dataset %>%
    filter((refc + altc) / 205 > prop)
} 

nrow(proportion_sequenced(snps3L, 0.75)) # 98195
nrow(proportion_sequenced(snps3L, 0.80)) # 97604
nrow(proportion_sequenced(snps3L, 0.90)) # 95625
nrow(proportion_sequenced(snps3L, 0.95)) # 90725

genotyped95 <- proportion_sequenced(snps3L, 0.95)
# Add maf as a feature for the dataset
genotyped95 <- genotyped95 %>%
  mutate(maf = ifelse(altc < refc, altc / (refc+altc), refc / (refc+altc)))


## Q2. Graph the distribution of SNP allele frequencies, graph the distribution of coverage of SNPs ("cov") and the distribution of number of lines genotyped for each snp
ggplot(genotyped95, aes(x = maf)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(title = 'Distribution of maf')

# distribution of coverage
ggplot(genotyped95, aes(x = cov)) +
  geom_histogram(bins = 31, aes(y = ..density..)) +
  theme_minimal() +
  labs(x = 'coverage', y = 'proportion', title = 'Histogram of the Distribution of Coverage')

# distribution of number of lines
ggplot(genotyped95, aes(x = refc + altc)) +
  geom_bar() +
  theme_minimal() +
  labs(x = 'Number Sequenced', title = 'Distribution of Sequenced Flies Where More Than 95 is Sequenced')

# Association between strong and weak alleles
# Number of Weak and Strong alleles in each 5% interval of maf
weak_strong = genotyped95 %>%
  group_by(range = as.integer(maf * 100) %/% 5) %>%
  summarise(ats = sum(ref %in% c('A', 'T') & alt %in% c('A', 'T')),  
            cgs = sum(ref %in% c('G', 'C') & alt %in% c('G','C')))

# collaps the last category with only those genes with maf exactly 50% with the category for containing 45-50
weak_strong[10, 2:3] = weak_strong[10, 2:3] + weak_strong[11, 2:3]
weak_strong = weak_strong[1:10, ] # drop the last column
weak_strong$range = weak_strong$range * 5

ggplot(weak_strong, aes(x = range)) +
  geom_col(aes(y = ats), fill = 'blue') +
  geom_col(aes(y = cgs), fill = 'red') +
  theme_minimal() +
  labs(x = 'maf', y = 'count', title = 'Maf of weak vs. strong alleles', subtitle = 'Blue: Weak Alles, Red: Strong Alleles')

chisq_test = chisq.test(as.matrix(x = weak_strong[, 2:3]))
chisq_test 
abs(chisq_test$observed - chisq_test$expected)
#>another way to check statistical association visually could be to plot
#>weak vs strong counts by bins and you can see it fits to a straight line:
# ggplot(data = weak_strong[, 2:3], mapping = aes(x = ats, y = cgs)) +
#     geom_point() +  #Each point represent counts in range
#     geom_smooth(method = lm)

# Since our null hypothesis - that weak and strong alleles are independent of maf - has a p-value of 0.024 we can reject it alpha = 0.05

## Q3. If coverage was homogenous throughout the genome (by that we mean that on average the coverage is the same for any given posiiton),
# - what probability distribution is expected to capture well the coverage ? 

# This is best modelled using a poison distribution to test whether the 'succes' of obtaining a low 'low coverage' datapoint is independent of the space of the gene


## Q4. Make goodness of fit test of the coverage data: is the theoretical distribution proposed in Q3 a good fit for the data?

# Divide the data into 10 parts each containing the same number of genes
positions = c(0, quantile(x = genotyped95$pos, probs = seq(from = 0.1, to = 1, by = 0.1)))
low_cov = quantile(x = genotyped95$cov, probs = 0.05) # The coverage we will define as being 'low' coverage
low_cov_prop = mean(genotyped95$cov <= low_cov) # proportion being low-coverage (not exactly 5% since we have repeated values)

# The expected number of low coverage genes in each bucket. It is the same since under the null they have an equal probability across the gene and they are the same size
expected = rep(nrow(genotyped95) * low_cov_prop / 10, 10) 
observed = rep(-7, 10)
for (i in seq(1, 10)) {
  # The number of genes between position i and i+1 of which cov is 'low'
  observed[i] = sum(genotyped95[between(genotyped95$pos, positions[i], positions[i+1]), ]$cov <= low_cov)
}

chisq_statistic = sum((observed - expected)^2 / expected)
pchisq(chisq_statistic, df = 10 - 1 - 1, lower.tail = F) # 10 categories and one estimated parameter

# Our p value is ridiculously small, so the amount of 'low coverage' is not evenly distributed across the gene

# The number of 'low coverage' genes in each of the equal-sized buckets over the gene, the line representing the expected number of genes in the bucket
ggplot(data.frame(x = 1:10, y = observed), aes(x, y)) +
  geom_col() +
  geom_hline(yintercept = expected[1], color = 'red') +
  theme_minimal() + 
  labs(x = 'Bins of each 10% of the data', y = '#low coverage genes', title = 'Distribution of low-coverage genes throughout the geneome', subtitle = 'The line represents the expected number of low-coverage genes in each bin')

# How the equal-sized buckets are distributed across the gene
ggplot(genotyped95, aes(pos)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = positions, color = 'red') +
  theme_minimal()

## Q5. Pick yourself a small (100 kb) region on chromosome 3L


smallRegion <- genotyped95 %>%
  filter(between(pos, 450000, 450000 + 100000))

##  Graph the frequency of "0" alleles along that region
# Frequency of the reference gene along the part of the chromosome
ggplot(smallRegion, aes(pos, refc / (refc+altc))) +
  geom_smooth() +
  theme_minimal() +
  labs(y = 'Proportion expressing the reference', title = 'The proportion expressing the reference throughout the genome')


 ## Calculate the statistical association between allelic states of each pair of neighboring SNPS (this association is also called LD)

# Translate lines into matrix
# factor gets translated to int: - becomes 1, 2 becomes 3, 0 becomes 2
lines = as.matrix(sapply(smallRegion[, seq(10, ncol(smallRegion)-1)], as.integer))
lines[lines == 1] = NA
lines[lines == 2] = 0
lines[lines == 3] = 2
#>I think there is a problem with this. Wouldn't it be like this? (by what you have commented):
# lines[lines == NA] = 1
# lines[lines == 2] = 3
# lines[lines == 0] = 2

#>Nevertheless, when you are calculating n02, n20, n22 you check for 0s and 2s so I guess you should do:
# lines[lines == NA] = 1
# or
# lines[is.na(lines)] <- 1
#>The graph shows nicely for me like that

genes <- nrow(smallRegion)
maxDist = 1000 # the maximum distance between genes we look at
index = 1
# allocate the maximum possibe space required upfront
lds = matrix(data = NA, nrow = genes * min(genes-1, maxDist * 2) / 2, ncol = 2)
colnames(lds) <- c('distance', 'ld')
for (i in seq(1,nrow(smallRegion)-1)) {
  for (j in seq(i+1, nrow(smallRegion))) {
    dist = abs(smallRegion[i, 'pos'] - smallRegion[j, 'pos']) # diff between pos of snp1 and snp2
    if (dist > maxDist) { # only look at pair of snips which are closer than maxDist apart
      break
    }
    snp1 = lines[i, ]
    snp2 = lines[j, ]
    na = sum(is.na(snp1 == snp2)) # find the number of indexes which has indexes in at least one of the columns
    N = 205 - na
    
    # calculate the number of lines which has either combination of genes for the 2 snips
    # n00 = length(intersect(which(snp1 == 0), which(snp2 == 0)))
    n02 = length(intersect(which(snp1 == 0), which(snp2 == 2)))
    n20 = length(intersect(which(snp1 == 2), which(snp2 == 0)))
    n22 = length(intersect(which(snp1 == 2), which(snp2 == 2)))
    
    # if (n00 + n02 + n20 + n22 != N) { print('something is wrong') }
      
    f22 = n22 / N
    f1 = (n22 + n20) / N
    f2 = (n22 + n02) / N
    
    ld = f22 - f1 * f2
    
    lds[index, 'distance'] = dist
    lds[index, 'ld'] = ld
    index = index + 1
  }
}


# only save those similarities which were actually computed
firstNa = min(which(is.na(lds[,'distance'])))
lds = lds[1:(firstNa-1), ]

lds = lds[order(lds[,'distance']), ] # just for fun

# Calculate the mean for each distance
lds_average = matrix(data = NA, nrow = maxDist, ncol = 2)
colnames(lds_average) <- c('distance', 'ld')
for (i in seq(1,1000)) {
  lds_average[i,'distance'] = i
  lds_average[i,'ld'] = mean(lds[lds[,'distance'] == i, 'ld'])
}

ggplot(as.data.frame(lds_average), aes(distance, ld)) +
  geom_smooth() + 
  geom_point() +
  theme_minimal() +
  labs(x = 'Distance between pairs', title = 'LD vs. distance between SNP\'s')

#>This last plot doesn't show like its suppose to be (it shows a straight line y = 0 for me)
#>I guess there is a problem in lines 158-160 so when its calculating n02, n20, n22 it doesn't find 0s and 2s

# Looks quite similar to the graph of the paper. There is some association between close SNPs compared to distant but exactly how they compare is difficult to say since they have calculated R^2.
