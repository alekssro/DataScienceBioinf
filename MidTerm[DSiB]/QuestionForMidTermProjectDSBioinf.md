# Questions on the DRGP dataset

Reference paper describes the data obtained here.

Paper source: The Drosophila melanogaster Genetic Reference Panel, Nature 2012

Data source

 * http://dgrp2.gnets.ncsu.edu/  (the project)
 * http://dgrp2.gnets.ncsu.edu/data.html  (the data source tabular 3)

A zipped version of the data is available here: https://www.dropbox.com/s/257g23neey8gsny/dgrp2.tgeno.zip?dl=0


## Q0. Import the data (warning the file is big, consider using skip and nmax when reading in data)
Describe briefly the structure of the data and what the different variables mean
- how many different chromosomes do you have?

## Q1. Extract the first 100,000 SNPs located on chromosome 3L,
- Locate the genomic position of the first and last SNP
- filter out all SNPs for which less than 75%, 80%, 90%, 95%  of the individuals could be genotyped, report how many SNPs are left based on the filter you apply

In all subsequent questions, restrict your attention to SNP genotyped in at least 95% of individuals

## Q2. Graph the distribution of SNP allele frequencies, graph the distribution of coverage of SNPs ("cov") and the distribution of number of lines genotyped for each snp

 - Is there a statistical association the different types of SNPs and the allele frequency of SNPs?

The type of snp is defined in the last part of the "id" column.

We define a SNP as either rare (maf < 0.05) or common (maf >= 0.05).

We define coverage as either "low coverage" (coverage in the lowest 5% quantile) or "high coverage".

## Q3. If coverage was homogenous throughout the genome (by that we mean that on average the coverage is the same for any given posiiton),
 - what probability distribution is expected to capture well the coverage ?

## Q4. Make goodness of fit test of the coverage data: is the theoretical distribution proposed in Q3 a good fit for the data?

## Q5. Pick yourself a small (100 kb) region on chromosome 3L
 - Graph the frequency of "0" alleles along that region
 - Calculate the statistical association between allelic states of each pair of neighboring SNPS (this association is also called LD)
 - Explore visually how LD varies with the physical distance beween SNPs (measured in bp)
 - Compare and discuss with the Fig 1 of the paper.
  
  