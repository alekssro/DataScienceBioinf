Meta data from Medicago truncatula dataset
Source (after trimming a few extra columns): supp data SD02 of Branca et al
[DataSource](http://www.pnas.org/content/suppl/2011/09/22/1104032108.DCSupplemental)

100 kb window genomic summary statistics
column headings:
X window number (irrelevant)
LdH_SNPS: Number of SNPs in a window used to calculate LD and recombination rate
LdH_start: start of window
LdH_win_size: physical windows size
rho_kb: estimated population scaled recombination rate per kb
rho_theta: rho divided by theta
bases: number of bases with data in more than 20 accessions
mutations: number of mutations
SNPs: Single Nucleotide Polymorphisms singletons
TajimasD: statistic describing how many rare SNP are found relative to common ones. Zero is the expectation under selective neutrality in a stable population. An excess of rare alleles yields negative Tajima's D value and conversely an excess of frequent SNP yields positive values.

qp.site (also called Theta_pi) scaled mutation as estimated through mean nucleotide diversity per site

qw.site (also called Theta_w) scaled mutation as estimated through the mean number of segregating sites (standardized per site)

dist_cent: distance from the start of the window to the centromere (negative values are used for the left arm)
prop_dist: distance to the centromere scaled by the length of the chromosomal arm
gene_dens: gene density estimated as the proportion of window that is predicted to be coding.
