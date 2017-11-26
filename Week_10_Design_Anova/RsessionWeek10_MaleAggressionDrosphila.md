## Applying comparisons of means and the ANOVA framework to detect QTL in the Drosophila lines
# Source of the data
Shorter et al 2015 Genetic architecture of natural variation in Drosophila melanogaster aggressive behavior, PNAS
https://www.pnas.org/cgi/doi/10.1073/pnas.1510104112

# Background
Measurements of male aggression are available for 200 lines out of the 205 lines of the DRGP pannel.
So in principle we can combinethe SNP data with this scores to examine if some SNPs in the genome can explain / predict aggression scores
Results of the SNP giving the top associations are given as supplementary table of the PNAS publication
(Dataset S2. Genome wide association analyses for aggressive behavior.), available as Excel file

## Getting the aggression score data

Pheno=read.csv("aggression.male.csv", header =F)
names(Pheno)=c("LineId","Aggression")

If you want to "know"watch more about male aggression in Drosphila:
https://www.youtube.com/watch?v=uJaDoTigvEI


1. Vizualize and summarize the distribution of male agression scores among the lines. Read in the supplementary info how this was scored.

2. Seclect the 4 SNPs reported in the publication as located on chr 3L and showing the most significant association between aggression and SNP genoytpe
For each SNP:
Measure and make a grpah illustratting the mean difference in aggression score for each group of flies.
A goup is here defined by the SNP genotype of each line (i.e. you are 0 or 2, note that some lines could not be genotyped at certain SNPs).
test the null hypothesis that agression scores have identical means between both groups
Discuss / Justify which test you used and how you obtained a p-value for the test.
Use an ANOVA with a single factor (fixed factor genotype) and calculate R^2 of the ANOVA model.
Tip in R the lm() function allows you to do ANOVAs.
