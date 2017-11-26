# Aspirin data from book ..  
library(abd)
data(Aspirin)
plot(Aspirin)

# Exercise 33 ch 10 
# https://www.ncbi.nlm.nih.gov/pubmed/19923577
kuru=matrix(data=c(13,77,14,22,12,14), nrow =2, ncol =3, byrow = T) ## Genotypes at codon 129 of the prion protein gene (PRNP) 
rownames(kuru)=c("Elderly", "Young")
colnames(kuru) = c("MM", "MV", "VV")
mosaicplot(t(kuru), main="", col=c("red","yellow"), cex=2)
## The Chi^2 goodness of fit test 
chisq.test(kuru)
myX2Test=chisq.test(kuru)
names(myX2Test)
myX2Test
myX2Test$observed
myX2Test$expected
myX2Test$residuals
myX2Test$method

fisher.test(kuru)
