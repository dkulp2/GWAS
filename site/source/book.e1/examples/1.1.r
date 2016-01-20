# Reading in FAMuSS data
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 1.1 (Identifying the minor allele and its frequency):
attach(fms)
GenoCount <- table(actn3_rs540874)
GenoCount
NumbObs <- sum(!is.na(actn3_rs540874))
GenoFreq <- as.vector(GenoCount/NumbObs)
GenoFreq
FreqA <- (2*GenoFreq[1] + GenoFreq[2])/2
FreqA
FreqG <- (GenoFreq[2] + 2*GenoFreq[3])/2
FreqG

# Alternatively, using genetics package:
install.packages("genetics")
library(genetics) 
Geno <- genotype(actn3_rs540874,sep="")
summary(Geno)
