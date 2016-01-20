# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 3.3 (Measuring LD based on r^2 and the \chi^2-statistic):
attach(fms)
library(genetics)
Actn3Snp1 <- genotype(actn3_r577x,sep="")
Actn3Snp2 <- genotype(actn3_rs540874,sep="")
LD(Actn3Snp1,Actn3Snp2)$"R^2"
LD(Actn3Snp1,Actn3Snp2)