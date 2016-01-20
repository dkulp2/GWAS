# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Necessary code from Example 3.1:
attach(fms)
library(genetics)
Actn3Snp1 <- genotype(actn3_r577x,sep="")
Actn3Snp2 <- genotype(actn3_rs540874,sep="")

# Example 3.2 (Measuring pairwise LD for a group of SNPs):
Actn3Snp3 <- genotype(actn3_rs1815739,sep="")
Actn3Snp4 <- genotype(actn3_1671064,sep="")
Actn3AllSnps <- data.frame(Actn3Snp1,Actn3Snp2,Actn3Snp3,Actn3Snp4)
LD(Actn3AllSnps)$"D'"
install.packages("LDheatmap")
library(LDheatmap)
LDheatmap(Actn3AllSnps, LDmeasure="D'")