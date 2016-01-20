# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 3.1 (Measuring LD using D-prime):
library(genetics)
attach(fms)
actn3_r577x[1:10]
actn3_rs540874[1:10]
Actn3Snp1 <- genotype(actn3_r577x,sep="")
Actn3Snp2 <- genotype(actn3_rs540874,sep="")
Actn3Snp1[1:10]
LD(Actn3Snp1,Actn3Snp2)$"D'"
Esr1Snp1 <- genotype(esr1_rs1801132,sep="")
LD(Actn3Snp1,Esr1Snp1)$"D'"
