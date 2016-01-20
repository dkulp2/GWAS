# Reading in HGDP data:
hgdp <- read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt", header=T, sep="\t")

# Example 3.8 (HWE and geographic origin):
attach(hgdp)
table(Geographic.area)
library(genetics)
Akt1Snp1 <- genotype(AKT1.C0756A, sep="")
HWEGeoArea <- tapply(Akt1Snp1,INDEX=Geographic.area,HWE.chisq)  
HWEGeoArea$"Central Africa"
HWEGeoArea$"South America"