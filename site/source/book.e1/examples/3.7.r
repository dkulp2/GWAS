# Reading in HGDP data:
hgdp <- read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt", header=T, sep="\t")

# Example 3.7 (Testing for HWE using Fishers exact test):
attach(hgdp)
Akt1Snp1Maya <- AKT1.C0756A[Population=="Maya"]
ObsCount <- table(Akt1Snp1Maya)
ObsCount
Nobs <- sum(ObsCount)
FreqC <- (2 * ObsCount[3] + ObsCount[2])/(2*Nobs)
ExpCount <- c(Nobs*(1-FreqC)^2, 2*Nobs*FreqC*(1-FreqC),Nobs*FreqC^2)
ExpCount
n11 <- ObsCount[3]
n12 <- ObsCount[2]
n22 <- ObsCount[1]
n1 <- 2*n11+n12
Num <- 2^n12 * factorial(Nobs)/prod(factorial(ObsCount))
Denom <- factorial(2*Nobs) / (factorial(n1)*factorial(2*Nobs-n1))
FisherP1 <- Num/Denom
FisherP1

library(genetics)
Akt1Snp1Maya <- genotype(AKT1.C0756A[Population=="Maya"], sep="")
HWE.exact(Akt1Snp1Maya)