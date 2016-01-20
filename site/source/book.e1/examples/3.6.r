# Reading in HGDP data:
hgdp <- read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt", header=T, sep="\t")

# Example 3.6 (Testing for HWE using Pearsons \chi^2-test):
attach(hgdp)
Akt1Snp1 <- AKT1.C0756A
ObsCount <- table(Akt1Snp1)
Nobs <- sum(ObsCount)
ObsCount
FreqC <- (2 * ObsCount[3] + ObsCount[2])/(2*Nobs)
ExpCount <- c(Nobs*(1-FreqC)^2, 2*Nobs*FreqC*(1-FreqC),Nobs*FreqC^2)
ExpCount
ChiSqStat <- sum((ObsCount - ExpCount)^2/ExpCount)
ChiSqStat
qchisq(1-0.05,df=1)

library(genetics)
Akt1Snp1 <- genotype(AKT1.C0756A, sep="")
HWE.chisq(Akt1Snp1)