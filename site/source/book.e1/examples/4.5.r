# Reading in Virco data:
virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", header=T, sep=",")

# Necessary code from Example 4.1:
attach(virco)
PrMut <- virco[,23:121]!="-" & virco[,23:121]!="."
NObs <- dim(virco)[1]
PrMutSub <-data.frame(PrMut[ , apply(PrMut,2,sum) > NObs*.05])
Trait <- IDV.Fold - NFV.Fold
TtestP <- function(Geno){
	return(t.test(Trait[Geno==1],
		Trait[Geno==0], na.rm=T)$"p.value")
	}
Pvec <- apply(PrMutSub, 2, TtestP)
Pvec <- as.vector(Pvec)

# Example 4.5 (Calculation of the q-value)
library(qvalue)
sort(qvalue(Pvec, lambda=0)$qvalues)
sort(qvalue(Pvec, pi0.method="bootstrap")$qvalues)
qvalue(Pvec, pi0.method="bootstrap")$pi0
