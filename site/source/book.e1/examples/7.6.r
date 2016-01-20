# Reading in Virco data:
virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", header=T, sep=",")

# Example 7.6 (Monte Carolo logic regression):
library(LogicReg)
attach(virco)
Trait <- SQV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
VircoLogicRegMCMC <- logreg(resp=Trait.c, bin=VircoGeno.c, select=7)
plot(sort(VircoLogicRegMCMC$single), xlab="Sorted SNPs", ylab="Number of selected models")
names(VircoGeno)[order(VircoLogicRegMCMC$single)]

