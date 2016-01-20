# Reading in Virco data:
virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", header=T, sep=",")

# Example 7.5 (Application of logic regression):
attach(virco)
Trait <- NFV.Fold - IDV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
install.packages("LogicReg")
library(LogicReg)
VircoLogicReg <- logreg(resp=Trait.c, bin=VircoGeno.c, select=1)
plot(VircoLogicReg)
VircoLogicReg
VircoLogicRegMult <- logreg(resp=Trait.c, bin=VircoGeno.c, select=2, ntrees=2, nleaves=8)
plot(VircoLogicRegMult)
VircoLogicRegMult
