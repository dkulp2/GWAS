# Reading in Virco data:
virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", header=T, sep=",")

# Example 7.7 (An application of MARS):
attach(virco)
Trait <- NFV.Fold - IDV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
install.packages("earth")
library(earth)
VircoMARS <- earth(Trait.c~., data=VircoGeno.c, degree=2)
summary(VircoMARS)
evimp(VircoMARS)
