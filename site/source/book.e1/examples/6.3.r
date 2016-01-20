# Reading in Virco data:
virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", header=T, sep=",")

# Necessary code from Example 6.2:
attach(virco)
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")

# Example 6.3 (Generating a regression tree):
library(rpart)
Trait <- NFV.Fold - IDV.Fold
Tree <- rpart(Trait~., method="anova", data=VircoGeno)
Tree
