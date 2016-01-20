# Reading in Virco data:
virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", header=T, sep=",")

# Example 6.2 (Creating a classification tree):
attach(virco)
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait <- as.factor(IDV.Fold > NFV.Fold)
library(rpart)
ClassTree <- rpart(Trait~., method="class", data=VircoGeno)
ClassTree
plot(ClassTree)
text(ClassTree)
rpart(Trait~., method="class", parms=list(split='information'), data=VircoGeno)
rpart(Trait~., method="class", parms=list(split='gini'), control=rpart.control(minsplit=150, minbucket=50), data=VircoGeno)

