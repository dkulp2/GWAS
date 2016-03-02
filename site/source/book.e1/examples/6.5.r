# Reading in Virco data:
virco <- read.csv("http://stat-gen.org/book.e1/data/Virco_data.csv", header=T, sep=",")

# Example 6.5 (Cost-complexity pruning):
attach(virco)
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"])
library(rpart)
Tree <- rpart(APV.Fold~.,data=VircoGeno)
Tree
printcp(Tree)
plotcp(Tree)
pruneTree <- prune(Tree,cp=.03)
pruneTree
