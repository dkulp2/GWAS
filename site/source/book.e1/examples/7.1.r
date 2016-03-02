# Reading in Virco data:
virco <- read.csv("http://stat-gen.org/book.e1/data/Virco_data.csv", header=T, sep=",")

# Example 7.1 (An application of random forests):
install.packages("randomForest")
library(randomForest)
attach(virco)
Trait <- NFV.Fold - IDV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
RegRF <- randomForest(VircoGeno.c, Trait.c, importance=TRUE)
RegRF
varImpPlot(RegRF,main="")
