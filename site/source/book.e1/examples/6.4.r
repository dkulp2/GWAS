# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 6.4 (Categorical and ordinal predictors in a tree):
attach(fms)
Trait <- NDRM.CH
library(rpart)
RegTree <- rpart(Trait~resistin_c30t+resistin_c398t+
	resistin_g540a+resistin_c980g+resistin_c180g+
	resistin_a537c, method="anova")
RegTree
RegTreeOr <- rpart(Trait~as.numeric(resistin_c30t)+
	as.numeric(resistin_c398t)+as.numeric(resistin_g540a)+
	as.numeric(resistin_c980g)+as.numeric(resistin_c180g)+
	as.numeric(resistin_a537c), method="anova")
RegTreeOr
