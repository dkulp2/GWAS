# Reading in FAMuSS data
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 2.3 (Cochran-Armitage (C-A) trend test for association):
install.packages("coin")
library(coin)
attach(fms)
Geno <- esr1_rs1042717
Trait <- as.numeric(pre.BMI>25)
GenoOrd <- ordered(Geno)
independence_test(Trait~GenoOrd,teststat="quad",
 	scores=list(GenoOrd=c(0,1,2)))
