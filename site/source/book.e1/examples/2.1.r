# Reading in FAMuSS data
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 2.1 (Chi-squared test for association):
attach(fms)
NamesEsr1Snps <- names(fms)[substr(names(fms),1,4)=="esr1"]
fmsEsr1 <- fms[,is.element(names(fms),NamesEsr1Snps)]
Trait <- as.numeric(pre.BMI>25)
newFunction <- function(Geno){
	ObsTab <- table(Trait,Geno)
	return(chisq.test(ObsTab)$p.value)
	}
apply(fmsEsr1,2,newFunction)

Geno <- fmsEsr1[,2]
ObsTab <- table(Trait,Geno)
ObsTab

