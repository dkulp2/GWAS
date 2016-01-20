# Reading in FAMuSS data
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Necessary code from Example 2.1:
attach(fms)
NamesEsr1Snps <- names(fms)[substr(names(fms),1,4)=="esr1"]
fmsEsr1 <- fms[,is.element(names(fms),NamesEsr1Snps)]
Trait <- as.numeric(pre.BMI>25)

# Example 2.2 (Fishers exact test for association):
newFunction <- function(Geno){
	ObsTab <- table(Trait,Geno)
	return(fisher.test(ObsTab)$p.value)
	}
	
apply(fmsEsr1,2,newFunction)
