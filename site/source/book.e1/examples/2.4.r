# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 2.4 (Two-sample tests of association for a quantitative trait)
attach(fms)
NamesResistinSnps <- names(fms)[substr(names(fms),1,8)=="resistin"]
fmsResistin <- fms[,is.element(names(fms),NamesResistinSnps)]
library(genetics)
TtestPval <- function(Geno){
	alleleMajor <- allele.names(genotype(Geno, sep="", reorder="freq"))[1]
	GenoWt <- paste(alleleMajor, alleleMajor, sep="")
	GenoBin <- as.numeric(Geno!=GenoWt)[!is.na(Geno)] 
	Trait <- NDRM.CH[!is.na(Geno)]
	return(t.test(Trait[GenoBin==1],Trait[GenoBin==0])$p.value)
	}
apply(fmsResistin,2,TtestPval)
Geno <- fms$"resistin_c180g"
table(Geno)
GenoWt <- names(table(Geno))[table(Geno)==max(table(Geno))]
GenoWt
GenoBin <- as.numeric(Geno!=GenoWt)[!is.na(Geno)] 
Trait <- NDRM.CH[!is.na(Geno)]
t.test(Trait[GenoBin==1],Trait[GenoBin==0])

