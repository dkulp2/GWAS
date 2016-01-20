# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 7.4 (MIRF):
install.packages("mirf")
library(mirf)
attach(fms)
genoSNPnames <- c("actn3_r577x","actn3_rs540874","actn3_rs1815739","actn3_1671064","resistin_c30t","resistin_c398t",
	"resistin_g540a","resistin_c980g","resistin_c180g","resistin_a537c")
FMSgeno <- fms[,is.element(names(fms),genoSNPnames)][!is.na(NDRM.CH),]
Geno <- sepGeno(FMSgeno)
Trait <- NDRM.CH[!is.na(NDRM.CH)]
mirf(geno=Geno$geno, y=Trait, gene.column=c(8,12), SNPnames=genoSNPnames, M=10)
