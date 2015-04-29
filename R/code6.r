#Code Snippet 6: Association analysis of typed SNPs (using parallel processing)

load("GWAStutorial_pcs.RData")

# Create data frame from genotype file and merge with clincal data to create
# phenotype table
pcs<-data.frame(FamID=row.names(genotype),pcs)
phenoSub<-merge(clinical,pcs,by.x="FamID",by.y="FamID")

# We will do a rank-based inverse normal transformation of hdl
require(GenABEL)
phenoSub$hdl.rnt<-rntransform(phenoSub$hdl, family="gaussian")

# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(phenoSub$hdl, main="Histogram of HDL", xlab="HDL")
hist(phenoSub$hdl.rnt, main="Histogram of Tranformed HDL", xlab="Transformed HDL")

# Remove hdl column from table
phenoSub$hdl<-NULL
phenoSub$CAD<-NULL

# Rename columns to match names necessary for GWAS() function
colnames(phenoSub)[colnames(phenoSub)=="hdl.rnt"]<-"phenotype"
colnames(phenoSub)[colnames(phenoSub)=="FamID"]<-"id"

# Include only subjects with hdl data
phenoSub<-phenoSub[!is.na(phenoSub$phenotype),]
# 1309 subjects included with phenotype data

library(parallel)

# Read in GWAA() function
source("GWAAfunction.R")

# Run GWAS analysis
# Note: This function writes a file, but does not produce an R object
start<-Sys.time()
GWAA(genodata=genotype, phenodata=phenoSub, filename="GWAStutorialout.txt", nSplits=10)
end<-Sys.time()
end-start
