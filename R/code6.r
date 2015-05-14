# ---- code6 ----
# Association analysis of typed SNPs (using parallel processing)

source("globals.R")

# load derived data from previous snippets
load(working.data.fname)

library(snpStats)
library(plyr)

# ---- code6-a ----
library(GenABEL)
library(parallel)
source("GWAA.R")

# Create data frame from genotype file and merge with clincal data to create
# phenotype table
pcs.df <- data.frame(FamID=row.names(genotype),pcs)
phenoSub <- merge(clinical,pcs.df)      # data.frame => [ FamID CAD sex age hdl pc1 pc2 ... pc10 ]

# We will do a rank-based inverse normal transformation of hdl
phenoSub$phenotype <- rntransform(phenoSub$hdl, family="gaussian")

# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(phenoSub$hdl, main="Histogram of HDL", xlab="HDL")
hist(phenoSub$phenotype, main="Histogram of Tranformed HDL", xlab="Transformed HDL")

# Remove hdl column from table
phenoSub$hdl <- NULL
phenoSub$CAD <- NULL

# Rename columns to match names necessary for GWAS() function
phenoSub <- rename(phenoSub, replace=c(FamID="id"))

# Include only subjects with hdl data
phenoSub<-phenoSub[!is.na(phenoSub$phenotype),]
# 1309 subjects included with phenotype data

print(head(phenoSub))

# ---- code6-b ----

# Run GWAS analysis
# Note: This function writes a file, but does not produce an R object
start <- Sys.time()
GWAA(genodata=genotype, phenodata=phenoSub, filename=gwaa.fname)
end <- Sys.time()
print(end-start)

# ---- code6-end ----

# Add phenosub to saved results
save(genotype, genoBim, clinical, pcs, imputed, target, rules, phenoSub, file=working.data.fname)
