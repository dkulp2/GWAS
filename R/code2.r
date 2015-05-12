#Code Snippet 2: SNP level filtering

source("globals.R")

# load data created in previous snippet
load(genotype.subset.fname)

##################

# Setting thresholds
call <- 0.95
minor <- 0.01

library(snpStats)

# Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)

# Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use),"SNPs will be removed due to low MAF or call rate.\n") #203287 SNPs will be removed

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]

print(dim(genotype))                           # 658186 SNPs remain

##################

# Write subsetted genotype data and derived results for future use
save(genotype, snpsum.col, genoBim, clinical, file=genotype.subset.fname)
