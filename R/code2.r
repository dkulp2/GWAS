# ---- code2 ----
# SNP level filtering

source("globals.R")

# load data created in previous snippet
load(working.data.fname(1))

library(snpStats)

# ---- code2-a ----
# Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)
print(head(snpsum.col))

# ---- code2-b ----
# Setting thresholds
call <- 0.95
minor <- 0.01

# Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use),"SNPs will be removed due to low MAF or call rate.\n") #203287 SNPs will be removed

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]

print(genotype)                           # 658186 SNPs remain

# ---- code2-end ----

# Write subsetted genotype data and derived results for future use
save(genotype, snpsum.col, genoBim, clinical, file=working.data.fname(2))
