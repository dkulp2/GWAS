#Code Snippet 2: SNP level filtering

source("globals.R")

# load data created in previous snippet
load(genotype.subset.fname)

##################

# Setting thresholds
call <- 0.9
minor <- 0.01
hardy <- 10^-6

library(snpStats)

# Create SNP summary statistics (MAF, Call rate, Hardy-Weinberg, etc.)
snpsum.col <- col.summary(genotype)

# Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use),"SNPs will be removed due to low MAF or call rate.\n") #164080 SNPs will be removed

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col<-snpsum.col[use,]

# Filtering on HWE in CAD controls only
CADcontrols <- clinical[ clinical$CAD==0, 'FamID' ]
snpsum.colCont <- col.summary( genotype[CADcontrols,] )
HWEuse <- with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)

HWEuse[is.na(HWEuse)] <- FALSE          # Remove NA's as well
cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n")  # 2208 SNPs removed

# Subset genotype and SNP summary data for SNPs that pass HWE criteria
genotype <- genotype[,HWEuse]
snpsum.col <- snpsum.col[HWEuse,]

print(dim(genotype))                           # 695185 SNPs remain

##################

# Write subsetted genotype data and derived results for future use
save(genotype, snpsum.col, genoBim, clinical, file=genotype.subset.fname)
