#Code Snippet 2: SNP level filtering

# Setting thresholds
call<-0.9
minor<-0.01
hardy<-10^-6

# Filtering on MAF and call rate

# Find SNP summary statistics
snpsum.col<-col.summary(genotype)

use <-  (!is.na(snpsum.col$MAF)&snpsum.col$MAF > minor) & snpsum.col$Call.rate >= call
use[is.na(use)] <- FALSE
nrow(genoBim)-sum(use) #164080 SNPs will be removed

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col<-snpsum.col[use,]

# Filtering on HWE in CAD controls only
snpsum.colCont<-col.summary(genotype[row.names(genotype)%in%clinical[clinical$CAD==0,"FamID"],])
HWEuse<-(!is.na(snpsum.colCont$z.HWE)) & (abs(snpsum.colCont$z.HWE) < abs(qnorm((hardy/2))))
HWEuse[is.na(HWEuse)] <- FALSE
ncol(genotype)-sum(HWEuse) # 2208 SNPs removed

# Subset genotype and SNP summary data for SNPs that pass HWE criteria
genotype <- genotype[,HWEuse]
snpsum.col<-snpsum.col[HWEuse,]
dim(genotype)
#695185 SNPs remain

rm(snpsum.colCont)

# Write subsetted genotype file
save(genotype, file="subsetted_genotype.RData")
load("subsetted_genotype.RData")
