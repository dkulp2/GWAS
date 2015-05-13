# ---- code5 ----
# Genotype imputation

source("globals.R")

# load data created in previous snippets
load(working.data.fname)     # loads genotype, genoBim and clinical

# ---- code5-a ----

library(chopsticks)                     # HapMap routines

# read in hapmap data for chromosomes 16
HapMap16 <- read.HapMap.data(sprintf(hapmap.url, 16))

detach("package:chopsticks", unload=TRUE) # remove to unmask snpStats functions
library(snpStats)

# Impute SNPs

# Subset for SNPs on given chromosome
presSnps <- colnames(genotype)
presDatChr <- genoBim[genoBim$SNP %in% presSnps & genoBim$chr==chrNum, ]
targetSnps <- presDatChr$SNP

# Obtain hapmap data for given chromosome
hapMatrix <- new("SnpMatrix",HapMap16$snp.data)
support <- HapMap16$snp.support

# Subset hapmap data for our SNPs
# "missing" and "present" are snpMatrix objects needed for imputation rules
is.present <- colnames(hapMatrix) %in% targetSnps
missing <- hapMatrix[,!is.present]
present <- hapMatrix[,is.present]

# Obtain positions of SNPs to be used for imputation rules
pos.pres <- support$Position[is.present]
pos.miss <- support$Position[!is.present]

# Calculate and store imputation rules using snp.imputation()
rules <- snp.imputation(present,missing,pos.pres,pos.miss)

# Obtain 'best guess' genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules, target, as.numeric=FALSE)
print(dim(imputed))  # 21342 SNPs were imputed

# ---- code5-end ----

rm(hapMatrix)
rm(missing)
rm(present)
rm(HapMap16)

# Add new imputed, target and rules data to saved results
save(genotype, genoBim, clinical, pcs, imputed, target, rules, file=working.data.fname)
