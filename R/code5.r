#Code Snippet 5: Genotype imputation

source("globals.R")

# load data created in previous snippets
load(genotype.subset.fname)     # loads genotype, genoBim and clinical

##################

# read in hapmap data for chromosomes 16
chrNum<-16

library(chopsticks)                     # HapMap routines

# Imputation of non-typed HapMap SNPs
presSnps <- colnames(genotype)

HapMap16 <- read.HapMap.data(sprintf(hapmap.url, chrNum))

detach("package:chopsticks", unload=TRUE) # remove to unmask snpStats functions
library(snpStats)

# Impute SNPs on given chromosome

# Subset for SNPs on given chromosome
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

rm(hapMatrix)
rm(missing)
rm(present)
rm(HapMap16)

##################

# Add new imputed, target and rules data to saved results
save(genotype, genoBim, clinical, pcs, imputed, target, rules, file=genotype.subset.fname)
