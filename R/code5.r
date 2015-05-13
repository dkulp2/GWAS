#Code Snippet 5: Genotype imputation

source("globals.R")

# load data created in previous snippets
load(genotype.subset.fname)     # loads genotype, genoBim and clinical

##################

library(snpStats)

# read in 1000g data for chromosomes 16
chrNum <- 16

# Read in 1000g data for given chromosome 16
thougeno <- read.pedfile(onethou.fn$ped, snps = onethou.fn$info, which=1)

# Obtain genotype data for given chromosome
genoMatrix <- thougeno$genotypes
support <- thougeno$map
colnames(support)<-c("SNP", "Position", "A1", "A2")

# Imputation of non-typed 1000g SNPs
presSnps <- colnames(genotype)

# Subset for SNPs on given chromosome
presDatChr <- genoBim[genoBim$SNP %in% presSnps & genoBim$chr==chrNum, ]
targetSnps <- presDatChr$SNP

# Subset 1000g data for our SNPs
# "missing" and "present" are snpMatrix objects needed for imputation rules
is.present <- colnames(genoMatrix) %in% targetSnps
missing <- genoMatrix[,!is.present]
present <- genoMatrix[,is.present]

# Obtain positions of SNPs to be used for imputation rules
pos.pres <- support$Position[is.present]
pos.miss <- support$Position[!is.present]

# Calculate and store imputation rules using snp.imputation()
rules <- snp.imputation(present, missing, pos.pres, pos.miss)

# Remove failed imputations
rules <- rules[can.impute(rules)]
cat("Imputation rules for",length(rules),"SNPs were estimated\n")  # Imputation rules for 197888 SNPs were estimated

# Quality control for imputation certainty and 
# Set thresholds
r2threshold <- 0.7
minor <- 0.01

# Filter on imputation certainty and MAF
rules <- rules[imputation.r2(rules)>=r2threshold]

cat(length(rules),"imputation rules remain after uncertain impuations were removed\n")  # 162565 imputation rules remain after uncertain impuations were removed

rules <- rules[imputation.maf(rules)>=minor]
cat(length(rules),"imputation rules remain after MAF filtering\n")  # 162565 imputation rules remain after MAF filtering


# Obtain 'best guess' genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules, target, as.numeric=FALSE)
print(dim(imputed))  # 162565 SNPs were imputed

rm(genoMatrix)
rm(missing)
rm(present)

##################

# Add new imputed, target and rules data to saved results
save(genotype, genoBim, clinical, pcs, imputed, target, rules, file=genotype.subset.fname)
