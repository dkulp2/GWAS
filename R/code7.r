#Code Snippet 7: Association analysis of imputed SNPs

source("globals.R")

# load derived data from previous snippets
load(genotype.subset.fname)

##################

library(snpStats)
library(plyr)

# Carry out association testing for imputed SNPs using single.snp.tests()
imp <- single.snp.tests(phenotype=hdl, stratum=sex, data=clinical, snp.data=target, rules=rules)

# Obtain p values for imputed SNPs
pVals <- p.value(imp,df=1)
pVals <- pVals[!is.na(pVals)]
results <- data.frame(SNP=names(pVals), p.value=pVals)

#Write a file containing the results
write.csv(results, impute.out.fname, row.names=FALSE)

# Merge imputation testing results with support to obtain coordinates
# Read in file containing SNPs mapped to protein coding genes within 5kb
genes<-read.csv(protein.coding.snps.fname)
imputeOut<-merge(results, genes)

# Find the -log_10 of the p-values
imputeOut$Neg_logP <- -log10(imputeOut$p.value)

# Order by p-value
imputeOut <- arrange(imputeOut, p.value)
print(head(imputeOut))

# Subset for CETP SNPs
impCETP <- imputeOut[imputeOut$gene=="CETP",c("SNP", "chr", "position","p.value", "Neg_logP")]
impCETP$type <- "imputed"

impCETPgeno <- imputed[, colnames(imputed) %in% impCETP$SNP ]

##################

save(genotype, genoBim, clinical, pcs, imputed, target, rules, phenoSub, impCETP, impCETPgeno, file=genotype.subset.fname)
