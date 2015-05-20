#Code Snippet 7: Association analysis of imputed SNPs

source("globals.R")

# load derived data from previous snippets
load(working.data.fname)

##################

library(snpStats)
library(plyr)

# Carry out association testing for imputed SNPs using snp.rhs.tests()
rownames(phenoSub) <- phenoSub$id

imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, family = "Gaussian", data = phenoSub, snp.data = target, rules = rules)

# Obtain p values for imputed SNPs
pVals <- p.value(imp)
results <- data.frame(SNP = imp@snp.names, p.value = pVals)
results <- results[!is.na(results$p.value),]


#Write a file containing the results
write.csv(results, impute.out.fname, row.names=FALSE)

# Merge imputation testing results with support to obtain coordinates
imputeOut<-merge(results, support[, c("SNP", "position")])
imputeOut$chr <- 16

imputeOut$type <- "imputed"

# Find the -log_10 of the p-values
imputeOut$Neg_logP <- -log10(imputeOut$p.value)

# Order by p-value
imputeOut <- arrange(imputeOut, p.value)
print(head(imputeOut))

# Read in file containing protein coding genes coords
genes <- read.csv(protein.coding.coords.fname)

# Subset for CETP SNPs
source("map2gene.R")
impCETP <- map2gene("CETP", coords = genes, SNPs = imputeOut)

impCETPgeno <- imputed[, colnames(imputed) %in% impCETP$SNP ]

##################

save(genotype, genoBim, clinical, pcs, imputed, target, rules,
     phenoSub, genes, impCETP, impCETPgeno, imputeOut, file = working.data.fname)

