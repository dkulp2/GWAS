#Code Snippet 8: Data Integration

source("globals.R")

# load data derived in previous snippets
load(working.data.fname)

##################

library(plyr)

# (1) Read in GWAS output, (2) add position and chromosome number, (3) order by significance
GWASout <- read.table(gwaa.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Find the -log_10 of the p-values
GWASout$Neg_logP <- -log10(GWASout$p.value)

# Merge output with genoBim to obtain coordinates
GWASout <- merge(GWASout, genoBim[,c("SNP", "chr", "position")])
rm(genoBim)

# Order SNPs by signal
GWASout <- arrange(GWASout, -Neg_logP)
print(head(GWASout))

# Combine typed and imputed
GWASout$type<-"typed"
GWAScomb<-rbind.fill(GWASout, imputeOut)

# Subset for CETP SNPs
source("map2gene.R")
typCETP <- map2gene("CETP", coords = genes, SNPs = GWASout)

# Combine CETP SNPs from imputed and typed analysis
CETP <- rbind.fill(typCETP, impCETP)
write.csv(CETP, CETP.fname, row.names=FALSE)

##################

save(genotype, clinical, pcs, imputed, target, rules, phenoSub, genes,
     impCETP, impCETPgeno, GWASout, GWAScomb, CETP, file=genotype.subset.fname)
