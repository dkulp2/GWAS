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

# Merge SNP output with file containing protein coding SNP information
GWASgenes <- merge(GWASout, genes[,c("SNP", "gene")], sort=FALSE) # Merge gene file with GWASout
CETP <- GWASgenes[GWASgenes$gene=="CETP",] # Subset CETP SNPs

# Combine CETP SNPs from imputed and typed analysis
CETP <- CETP[,c("SNP", "chr", "position","p.value", "Neg_logP")]
CETP$type <- "typed"
CETP <- rbind(CETP, impCETP)
CETP <- arrange(CETP, p.value)
write.csv(CETP, CETP.fname, row.names=FALSE)
