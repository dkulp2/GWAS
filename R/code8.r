# ---- code8 ----
# Data Integration

source("globals.R")

# load data derived in previous snippets
load(working.data.fname(7))

library(plyr)
source("map2gene.R")

# ---- code8-a ----
# Read in GWAS output that was produced by GWAA function
GWASout <- read.table(gwaa.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Find the -log_10 of the p-values
GWASout$Neg_logP <- -log10(GWASout$p.value)

# Merge output with genoBim by SNP name to add position and chromosome number
GWASout <- merge(GWASout, genoBim[,c("SNP", "chr", "position")])
rm(genoBim)

# Order SNPs by significance
GWASout <- arrange(GWASout, -Neg_logP)
print(head(GWASout))

# ---- code8-b ----
# Combine typed and imputed
GWASout$type <- "typed"

GWAScomb<-rbind.fill(GWASout, imputeOut)
head(GWAScomb)
tail(GWAScomb)

# Subset for CETP SNPs
typCETP <- map2gene("CETP", coords = genes, SNPs = GWASout)

# Combine CETP SNPs from imputed and typed analysis
CETP <- rbind.fill(typCETP, impCETP)[,c("SNP","p.value","Neg_logP","chr","position","type","gene")]
print(CETP)
# ---- code8-end ----


write.csv(CETP, CETP.fname, row.names=FALSE) # save for future use

save(genotype, clinical, pcs, imputed, target, rules, phenoSub, support, genes,
     impCETP, impCETPgeno, imputeOut, GWASout, GWAScomb, CETP, file=working.data.fname(8))
