#Code Snippet 9: Visualizing and QC of GWA findings

source("globals.R")

# load derived data from previous snippets
load(genotype.subset.fname)

##################

library(GenABEL)
library(LDheatmap)
library(rtracklayer)

# Create Manhattan Plot
source("GWAS_ManhattanFunction.R")
par(mfrow=c(1,1))
GWAS_Manhattan(GWASout)

# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[,c("id","phenotype")]
GWAA(genodata=genotype, phenodata=phenoSub2, filename=gwaa.unadj.fname, nSplits=10)
GWASoutUnadj <- read.table(gwaa.unadj.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Create QQ plots for adjusted and unadjusted model outputs
estlambda(GWASout$t.value^2,plot=TRUE,method="median")
estlambda(GWASoutUnadj$t.value^2,plot=TRUE,method="median")

# Combine genotypes and imputed genotypes for CETP region
subgen <- genotype[, colnames(genotype) %in% CETP$SNP ]   # CETP subset of snpMatrix from typed SNPs
subgen <- cbind(subgen, impCETPgeno)                      # CETP subset of snpMatrix from imputed SNPs 

# Order CETP SNPs by position
CETP <- arrange(CETP, position)
subgen <- subgen[, order(match(colnames(subgen),CETP$SNP)) ]

# Create LDheatmap
ld <- ld(subgen,subgen, stats="R.squared") # Find LD map of CETP SNPs

ll <- LDheatmap(ld,CETP$position,flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination, and scatterplot
llplusgenes <- LDheatmap.addGenes(ll, chr="chr16", genome="hg19")
llGenesRecomb <- LDheatmap.addRecombRate(llplusgenes, chr="chr16", genome="hg19", view = "full")
llGenesRecombScatter <- LDheatmap.addScatterplot(llGenesRecomb,CETP$Neg_logP, ylab="-log10 p-value")
