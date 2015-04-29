#Code Snippet 9: Visualizing and QC of GWA findings

# Create Manhattan Plot
source("GWAS_ManhattanFunction.R")
par(mfrow=c(1,1))
GWAS_Manhattan(GWASout)

# Create QQ plots for adjusted and unadjusted model outputs
estlambda(GWASout$t.value^2,plot=TRUE,method="median")

# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[,c("id","phenotype")]
GWAA(genodata=genotype, phenodata=phenoSub2, filename="GWAStutorialoutUnadj.txt", nSplits=10)
GWASoutUnadj<-read.table("GWAStutorialoutUnadj.txt", header=TRUE, colClasses=c("character", rep("numeric",4)))

estlambda(GWASoutUnadj$t.value^2,plot=TRUE,method="median")

# Create LDheatmap
CETP<-CETP[order(CETP$position),] # Order CETP SNPs by position
subgen<-genotype[,colnames(genotype)%in%CETP$SNP] # Subset snpMatrix object
subgen<-cbind(subgen, impCETPgeno)
subgen<-subgen[,order(match(colnames(subgen),CETP$SNP))]

ld <- ld(subgen,subgen, stats="R.squared") # Find LD map of CETP SNPs

require(LDheatmap)
require(rtracklayer)
ll<-LDheatmap(ld,CETP$position,flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination, and scatterplot
llplusgenes <- LDheatmap.addGenes(ll, chr="chr16", genome="hg19")
llGenesRecomb <- LDheatmap.addRecombRate(llplusgenes, chr="chr16", genome="hg19", view = "full")
llGenesRecombScatter<-LDheatmap.addScatterplot(llGenesRecomb,CETP$Neg_logP, ylab="-log10 p-value")
