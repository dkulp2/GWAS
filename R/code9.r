#Code Snippet 9: Visualizing and QC of GWA findings

source("globals.R")

# load derived data from previous snippets
load(working.data.fname)

##################

library(GenABEL)
library(LDheatmap)
library(rtracklayer)
library(plyr)
library(doParallel)
library(ggplot2)

# Create Manhattan Plot
source("GWAS_ManhattanFunction.R")
par(mfrow=c(1,1))
GWAS_Manhattan(GWAScomb)

# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[,c("id","phenotype")]

source("GWAA.R")
GWAA(genodata=genotype, phenodata=phenoSub2, filename=gwaa.unadj.fname)
GWASoutUnadj <- read.table(gwaa.unadj.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Create QQ plots for adjusted and unadjusted model outputs
estlambda(GWASout$t.value^2,plot=TRUE,method="median")
estlambda(GWASoutUnadj$t.value^2,plot=TRUE,method="median")

# Add "rs247617" to CETP
CETP <- rbind.fill(GWASout[GWASout$SNP == "rs247617",], CETP)

# Combine genotypes and imputed genotypes for CETP region
subgen <- genotype[, colnames(genotype) %in% CETP$SNP ]   # CETP subset of snpMatrix from typed SNPs
subgen <- cbind(subgen, impCETPgeno)                      # CETP subset of snpMatrix from imputed SNPs

# Subset SNPs for only certain genotypes
certain <- NULL
for(i in 1:ncol(subgen)){
  certain[i] <- sum(!unique(as(subgen[,i], "numeric")) %in% c(NA, 0, 1, 2)) == 0}
subgen <- subgen[,certain]

# Subset and order CETP SNPs by position
CETP <- CETP[CETP$SNP %in% colnames(subgen),]
CETP <- arrange(CETP, position)
subgen <- subgen[, order(match(colnames(subgen),CETP$SNP)) ]

# Create LDheatmap
ld <- ld(subgen, subgen, stats="R.squared") # Find LD map of CETP SNPs

ll <- LDheatmap(ld, CETP$position, flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination
llplusgenes <- LDheatmap.addGenes(ll, chr = "chr16", genome = "hg19", genesLocation = 0.01)

# Add scatterplot

scatter2<-ggplotGrob(
  {
    qplot(position, Neg_logP, data = CETP, xlab="", ylab = "Negative Log P-value", xlim = range(CETP$position), asp = 1/10, color = factor(type), colour=c("#000000", "#D55E00")) + 
      theme(axis.text.x = element_blank(),
            axis.title.y = element_text(size = rel(0.75)), legend.position = "none", 
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) +
      scale_color_manual(values = c("red", "black"))
  }
)

plot.new()
llQplot2<-LDheatmap.addGrob(llplusgenes, rectGrob(gp = gpar(col = "white")),height = .2)
pushViewport(viewport(x = 0.485, y= 0.76, width = 1,height = .4))
grid.draw(scatter2)