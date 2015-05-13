#Code Snippet 4: Generating principal components

source("globals.R")

# load data created in previous snippets
load(genotype.subset.fname)

##################

library(SNPRelate)

# Find and record first 10 principal components
# pcs will be a N:10 matrix.  Each column is a principal component.
num.princ.comp <- 10

# Read in gds file for SNPRelate functions
genofile <- openfn.gds(gwas.fn$gds, readonly = TRUE)

#Prune SNPs for PCA

# Set LD threshold to 0.2
ld.thresh <- 0.2

set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, # Only analyze the filtered samples
                          snp.id = colnames(genotype)) # Only analyze the filtered SNPs
snpset.pca <- unlist(snpSUB, use.names=FALSE)
print(length(snpset.pca))  #72578 SNPs will be used in IBD analysis

# Find PCA
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.pca, num.thread=1)

# Create data frame of first ten principal comonents
pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : num.princ.comp],
                  stringsAsFactors = FALSE)
colnames(pcs)[2:11]<-paste("pc", 1:num.princ.comp, sep = "")

# Close GDS file
closefn.gds(genofile)

##################

# Store pcs for future reference with the rest of the derived data
save(genotype, genoBim, clinical, pcs, file=genotype.subset.fname)
