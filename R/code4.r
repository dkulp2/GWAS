#Code Snippet 4: Generating principal components

source("globals.R")

# load data created in previous snippets
load(genotype.subset.fname)

##################

library(snpStats)

# Find and record first 10 principal components
# pcs will be a N:10 matrix.  Each column is a principal component.
num.princ.comp <- 10
xxmat <- xxt(genotype, correct.for.missing = FALSE) # SHOULDN'T THIS BE APPLIED TO ONLY A SINGLE CHROMOSOME AT A TIME? (see xxt docs)
evv <- eigen(xxmat, symmetric = TRUE)
evecs <- evv$vectors[,1:num.princ.comp]
evals <- evv$values[1:num.princ.comp]
btr <- snp.pre.multiply(genotype, diag(1/sqrt(evals)) %*% t(evecs))
pcs <- snp.post.multiply(genotype, t(btr))
colnames(pcs) <- paste("pc",1:num.princ.comp,sep="")

rm(xxmat)

##################

# Store pcs for future reference with the rest of the derived data
save(genotype, genoBim, clinical, pcs, file=genotype.subset.fname)
