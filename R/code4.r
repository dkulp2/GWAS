#Code Snippet 4: Generating principal components

# customize as needed
out.dir <- '/Volumes/genome/Research/GWAS' # same as input dir, but can be different
genotype.subset.fname <- sprintf("%s/subsetted_genotype.RData", out.dir)
pcs.fn <- sprintf("%s/GWAStutorial_pcs.RData", out.dir)

##################

# load data created in previous snippets
load(genotype.subset.fname)

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

# Create a file of the first 10 principal components
save(pcs, file=pcs.fn)
