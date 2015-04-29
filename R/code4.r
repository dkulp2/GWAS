#Code Snippet 4: Generating principal components

# Find and record first 10 principal components
num.princ.comp<-10
xxmat <- xxt(genotype, correct.for.missing = FALSE)
evv <- eigen(xxmat, symmetric = TRUE)
pcs <- evv$vectors[,1:num.princ.comp]
evals <- evv$values[1:num.princ.comp]
btr <- snp.pre.multiply(genotype, diag(1/sqrt(evals)) %*% t(pcs))
pcs <- snp.post.multiply(genotype, t(btr))
colnames(pcs) <- paste("pc",1:num.princ.comp,sep="")
rm(xxmat)

# Create a file of the first 10 principal components
save(pcs, file="GWAStutorial_pcs.RData")