#Code Snippet 7: Association analysis of imputed SNPs

# Carry out association testing for imputed SNPs using single.snp.tests()
rownames(clinical) <- clinical$FamID
imp <- single.snp.tests(phenotype=hdl,stratum=sex,data=clinical,snp.data=target,rules=rules)

# Obtain p values for imputed SNPs
pVals <- p.value(imp,df=1)
pVals <- pVals[!is.na(pVals)]
results <- data.frame(names(pVals), pVals)
colnames(results) <- c("SNP", "p.value")
rownames(results) <- NULL # 19352 SNPs were tested

#Write a file containing the results
write.csv(results, "GWAStutorial_imputationOut.csv", row.names=FALSE)

# Merge imputation testing results with support to obtain coordinates
# Read in file containing SNPs mapped to protein coding genes within 5kb
genes<-read.csv("protein_coding_SNPs.csv")
colnames(genes)<-c("SNP", "gene","chr", "position")
imputeOut<-merge(results, genes)

# Find the -log_10 of the p-values
imputeOut$Neg_logP<--log10(imputeOut$p.value)

# Order by p-value
imputeOut<-imputeOut[order(imputeOut$p.value),]
head(imputeOut)

# Subset for CETP SNPs

impCETP<-imputeOut[imputeOut$gene=="CETP",c("SNP", "chr", "position","p.value", "Neg_logP")]
impCETP$type<-"imputed"
impCETPgeno<-imputed[,colnames(imputed)%in%impCETP$SNP]
