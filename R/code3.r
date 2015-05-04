#Code Snippet 3: Sample level filtering

# customize as needed
out.dir <- '/Volumes/genome/Research/GWAS' # same as input dir, but can be different
genotype.subset.fname <- sprintf("%s/subsetted_genotype.RData",out.dir)

# Setting thresholds
sampcall <- 0.9
hetcutoff <- 0.1
ld.thresh <- 0.2
kin.thresh <- 0.1

##################

library(SNPRelate)                      # Estimating relatedness
library(plyr)

# load data created in previous snippets
load(genotype.subset.fname)

# Create sample statistics (Call rate, Heterozygosity)
snpsum.row <- row.summary(genotype)

# Calculating F-stat
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
hetF<-1-(hetObs/hetExp)

sampleuse <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall & abs(hetF)<=hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE    # remove NA's as well
cat(nrow(genotype)-sum(sampleuse), "subjects will be removed due to low sample call rate or high F statistic.\n") #0 subjects removed

# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical<- clinical[ rownames(genotype), ]

# Checking for Relatedness

# Create gds file, required for SNPRelate functions
snpgdsBED2GDS(gwas.fn$bed, gwas.fn$fam, gwas.fn$bim, gwas.fn$gds)
genofile <- snpgdsOpen(gwas.fn$gds)

#Prune SNPs for IBD analysis
set.seed(1000)
geno.sample.ids <- paste(rownames(genotype),1,sep="-") # Samples all have an extra "-1" in GDS.
snpSUB <- snpgdsLDpruning(genofile, ld.threshold=ld.thresh,
                          sample.id=geno.sample.ids, # Only analyze the filtered samples
                          snp.id=colnames(genotype)) # Only analyze the filtered SNPs
snpset.ibd <- unlist(snpSUB, use.names=FALSE)
length(snpset.ibd)  #74870 SNPs will be used in IBD analysis

# Find IBD coefficients using Method of Moments procedure.  Include pairwise kinship.
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
                    sample.id=geno.sample.ids,
                    snp.id=snpset.ibd,
                    num.thread=1)
ibdcoeff <- snpgdsIBDSelection(ibd)     # Pairwise sample comparison
ibdcoeff$ID1 <- sub("-1", "", ibdcoeff$ID1) # Remove the extra "-1" from the sample IDs
ibdcoeff$ID2 <- sub("-1", "", ibdcoeff$ID2) # ""

# Check if there are any candidates for relatedness
ibdcoeff <- ibdcoeff[ ibdcoeff$kinship >= kin.thresh, ]
nrow(ibdcoeff)

# iteratively remove samples with high kinship starting with the sample with the most pairings
related.samples <- NULL
while ( nrow(ibdcoeff) > 0 ) {

    # count the number of occurrences of each and take the top one
    sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
    rm.sample <- sample.counts[1, 'x']
    cat("Removing sample", as.character(rm.sample), 'too closely related to', sample.counts[1, 'freq'],'other samples.\n')

    # remove from ibdcoeff and add to list
    ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
    related.samples <- c(as.character(rm.sample), related.samples)
}

# filter genotype and clinical to include only unrelated samples
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- paste(rownames(genotype),1,sep="-") # Samples all have an extra "-1" in GDS.

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", kin.thresh,"\n") 
dim(genotype) #All 1401 subjects remain

# Checking for ancestry (MDS)

# Find IBS pairwise distance matrix
ibs <- snpgdsIBS(genofile,sample.id=geno.sample.ids,  snp.id=snpset.ibd, num.thread=1)

#nPerform MDS on IBS pairwise distance matrix
loc <- cmdscale(1 - ibs$ibs, k = 2)
MDS1 <- loc[, 1]
MDS2 <- loc[, 2]

# Plot the first two MDS components
plot(MDS2, MDS1, xlab="MDS Component 2", ylab="MDS Component 1", main="Multidemensional Scaling Analysis for Ancestry")

# Overwrite old genotype with new filtered version
save(genotype, genoBim, clinical, file=genotype.subset.fname)
