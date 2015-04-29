#Code Snippet 1: Reading data into R

# Download and open snpStats package
# Note: the library() function must be used at the beginning of each new R session
source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
library(snpStats)

# Read in PLINK filesto create snpMatrix objects
geno<-read.plink("GWAStutorial.bed","GWAStutorial.bim","GWAStutorial.fam", na.strings=("-9"))

# Obtain the genotypes table from snpMatrix objects
# Note: Phenotypes and covariates are being read from the clinical data file (see below)
genotype<-geno$genotype
dim(genotype)
#861473 SNPs read in for 1401 subjects

#Obtain the SNP information table from snpMatrix objects
genoBim<-geno$map
colnames(genoBim)<-c("chr", "SNP", "gen.dist", "position", "A1", "A2")

# Remove raw file to open up memory
rm(geno)

# Read in clinical file
clinical<-read.csv("GWAStutorial_clinical.csv", colClasses=c("character", rep("factor",2), rep("numeric", 2)))

# Subset genotype for subject data
genotype<-genotype[rownames(genotype)%in%clinical$FamID,]
# All 1401 subjects contain both clinical and genotype data
