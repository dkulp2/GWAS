#Code Snippet 5: Genotype imputation

# Imputation of non-typed HapMap SNPs
presSnps <- colnames(genotype)

# read in hapmap data for chromosomes 16
chrNum<-16

biocLite("chopsticks")
require(chopsticks)
readHapMap <- function(chrNum){
    hapmap <- read.HapMap.data(paste("ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/
  	hapmap_format/polymorphic/genotypes_chr",chrNum,"_CEU_phase3.2_nr.b36_fwd.txt.gz",sep=""))
    list(hapmap)
}
HapMap16 <- readHapMap(chrNum)
detach("package:chopsticks", unload=TRUE)

# Impute SNPs on given chromosome

# Subset for SNPs on given chromosome
presDat <- genoBim[which(genoBim$SNP%in%presSnps),]
presDatChr <- presDat[which(presDat$chr==chrNum),]

# Obtain hapmap data for given chromosome
hapMatrix <- new("SnpMatrix",HapMap16[[1]]$snp.data)
support <- HapMap16[[1]]$snp.support

# Subset hapmap data for our SNPs
# "missing" and "present" are snpMatrix objects needed for imputation rules
presDatChr$SNP <- as.character(presDatChr$SNP)
missing <- which(!colnames(hapMatrix)%in%presDatChr$SNP)
missing <- hapMatrix[,missing]
present <- which(colnames(hapMatrix)%in%presDatChr$SNP)
present <- hapMatrix[,present]

# Obtain positions of SNPs to be used for imputation rules
pos.pres <- support$Position[colnames(hapMatrix)%in%presDatChr$SNP]
pos.miss <- support$Position[!colnames(hapMatrix)%in%presDatChr$SNP]

# Calculate and store imputation rules using snp.imputation()
rules <- snp.imputation(present,missing,pos.pres,pos.miss)

# Obtain 'best guess' genotypes of imputed snps
targetSnps <- c(presDatChr$SNP)
target <- genotype[,targetSnps]
imputed <- impute.snps(rules, target, as.numeric=FALSE)
# 21342 SNPs were imputed

rm(hapMatrix)
rm(missing)
rm(present)
rm(HapMap16)
