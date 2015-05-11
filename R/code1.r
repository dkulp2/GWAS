#Code Snippet 1: Reading data into R

source("globals.R")

##################

library(snpStats)                       # loading SNPs and performing QC

# Read in PLINK filesto create snpMatrix objects
geno <- read.plink(gwas.fn$bed,gwas.fn$bim,gwas.fn$fam, na.strings=("-9"))

# geno is now a list of three large objects as follows:
# geno$genotypes
#  A SnpMatrix with  1401 rows and  861473 columns
#  Row names:  10002 ... 11596              <<< samples
#  Col names:  rs10458597 ... rs5970564     <<< SNPs
#
# geno$fam                                  <<< sample details
#        pedigree member father mother sex affected
#  10002    10002      1      0      0   1        1
#  10004    10004      1      0      0   2        1
#  10005    10005      1      0      0   1        2
#  10007    10007      1      0      0   1        1
#  10008    10008      1      0      0   1        2
#  10009    10009      1      0      0   1        2
#  ...
#
# geno$map                                  <<< SNP details
#             chromosome   snp.name cM position allele.1 allele.2
#  rs10458597          1 rs10458597  0   564621     <NA>        C
#  rs12565286          1 rs12565286  0   721290        G        C
#  rs12082473          1 rs12082473  0   740857        T        C
#  rs3094315           1  rs3094315  0   752566        C        T
#  rs2286139           1  rs2286139  0   761732        C        T
#  rs11240776          1 rs11240776  0   765269        G        A
#  ...

# Obtain the genotypes table from snpMatrix objects
# Note: Phenotypes and covariates will be read from the clinical data file, below
genotype <- geno$genotype
print(dim(genotype))                  # 861473 SNPs read in for 1401 subjects

#Obtain the SNP information table from snpMatrix objects
genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
#chr16.snps <- genoBim[genoBim$chr == '16','SNP']

# Remove raw file to open up memory
rm(geno)

# Read in clinical file
clinical <- read.csv(sprintf("%s/GWAStutorial_clinical.csv", data.dir),
                     colClasses=c("character", "factor", "factor", "numeric", "numeric"))
rownames(clinical) <- clinical$FamID

# Subset genotype for subject data (and for chr16 -- DEBUG ONLY, REMOVE BEFORE PUB)
# genotype <- genotype[clinical$FamID, chr16.snps ]
# All 1401 subjects contain both clinical and genotype data

##################

# Write genotype, genoBim, clinical for future use
save(genotype, genoBim, clinical, file=genotype.subset.fname)
