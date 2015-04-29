#Code Snippet 3: Sample level filtering

# Setting thresholds
sampcall<-0.9
hetcutoff<-0.1

# Find sample statistics
snpsum.row <- row.summary(genotype)

# Calculating F-stat
MAF<-snpsum.col$MAF
callmatrix<-!is.na(genotype)
hetExp<-as.vector(callmatrix%*%(2*MAF*(1-MAF)))
hetObs<-snpsum.row$Heterozygosity*(ncol(genotype))*snpsum.row$Call.rate
hetF<-1-(hetObs/hetExp)

sampleuse<-!is.na(snpsum.row$Call.rate) & snpsum.row$Call.rate > sampcall & abs(hetF)<=hetcutoff

sampleuse[is.na(sampleuse)] <- FALSE
nrow(genotype)-sum(sampleuse) #0 subjects removed

# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical<-clinical[clinical$FamID%in%rownames(genotype),]

# Checking for Relatedness
biocLite("SNPRelate")
require(SNPRelate)

# Create gds file, required for SNPRelate functions
snpgdsBED2GDS("GWAStutorial.bed", "GWAStutorial.fam", "GWAStutorial.bim", "GWAStutorial.gds")
genofile <- openfn.gds("GWAStutorial.gds")

#Prune SNPs for IBD analysis
set.seed(1000)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold=0.2,sample.id=paste(rownames(genotype),1,sep="-"),
snp.id=colnames(genotype))
snpset.ibd <- unlist(snpSUB)
#74870 SNPs will be used in IBD analysis

# Find IBD coefficients using Method of Moments procedure
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE, sample.id=paste(rownames(genotype),1,sep="-"),
snp.id=snpset.ibd ,num.thread=1)
ibdcoeff <- snpgdsIBDSelection(ibd)
ibdcoeff$ID1<-sub("-1", "", ibdcoeff$ID1)
ibdcoeff$ID2<-sub("-1", "", ibdcoeff$ID2)

# Check if there are any candidates for relatedness
sum(ibdcoeff$kinship>=0.1)
ibdposi<-ibdcoeff
removeID<-NULL
while(sum(ibdposi$kinship>=0.1)>0){
    print(sum(ibdposi$kinship>=0.1))
    ibdposi<-ibdposi[ibdposi$kinship>=0.1,]
    idvec<-data.frame(table(c(ibdposi$ID1, ibdposi$ID2)))
    idvec<-idvec[order(idvec$Freq, decreasing=TRUE),]
    removeID<-c(removeID,as.character(idvec$Var1[1]))
    ibdposi<-ibdposi[!ibdposi$ID1%in%removeID & !ibdposi$ID2%in%removeID,]}

length(removeID) #0 Subjects removed due to correlation coefficient >=0.1

genotype <- genotype[!rownames(genotype)%in%removeID,]
clinical<-clinical[clinical$FamID%in%rownames(genotype),]
dim(genotype) #All 1401 subjects remain

# Checking for ancestry (MDS)

# Find IBS pairwise distance matrix
ibs <- snpgdsIBS(genofile,sample.id=paste(rownames(genotype),1,sep="-"),  snp.id=snpset.ibd, num.thread=1)

#nPerform MDS on IBS pairwise distance matrix
loc <- cmdscale(1 - ibs$ibs, k = 2)
MDS1 <- loc[, 1]
MDS2 <- loc[, 2]

# Plot the first two MDS components
plot(MDS2, MDS1, xlab="MDS Component 2", ylab="MDS Component 1", main="Multidemensional Scaling Analysis for Ancestry")

save(genotype, file="subsetted_genotype.RData")
load("subsetted_genotype.RData")
