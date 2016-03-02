# Reading in FAMuSS data:
fms <- read.delim("http://stat-gen.org/book.e1/data/FMS_data.txt", header=T, sep="\t")

# Necessary code from Example 3.9:
attach(fms)
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)]
FMSgenoNum <- data.matrix(FMSgeno)
FMSgenoNum[is.na(FMSgenoNum)] <- 4
DistFmsGeno <- as.matrix(dist(FMSgenoNum))

# Exammple 3.10 (Multidimensional scaling (MDS) for identifying population substructure):
plot(cmdscale(DistFmsGeno),xlab="C1",ylab="C2")
abline(v=0,lty=2)	
abline(h=4,lty=2)
