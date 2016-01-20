# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 3.9 (Generating a similarity matrix):
attach(fms)
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
NamesAkt1Snps
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)]
FMSgenoNum <- data.matrix(FMSgeno)
FMSgenoNum[is.na(FMSgenoNum)] <- 4
DistFmsGeno <- as.matrix(dist(FMSgenoNum))
DistFmsGeno[1:5,1:5]
