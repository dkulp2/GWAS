# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",header=T,sep="\t")

# Example 7.2 (RF with missing SNP data - single imputation):
attach(fms)
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race) & !is.na(NDRM.CH)]
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)][Race=="Caucasian" & !is.na(Race) &!is.na(NDRM.CH),]
dim(FMSgeno)
round(apply(is.na(FMSgeno),2,sum)/dim(FMSgeno)[1],3)
library(randomForest)
FMSgenoRough <- na.roughfix(FMSgeno)
table(FMSgeno$"akt1_t22932c")
round(apply(is.na(FMSgenoRough),2,sum)/dim(FMSgeno)[1],3)
RandForRough <- randomForest(FMSgenoRough,Trait,importance=TRUE)
RandForRough$"importance"[order(RandForRough$"importance"[,1],decreasing=TRUE),]
