# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",header=T,sep="\t")

# Necessary code from Example 7.2
attach(fms)
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race) & !is.na(NDRM.CH)]
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)][Race=="Caucasian" & !is.na(Race) &!is.na(NDRM.CH),]
library(randomForest)
FMSgenoRough <- na.roughfix(FMSgeno)

# Example 7.3 (RF with missing SNP data - multiple imputation):
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race) & !is.na(NDRM.CH)]
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms), NamesAkt1Snps)][Race=="Caucasian" & !is.na(Race) & !is.na(NDRM.CH),]
FMSgenoMI <- rfImpute(FMSgeno, Trait)
RandForFinal <- randomForest(FMSgenoMI[,-1], Trait, importance=TRUE)
RandForFinal$"importance"[order(RandForFinal$"importance"[,1], decreasing=TRUE),]
table(FMSgenoMI$akt1_t10726c_t12868c, FMSgenoRough$akt1_t10726c_t12868c)
