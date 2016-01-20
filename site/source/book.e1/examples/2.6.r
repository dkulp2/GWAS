# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 2.6 (Linear regression)
attach(fms)
Geno <- actn3_r577x
Trait <- NDRM.CH
ModFull <- lm(Trait~Geno+Gender+Geno*Gender, na.action=na.exclude)
print(summary(ModFull))
ModReduced <- lm(Trait~Geno+Gender, na.action=na.exclude)
anova(ModReduced, ModFull)
NewDat <- data.frame(Geno=rep(c("CC","CT","TT"),2), 
	Gender=c(rep("Female",3), rep("Male",3)))
predict.lm(ModFull, NewDat, interval="prediction", level=0.95)
