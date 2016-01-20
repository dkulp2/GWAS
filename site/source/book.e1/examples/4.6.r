# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 4.6: Free step-down resampling adjustment
attach(fms)
Actn3Bin <- data.frame(actn3_r577x!="TT", actn3_rs540874!="AA",
				actn3_rs1815739!="TT", actn3_1671064!="GG")
Mod <- summary(lm(NDRM.CH~.,data=Actn3Bin))
Mod
TestStatObs <- Mod$coefficients[-1,3]
Tobs <- as.vector(sort(abs(TestStatObs)))
MissDat <- apply(is.na(Actn3Bin),1,any) | is.na(NDRM.CH)
Actn3BinC <- Actn3Bin[!MissDat,]
Ord <- order(abs(TestStatObs))
M <- 1000
NSnps <- 4
Nobs <- sum(!MissDat)
TestStatResamp <- matrix(nrow=M,ncol=NSnps)
for (i in 1:M){ 
		Ynew <- sample(Mod$residuals,size=Nobs,replace=T)
		ModResamp <- summary(lm(Ynew~.,data=Actn3BinC))
		TestStatResamp[i,] <- abs(ModResamp$coefficients[-1,3])[Ord]
		}
Qmat <- t(apply(TestStatResamp, 1, cummax))
Padj <- apply(t(matrix(rep(Tobs,M), NSnps)) < Qmat, 2, mean)
