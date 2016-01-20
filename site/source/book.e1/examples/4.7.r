# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Neccessary code from Example 4.6:
attach(fms)
Actn3Bin <- data.frame(actn3_r577x!="TT",actn3_rs540874!="AA",
				actn3_rs1815739!="TT",actn3_1671064!="GG")
Mod <- summary(lm(NDRM.CH~.,data=Actn3Bin))
MissDat <- apply(is.na(Actn3Bin),1,sum)>0 | is.na(NDRM.CH)
Actn3BinC <- Actn3Bin[!MissDat,]
Nobs <- sum(!MissDat)

# Example 4.7 (Null unrestricted bootstrap approach):
CoefObs <- as.vector(Mod$coefficients[-1,1])
B <- 1000
NSnps <- 4
Nobs <- sum(!MissDat)
TestStatBoot <- matrix(nrow=B,ncol=NSnps)
for (i in 1:B){ 
		SampID <- sample(1:Nobs,size=Nobs, replace=T)
		Ynew <- NDRM.CH[!MissDat][SampID]
		Xnew <- Actn3BinC[SampID,]
		CoefBoot <- summary(lm(Ynew~.,data=Xnew))$coefficients[-1,1]
		SEBoot <- summary(lm(Ynew~.,data=Xnew))$coefficients[-1,2]
		if (length(CoefBoot)==length(CoefObs)){
			TestStatBoot[i,] <- (CoefBoot-CoefObs)/SEBoot
			}
		}
for (cj in seq(2.7,2.8,.01)){
	print(cj)
	print(mean(apply(abs(TestStatBoot)>cj,1,sum)>=1,na.rm=T))
	}
	

