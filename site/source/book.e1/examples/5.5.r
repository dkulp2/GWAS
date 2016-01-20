# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Necessary code from Example 5.1 
attach(fms)
library(haplo.stats)
Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
 	substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
 	substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
 	substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
SNPnames <- c("actn3_r577x", "actn3_rs540874", "actn3_rs1815739",
	"actn3_1671064")
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]
HaploEM <- haplo.em(Geno.C,locus.label=SNPnames, control=haplo.em.control(min.posterior=1e-4))
HaploEM

# Necessary code from Example 5.4
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race)]

# Example 5.5 (Multiple imputation for haplotype effect estimation and testing)
Nobs <- sum(Race=="Caucasian", na.rm=T)
Nhap <- length(HaploEM$hap.prob)
D <- 1000
Est <- rep(0,D)
SE <- rep(0,D)
for (nimput in 1:D){
	Xmat <- matrix(data=0,nrow=Nobs,ncol=Nhap)
	for (i in 1:Nobs){
		IDSeq <- seq(1:sum(HaploEM$nreps))[HaploEM$indx.subj==i]
		if (length(IDSeq)>1){Samp <- sample(IDSeq,size=1,
			prob=HaploEM$post[IDSeq])}
		if (length(IDSeq)==1){Samp <- IDSeq}	
		Xmat[i,HaploEM$hap1code[Samp]] <-1
		Xmat[i,HaploEM$hap2code[Samp]] <-1
		}	
	h9 <- Xmat[,9]>=1
	Est[nimput] <- summary(lm(Trait~h9))$coefficients[2,1]
	SE[nimput] <- summary(lm(Trait~h9))$coefficients[2,2]
}
MeanEst <- mean(Est)
Wd <- mean(SE^2)
Bd <- (1/(D-1))*sum((Est-MeanEst)^2)
Td <- Wd + ((D+1)/D)*Bd
nu <- D-1*(1 + (1/(D+1))*(Wd/Bd))^2
1-pt(MeanEst/sqrt(Td),df=nu)
