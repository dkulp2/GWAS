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
HaploEM <- haplo.em(Geno.C, locus.label=SNPnames, control=haplo.em.control(min.posterior=1e-4))
HaploEM

Geno.AA <- Geno[Race=="African Am" & !is.na(Race),]
HaploEM2 <- haplo.em(Geno.AA, locus.label=SNPnames, control=haplo.em.control(min.posterior=1e-4))
HaploEM2

# New functions needed for Example 5.3

##########################################################################
# Description: This function creates a design matrix with i,j 
#		element equal to the conditional expectation 
#		of the number of copies of haplotype j for 
#		individual i based on the output from haplo.em()
# Input:	HaploEM (object resulting from haplo.em())
# Output:	XmatHap
##########################################################################
 
HapDesign <- function(HaploEM){
	Nobs <- length(unique(HaploEM$indx.subj)) # number of observations
	Nhap <- length(HaploEM$hap.prob)	# number of haplotypes
	XmatHap <- matrix(data=0,nrow=Nobs,ncol=Nhap)
	for (i in 1:Nobs){
		IDSeq <- seq(1:sum(HaploEM$nreps))[HaploEM$indx.subj==i]
		for (j in 1:length(IDSeq)){
			XmatHap[i,HaploEM$hap1code[IDSeq][j]] <- 
				XmatHap[i,HaploEM$hap1code[IDSeq][j]] + 
				HaploEM$post[IDSeq][j]
			XmatHap[i,HaploEM$hap2code[IDSeq][j]] <- 
				XmatHap[i,HaploEM$hap2code[IDSeq][j]] + 
				HaploEM$post[IDSeq][j]
			}
		}	
	return(XmatHap)
}

##########################################################################
# Description: This function creates a vector with jth element 
#		equal to the standard error of haplotype j 
#		based on the output from haplo.em()
# Input: 	HaploEM (object resulting from haplo.em())
# Output: 	HapSE
##########################################################################
 

HapFreqSE <- function(HaploEM){
	HapMat <- HapDesign(HaploEM)
	Nobs <- length(unique(HaploEM$indx.subj)) # number of observations
	Nhap <- length(HaploEM$hap.prob)	# number of haplotypes
	S.Full<-matrix(data=0, nrow=Nobs, ncol=Nhap-1)
	for(i in 1:Nobs){
		for(k in 1:(Nhap-1)){
		S.Full[i,k]<-HapMat[i,k]/HaploEM$hap.prob[k]-
			HapMat[i,Nhap]/HaploEM$hap.prob[Nhap]
		}
	}
	Score<-t(S.Full)%*%S.Full
	invScore<-solve(Score)
	HapSE<-c(sqrt(diag(invScore)), 
		sqrt(t(rep(1,Nhap-1))%*%invScore%*%rep(1,Nhap-1)))
	return(HapSE)
	}


# Example 5.3 (Testing hypotheses about haplotype frequencies within the EM framework)
FreqDiff <- HaploEM2$hap.prob[4] - HaploEM$hap.prob[4]
s1 <- HapFreqSE(HaploEM)[4] 
s2 <- HapFreqSE(HaploEM2)[4]
SE <- sqrt(s1^2 + s2^2)
CI <- c(FreqDiff - 1.96*SE, FreqDiff + 1.96*SE)
CI
