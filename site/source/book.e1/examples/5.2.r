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

# Example 5.2 (Calculating posterior haplotype probabilities)
HaploEM$nreps[1:5]
HaploEM$indx.subj[1:8]
HaploEM$hap1code[1:8]
HaploEM$hap2code[1:8]
HaploEM$post[1:8]
HapProb <- HaploEM$hap.prob
HapProb
p1 <- 2*prod(HapProb[c(3,8)])
p2 <- 2*prod(HapProb[c(4,7)])
p1 / (p1+p2)
p2 / (p1+p2)
