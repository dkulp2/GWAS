# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 5.1 (EM approach to haplotype frequency estimation)
attach(fms)
install.packages("haplo.stats")
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
