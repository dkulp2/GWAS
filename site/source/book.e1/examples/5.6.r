# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Necessary code from Example 5.1
attach(fms)
Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
 	substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
 	substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
 	substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]

# Example 5.6 (EM for estimation and testing of haplotype-trait association)
library(haplo.stats)
Geno.C <- setupGeno(Geno.C)
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race)]
Dat <- data.frame(Geno.C=Geno.C, Trait=Trait)
haplo.glm(Trait~Geno.C,data=Dat, allele.lev=attributes(Geno.C)$unique.alleles)
haplo.glm(Trait~Geno.C,data=Dat, allele.lev=attributes(Geno.C)$unique.alleles, 
	control=haplo.glm.control(haplo.base=9))
haplo.glm(Trait~Geno.C,data=Dat, allele.lev=attributes(Geno.C)$unique.alleles,
 	control=haplo.glm.control(haplo.effect="dominant"))

