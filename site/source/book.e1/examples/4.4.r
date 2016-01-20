# Reading in Virco data:
virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv",header=T,sep=",")

# Necessary code from Example 4.1:
attach(virco)
PrMut <- virco[,23:121]!="-" & virco[,23:121]!="."
NObs <- dim(virco)[1]
PrMutSub <-data.frame(PrMut[ , apply(PrMut,2,sum) > NObs*.05])
Trait <- IDV.Fold - NFV.Fold
TtestP <- function(Geno){
	return(t.test(Trait[Geno==1],
		Trait[Geno==0], na.rm=T)$"p.value")
	}
Pvec <- apply(PrMutSub, 2, TtestP)
Pvec <- as.vector(Pvec)

# Example 4.4 (Benjamini and Yekutieli adjustment)
BYp <- p.adjust(Pvec, method="BY")
sort(BYp)
names(PrMutSub)[BYp < 0.05]


