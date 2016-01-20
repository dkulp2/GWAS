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
	
# Example 4.3 (Benjamini and Hochberg Adjustment):
Pvec <- as.vector(Pvec)
m <- length(Pvec)
BHp <- sort(Pvec,decreasing=T)*m/seq(m,1)
sort(cummin(BHp))
BHp[order(Pvec,decreasing=T)] <- cummin(BHp)
names(PrMutSub)[BHp < 0.05]
sort(p.adjust(Pvec, method="BH"))



