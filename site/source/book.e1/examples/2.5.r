# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 2.5 (M-sample tests of association for a quantitative trait)
attach(fms)
Geno <- as.factor(resistin_c180g)
Trait <- NDRM.CH
AnovaMod <- lm(Trait~Geno, na.action=na.exclude)
summary(AnovaMod)
kruskal.test(Trait, Geno, na.action=na.exclude)
