# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 4.2 (Tukeys single-step method):
attach(fms)
Trait <- NDRM.CH
summary(lm(Trait~resistin_c180g))
TukeyHSD(aov(Trait~resistin_c180g))

