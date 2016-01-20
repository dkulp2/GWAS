# Reading in FAMuSS data
fmsURL <- "http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms <- read.delim(file=fmsURL, header=T, sep="\t")
