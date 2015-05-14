#Code Snippet 0: run this once interactively to download and install BioConductor and the snpStats package

source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
biocLite("SNPRelate")
biocLite("chopsticks")
biocLite("rtracklayer")
install.packages(c('plyr', 'GenABEL', 'LDheatmap', 'reshape'))
