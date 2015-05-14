# ---- code0 ----
# Run this once interactively to download and install BioConductor packages and other packages.

source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
biocLite("SNPRelate")
biocLite("chopsticks")
biocLite("rtracklayer")
install.packages(c('plyr', 'GenABEL', 'LDheatmap', 'reshape'))
