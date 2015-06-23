# ---- code0 ----
source("globals.R")

# ---- code0-a ----
# Run this once interactively to download and install BioConductor packages and other packages.

source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
biocLite("SNPRelate")
biocLite("rtracklayer")
biocLite("biomaRt")
install.packages(c('plyr', 'GenABEL', 'LDheatmap','doParallel', 'ggplot2', 'coin', 'igraph', 'devtools'))

library(devtools)
install_url("http://cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")


# ---- code0-b ----
# Download and unzip data needed for this tutorial

download.file(urlSupport, zipSupport.fn, method = "curl")
unzip(zipSupport.fn, exdir = data.dir)

