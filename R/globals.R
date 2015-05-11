# Modify data.dir to indicate the location of the GWAStutorial files
# Intermediate data files will also be stored in this same location unless you set out.dir

# customize as needed
data.dir <- '/Volumes/genome/Research/GWAS'
out.dir <- data.dir                     # may want to write to a separate dir to avoid clutter

genotype.subset.fname <- sprintf("%s/subsetted_genotype.RData",out.dir)

gwas.fn = lapply(c(bed='bed',bim='bim',fam='fam',gds='gds'), function(n) sprintf("%s/GWAStutorial.%s", data.dir, n))
gwaa.fname <- sprintf("%s/GWAStutorialout.txt", out.dir)
gwaa.unadj.fname <- sprintf("%s/GWAStutorialoutUnadj.txt", out.dir)
impute.out.fname <- sprintf("%s/GWAStutorial_imputationOut.csv", out.dir)
protein.coding.snps.fname <- sprintf("%s/protein_coding_SNPs.csv", out.dir)
CETP.fname <- sprintf("%s/CETP_GWASout.csv", out.dir)

hapmap.url <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/hapmap_format/polymorphic/genotypes_chr%s_CEU_phase3.2_nr.b36_fwd.txt.gz"

