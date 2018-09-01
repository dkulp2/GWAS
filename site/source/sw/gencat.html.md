# GenCAT - An R Package for Class-level Genetic Association Testing

### Overview

Methods for genetic association, such as genome-wide association studies (GWAS), are wideley used to test for SNP level relationships for traits and disease. Accounting for multiple-comparisons of such techniques typically calls for conservative Bonferroni adjustements for the number of SNPs included in the analysis (often ~ 1 million). Commonly referred to as the Minimum P approach, a *class level* signal is often given by the minimum p-value of all SNPs within a class, where class is defined as as a meaningful segment of the genome (such as a protein-coding gene). Therefore, only classes with a single strong SNP level signal are deemed significant. In this package we present tools for implementing an alternative method, known as **Genetic Class Association Testing (GenCAT)**, which utilizes statistical signal from every SNP accross an entire class for each test. This method (1) reduces the number of independent tests and thus provides less demanding thresholds and (2) detects consistent moderate signal within a class. It requires defining classes of interest, mapping SNPs to said classes, and utilizing the unique within-class correlation structure to compute class-level statistics. The GenCAT package provides R functions for mapping SNPs to pre-defined classes, running the GenCAT analysis, and visualizing the results with a customized manhattan plot.

GenCAT package details can be found on the [CRAN website](https://cran.r-project.org/web/packages/GenCAT/index.html).

