GWAA <- function(genodata=genotypes,  phenodata=phenotypes, filename=NULL,
                 append=FALSE, nCores=4, flip=TRUE, select.snps=NULL,
                 type="parallel", nSplits=4)
{
    #Check to that a filename was specified
    if(is.null(filename)) stop("Must specify a filename for output.")
    
    #Check that the genotype data is of class 'SnpMatrix'
    if( class(genodata)!="SnpMatrix") stop("Genotype data must of class 'SnpMatrix'.")
    
    #Check that there is a variable named 'phenotype' in phenodata table
    if( !"phenotype" %in% colnames(phenodata))  stop("Phenotype data must have column named 'phenotype'")
    
    #Check that there is a variable named 'id' in phenodata table
    if( !"id" %in% colnames(phenodata)) stop("Phenotype data must have column named 'id'.")
    
    #If a vector of SNPs is given, subset genotype data for these SNPs
    if(!is.null(select.snps)) genodata<-genodata[,which(colnames(genodata)%in%select.snps)]
    
    #Check that there are still SNPs in 'SnpMatrix' object
    if(ncol(genodata)==0) stop("There are no SNPs in the 'SnpMatrix' object.")
    
    #Check if type of parallelization is specified
    if(!type %in% c("parallel","doParallel")) stop("Parallelization type not-supported")
    
    #Print the number of SNPs to be checked
    print(paste(ncol(genodata), " SNPs included in analysis."))
    
    #If append=FALSE than we will overwrite file with column names
    if(!isTRUE(append)) {
        columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
        write.table(t(columns), filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
    }

    #This function will change which allele is counted (major or minor)
    flip<-function(x) {
        zero2 <- which(x==0)
        two0 <- which(x==2)
        one1 <- which(x==1)
        mat <- matrix(data=NA, nrow=nrow(x), ncol=ncol(x))
        rownames(mat) <- rownames(x)
        colnames(mat) <- colnames(x)
        mat[zero2] <- 2
        mat[two0] <- 0
        mat[one1] <- 1
        mat
    }

    #This function contains the linear model for GWAS
    GWAAModel<-function(SNPuse, genodatasub=genoNum) {
        tempSNP <- data.frame(id=row.names(genodatasub), snp=genodatasub[,SNPuse])
        dat <- merge(phenodata, tempSNP, by.x="id", by.y="id", all.x=TRUE)
        a <- summary(glm(phenotype~ . - id, family = gaussian, data=dat))
        out <- as.matrix(a$coefficients['snp',])
        write.table(cbind(SNPuse,t(out)), filename, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
        invisible(out)
    }

    #We will split up the analysis into 4 subsets if there are at least 500 SNPs in analysis
    if(ncol(genodata) > 500) {
        nSNPs <- ncol(genodata)
        genosplit <- ceiling(nSNPs/nSplits)

        mat <- matrix(1:nSplits, nSplits, 3)
        mat[,2] <- genosplit*1:nSplits
        mat[2:nSplits,1] <- genosplit*1:(nSplits-1)+1
        mat[nSplits,2] <- nSNPs

        if(type=="parallel"){

            splitGWAS <- function(matrow) {
                genoNum <- as(genodata[,matrow[1]:matrow[2]], "numeric")
                #By default we will flip the numeric values of genotypes to count minor allele

                if (isTRUE(flip)) genoNum<-flip(genoNum)

                rsList <- as.list(as.matrix(colnames(genoNum)))
                mclapply(rsList, GWAAModel, genodata=genoNum ,mc.cores=nCores)
                print(paste("The GWAS is", 100*matrow[3]/nSplits, "% finished"))
            }
      
            apply(mat,1,splitGWAS)
        }
        else { # doParallel

            splitGWAS <- function(matrow, GWAAModel, phenodata, genodata) {
                cl <- makeCluster(nCores)
                registerDoParallel(cl)
                genoNum <- as(genodata[,matrow[1]:matrow[2]], "numeric")

                #By default we will flip the numeric values of genotypes to count minor allele
                if (isTRUE(flip)) genoNum<-flip(genoNum)

                rsVec <- colnames(genoNum)
                foreach(i=rsVec) %dopar% GWAAModel(i, genoNum)
                print(paste("The GWAS is", 100*matrow[3]/nSplits, "% finished"))
                stopCluster(cl)
            }
      
            apply(mat,1,splitGWAS, GWAAModel, phenodata, genodata)
        }
    } 
    else { # small SNP count
        genoNum <- as(genodata, "numeric")
        if (isTRUE(flip)) genoNum<-flip(genoNum)

        rsList <- as.list(as.matrix(colnames(genoNum)))

        lapply(rsList,GWAAModel, genodata=genoNum)
    }
  
  
    return(print("Done."))
}

