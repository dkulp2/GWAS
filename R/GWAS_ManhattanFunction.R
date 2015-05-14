#Plot method
GWAS_Manhattan <- function(GWAS, col.snps=c("black","gray"), col.detected=c("blue"), col.imputed=c("red"), col.text="black", title="GWAS Tutorial Manhattan Plot", display.text=TRUE) {
    
    manhat <- GWAS[!grepl("[A-z]",GWAS$chr),]
    
    #sort the data by chromosome and then location
    manhat.ord <- manhat[order(as.numeric(manhat$chr),manhat$position),]
    manhat.ord <- manhat.ord[!is.na(manhat.ord$position),]
    
    ##Finding the maximum position for each chromosome
    max.pos <- sapply(1:21, function(i) { max(manhat.ord$position[manhat.ord$chr==i],0) })
    max.pos2 <- c(0, cumsum(max.pos))                  
    
    #Add spacing between chromosomes
    max.pos2 <- max.pos2 + c(0:21) * 10000000
    
    #defining the positions of each snp in the plot
    manhat.ord$pos <- manhat.ord$position + max.pos2[as.numeric(manhat.ord$chr)]
    
    #defining the coloring for the Manhattan plot
    manhat.ord$col[as.numeric(manhat.ord$chr)%%2==0] <- col.snps[1]
    manhat.ord$col[as.numeric(manhat.ord$chr)%%2==1] <- col.snps[2]
    
    text.pos <- sapply(c(1:22), function(i) { mean(manhat.ord$pos[manhat.ord$chr==i]) })
    
    # Plot the data
    plot(manhat.ord$pos[manhat.ord$type=="typed"]/1000000, manhat.ord$Neg_logP[manhat.ord$type=="typed"], pch=20, cex=.3, col= manhat.ord$col[manhat.ord$type=="typed"], xlab="Chromosome", ylab="Negative Log P-Value", axes=F, ylim=c(0,max(manhat$Neg_logP)+1))
    
    points(manhat.ord$pos[manhat.ord$type=="imputed"]/1000000, manhat.ord$Neg_logP[manhat.ord$type=="imputed"], pch=20, cex=.4, col = col.imputed)
    
    points(manhat.ord$pos[manhat.ord$type=="typed"]/1000000, manhat.ord$Neg_logP[manhat.ord$type=="typed"], pch=20, cex=.3, col = manhat.ord$col[manhat.ord$type=="typed"])
    
    
    axis(2)
    abline(h=0)
    
    SigNifSNPs <- as.character(GWAS[GWAS$Neg_logP > Lstringent.thresh & GWAS$type=="typed", "SNP"])
    
    #Add legend
    legend("topright",c("Bonferroni Corrected Threshold*", "Less Stringent Threshold**"),
           border="black", col=c("gray60", "gray60"), pch=c(0, 0), lwd=c(1,1),
           lty=c(1,2), pt.cex=c(0,0), bty="o", cex=0.7)
    
    #Add chromosome number
    text(text.pos/1000000, -.3, seq(1,22,by=1), xpd=TRUE, cex=1)
    
    #Add bonferroni line
    abline(h=bonferroni.thresh, untf = FALSE, col = "gray60")
    
    #Add "less stringent" line
    abline(h=Lstringent.thresh, untf = FALSE, col = "gray60", lty = 2 )
    
    #Plotting detected genes
    #Were any genes detected?
    if (length(SigNifSNPs)>0){

        sig.snps <- manhat.ord[,'SNP'] %in% SigNifSNPs
        
        points(manhat.ord$pos[sig.snps]/1000000,
               manhat.ord$Neg_logP[sig.snps],
             pch=15,col=col.detected, bg=col.detected,cex=0.5)
      
      text(manhat.ord$pos[sig.snps]/1000000,
           manhat.ord$Neg_logP[sig.snps],
           as.character(manhat.ord[sig.snps,1]), col=col.text, offset=1, adj=-.1, cex=.7)
    }
}
