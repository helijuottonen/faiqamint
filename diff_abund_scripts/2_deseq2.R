
# DESeq2

# DEseq2 analysis to find OTUs occurring differentially in two treatments

# DEseq2 manual: https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf

# library(vegan)
# library(phyloseq)
# library(tidyverse)
# library(cowplot)
library("DESeq2")

DAdeseq <- function(filt_otu, filt_meta, tr_name) {

# creating deseq object
# here countData = transposed OTU table
# colData = metadata
# design = the column of metadata that contains the treatment
filtdata.ds <- DESeqDataSetFromMatrix(countData= filt_otu, colData=filt_meta, design= ~Treatments)

# correcting the order of treatments so that the control is first
# if not corrected, the first in alphabetical order will be 'control'
filtdata.ds$Treatments <- factor(filtdata.ds$Treatments, levels = c("RAS", "aquaponics"))
# run the deseq analysis
filtdata.ds <- DESeq(filtdata.ds)
# write deseq results into an object
filtdata.ds.res <- results(filtdata.ds)
# check a summary of the results
summary(filtdata.ds.res)
# write the results into a file
write.csv2(filtdata.ds.res, paste("results/deseq_", tr_name, ".csv", sep = ""))

}


