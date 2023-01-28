# Heli Juottonen 2023 heli.juottonen@alumni.helsinki.fi
# R scripts for bacterial 16S rRNA gene sequence data analysis of Faiqa Atique's PhD manuscript on aquaponics with mint
#  https://github.com/helijuottonen/faiqamint

# For the article: xxx
# Atique F 12, Juottonen H 1, Kytöviita MM 1
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Institute of Bioeconomy, JAMK University of Applied Sciences, Finland

# Disclaimer: I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices.
# Use at your own risk!


# 0: Basic processing of OTU table and metadata:
# starting with shared file from mothur and metadata file


# 1. Start by setting up a project:
# Make a new folder for your project and give it a short informative name
# Inside the project folder, create a folder called 'data', and a folder called 'plots'.
# Move your data files (here shared, taxonomy, metadata) into the 'data' folder.
# Create a new R project: RStudio -> File -> New project -> Existing directory -> choose the project folder you made

# 2. Getting the OTU and metadata ready for further analyses

# loading necessary packages
library(vegan)
library(tidyverse)

# reading in OTU tables (shared file as it comes from mothur):
faotu <- read.table("data/fams2_all.unique.good.filter.unique.precluster.denovo.uchime.pick.opti_mcc.shared", header=T)

# reading in the metadata file:
# IMPORTANT: use 'read.csv' if the separator in the csv file is comma (,). Use 'read.csv2' is the separator is semicolon (;).
# note: the first column must have the sample names in the same format as they are in the shared file
# use this if the separator in the csv file is ;
fameta <- read.csv("data/ms2_metadata_090122.csv", header=T, row.names=1)

# converting otu table into dataframe
faotu.df <- as.data.frame(faotu)

# insert sample names in the OTU table as row names (needed in later commands):
row.names(faotu.df) <- faotu.df$Group

# remove unnecessary columns in the OTU table:
# done because the first three columns (=1:3) of the mothur shared files contain unnecessary values
# subset and -c means we are keeping else except the first 3 columns
faotu.df <- subset(faotu.df, select = -c(1:3))

# removing leaf samples and controls
fameta <- subset(fameta, Organ != "leave" & Organ != "cont" & Organ != "feed")

# Note, important: check what is in the negative controls at some point!

# making sure the OTU table and metadata contain the same samples
# this is done based on row names so otu table and metadata file must have identical formats of row names/sample names

fameta <- fameta[row.names(fameta) %in% row.names(faotu.df),]
faotu.df2 <- faotu.df[row.names(faotu.df) %in% row.names(fameta),]

# ordering otu table and metadata based on row names to make sure the order is the same
faotu.df2 <- faotu.df2[order(row.names(faotu.df2)),]
fameta <- fameta[order(row.names(fameta)),]

# adding read numbers of metadata object
fameta$reads <- rowSums(faotu.df2)

# subsampling the data (rarefying):

# defining a function for subsampling based on medians:
# picking the same number of reads from all samples: median of read numbers (=row sums)
# except samples that have less reads than the median: all reads are taken
# note: especially for comparing alpha diversity it would be better to have the exact same number of reads
# this is a compromise for situations when some samples have few reads
median_raref <- function(x) {
  x_sums <- rowSums(x)
  x_sums2 <- replace(x_sums, x_sums>(median(x_sums)), median(x_sums))
  set.seed(1732)
  x.r <- rrarefy(x, x_sums2)
  x.rdf <- as.data.frame(x.r)
  return(x.rdf)
}

median(rowSums(faotu.df2))
# 6454.5

# sumsampling to median number of reads:
faotu.r.all <- median_raref(faotu.df2)

# removing OTUs with less than 5 reads overall
faotu.r <- faotu.r.all[,colSums(faotu.r.all) > 5]

# converting to relative abundances to even out the effects of median rarefaction, removal of rare OTUs
faotu.st <- decostand(faotu.r, "normalize")

# OTU tables, metadata ready for further analyses 

