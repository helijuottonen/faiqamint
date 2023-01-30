# Heli Juottonen 2023 heli.juottonen@alumni.helsinki.fi
# R scripts for bacterial 16S rRNA gene sequence data analysis of Faiqa Atique's PhD manuscript on aquaponics with mint
# https://github.com/helijuottonen/faiqamint

# For the article: xxx
# Atique F 12, Juottonen H 1, Kytöviita MM 1
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Institute of Bioeconomy, JAMK University of Applied Sciences, Finland

# Disclaimer: I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices.
# Use at your own risk!

# 4: Redundancy analysis for correlation of bacterial community and fish weight

# RUN FIRST:
# 0_basic_processing.R

# OTU table used in this script: faotu.r from 0_basic_processing.R
# OTUs with less than 5 reads overall have been removed

# loading the needed packages
#library(vegan)
library(tidyverse)
library(phyloseq)
library(cowplot)

# removing samples with no fish weight data from the metadata file
fameta_fw <- subset(fameta, FishFW > 0)

# making sure the otu table and metadata contain the same samples
# this is done based on row names so both object must have identical format of row names
faotu_fw = faotu.r[row.names(faotu.r) %in% row.names(fameta_fw),]
fameta_fw2 = fameta_fw[row.names(fameta_fw) %in% row.names(faotu_fw),]

#ordering both files based on row names to make sure the order is the same
faotu_fw <- faotu_fw[order(row.names(faotu_fw)),]
fameta_fw2 <- fameta_fw2[order(row.names(fameta_fw2)),]

# subsetting metadata based on treatments
fameta_fw_m <- subset(fameta_fw2, Treatments=="aquaponics")
fameta_fw_wom <- subset(fameta_fw2, Treatments == "RAS")
fameta_fw_m <-  droplevels(fameta_fw_m)
fameta_fw_wom <-  droplevels(fameta_fw_wom)

# selecting the same samples from OTU tables
faotu_fw_m = faotu_fw[row.names(faotu_fw) %in% row.names(fameta_fw_m),]
faotu_fw_wom = faotu_fw[row.names(faotu_fw) %in% row.names(fameta_fw_wom),]

# standardizing OTU data to relative abundances because reads <5 have been removed
# = the read numbers are not the same anymore between samples
# using Hellinger transformation, good for OTU data and for RDA
faotu_fw_m.h <- decostand(faotu_fw_m, "hell")
faotu_fw_wom.h <- decostand(faotu_fw_wom, "hell")

# redundancy analysis (RDA) with fish weight as the only environmental variable
faotu_fw_m.rda <- rda(faotu_fw_m.h ~ FishFW + Condition(Organ), data=fameta_fw_m)
faotu_fw_wom.rda <- rda(faotu_fw_wom.h ~ FishFW + Condition(Organ), data=fameta_fw_wom)

# check results
# variation explained by fish biomass = constrained proportion
faotu_fw_m.rda
# Call: rda(formula = faotu_fw_m.h ~ FishFW + Condition(Organ), data =
#             fameta_fw_m)
# 
# Inertia Proportion Rank
# Total         0.57112    1.00000     
# Conditional   0.15780    0.27629    2
# Constrained   0.01903    0.03332    1
# Unconstrained 0.39429    0.69039   23
# Inertia is variance

faotu_fw_wom.rda
# Call: rda(formula = faotu_fw_wom.h ~ FishFW + Condition(Organ), data =
#             fameta_fw_wom)
# 
# Inertia Proportion Rank
# Total         0.571752   1.000000     
# Conditional   0.176468   0.308644    2
# Constrained   0.009252   0.016181    1
# Unconstrained 0.386032   0.675175   23
# Inertia is variance 

# signifigance testing

# based on:
# https://fromthebottomoftheheap.net/slides/advanced-vegan-webinar-2020/advanced-vegan#1
# https://www.r-bloggers.com/2014/11/analysing-a-randomised-complete-block-design-with-vegan/

set.seed(1732)

h_m <- how(blocks = fameta_fw_m$Organ, nperm = 999)
h_wom <- how(blocks = fameta_fw_wom$Organ, nperm = 999)

# checking the significance of the RDA model

# aquaponics (with mint)
anova(faotu_fw_m.rda, permutations = h_m)

# Permutation test for rda under reduced model
# Blocks:  fameta_fw_m$Organ 
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = faotu_fw_m.h ~ FishFW + Organ, data = fameta_fw_m)
# Df Variance      F Pr(>F)
# Model     3  0.17683 3.4382  0.311
# Residual 23  0.39429

# RAS (without mint)
anova(faotu_fw_wom.rda, permutations= h_wom)

# Permutation test for rda under reduced model
# Blocks:  fameta_fw_wom$Organ 
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = faotu_fw_wom.h ~ FishFW + Condition(Organ), data = fameta_fw_wom)
# Df Variance      F Pr(>F)
# Model     1  0.00925 0.5512  0.927
# Residual 23  0.38603               

